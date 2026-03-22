# 구현 계획: Fog 강화 + Superposition 수정 + 비-braille 경로 제거

---

## 문제 1: Depth Fog가 너무 약함

### 원인 분석 (복합적 — 단일 원인 아님)

#### A. 색상 대비 부족

현재 near/mid/far 3단계의 실제 RGB 차이:

**Protein 0 (olive) 예시**:
| Band | xterm | RGB | 명칭 |
|------|-------|-----|------|
| near | 142 (PROTEIN_NEAR) | (175, 175, 0) | 밝은 올리브 |
| mid | 100 (PROTEIN_COLORS) | (135, 135, 0) | 기존 올리브 |
| far | 58 (PROTEIN_DIM) | (95, 95, 0) | 어두운 올리브 |

**near↔mid RGB 차이: (40, 40, 0) — 너무 작음!**
**mid↔far RGB 차이: (40, 40, 0) — 역시 너무 작음!**

xterm-256 팔레트는 40단위(0, 95, 135, 175, 215, 255) 간격이라, 한 단계 차이는 인지하기 어려움.

#### B. 3단계가 너무 적음

near/mid/far 3단계로는 그라데이션이 아니라 거의 균일하게 보임. 특히 대부분의 원자가 mid band에 몰릴 가능성이 높음.

#### C. depth_base 범위가 focal_offset에 의해 압축될 수 있음

`focal_offset = clamp(2.5 * radius + 1.0, 2.0, 8.0)` — 정규화 후 구조물 반경이 ~1.0이면 focal_offset ≈ 3.5. 이 경우:
- 모든 z값이 `z + 3.5` → 범위 ≈ [2.5, 4.5]
- depth_range = 2.0이지만, 구조물의 실제 z 분포가 좁으면 대부분 같은 band에 진입

#### D. calibrate_depth_baseline이 1회성이며 화면 내 원자만 측정

`calibrate_depth_baseline_first_view()`는:
1. 초기 뷰에서 **화면 안에 투영되는 원자만** z범위 측정 (Screen.cpp:372-374)
2. 카메라 회전 시 **재보정하지 않음** (`depth_calibrated = true` 후 갱신 없음)
3. `minSpan = 0.3f`로 최소 범위 강제 — 실제 z분포가 좁아도 0.3으로 확장

**결과**: 정면 뷰에서 z범위가 0.2 미만이면 minSpan 0.3으로 확장되고, 90° 회전 후에는 z범위가 1.7까지 넓어지지만 재보정되지 않아 depth band가 완전히 무의미해짐.

#### E. Turquoise, Coral 등 고휘도 색상은 near가 거의 변하지 않음

| 색상 | mid 휘도 | near 휘도 | near/mid 비율 |
|------|---------|----------|--------------|
| Turquoise | 0.702 | 0.749 | **1.07x** — 사실상 동일 |
| Coral | 0.632 | 0.692 | **1.10x** |
| Average | 0.498 | 0.663 | **1.33x** — 전문 기준(2-4x)의 절반 이하 |

### 해결 방안

#### 1-A: 더 넓은 색상 대비

xterm-256에서 **2단계 이상** 차이나는 색상 선택:

**새 PROTEIN_NEAR_COLORS (더 밝게)**:
```
olive:     184 #d7d700 (215,215,0)   ← 기존 142 (175,175,0) 대비 +40
turquoise: 159 #afffff (175,255,255) ← 기존 116 (135,215,215) 대비 +40,+40,+40
navy:      105 #8787ff (135,135,255) ← 기존 63 (95,95,255) 대비 +40,+40
purple:    177 #d787ff (215,135,255) ← 기존 135 (175,95,255) 대비 +40,+40
pink:      225 #ffd7ff (255,215,255) ← 기존 219 (255,175,255) 대비 +40
coral:     223 #ffd7af (255,215,175) ← 기존 216 (255,175,135) 대비 +40,+40
brown:     220 #ffd700 (255,215,0)   ← 기존 178 (215,175,0) 대비 +40,+40
orange:    229 #ffffaf (255,255,175) ← 기존 221 (255,215,95) 대비 +40,+80
red:       210 #ff8787 (255,135,135) ← 기존 203 (255,95,95) 대비 +40,+40
```

**새 PROTEIN_FAR_COLORS (기존 DIM보다 더 어둡게 — 회색 쪽으로)**:
기존 PROTEIN_DIM_COLORS 대신 xterm grayscale 톤 사용:
```
olive far:     242 #6c6c6c (108,108,108) ← 색조 제거, 순수 회색
turquoise far: 243 #767676 (118,118,118)
navy far:      237 #3a3a3a (58,58,58)
purple far:    239 #4e4e4e (78,78,78)
pink far:      245 #8a8a8a (138,138,138)
coral far:     243 #767676 (118,118,118)
brown far:     240 #585858 (88,88,88)
orange far:    243 #767676 (118,118,118)
red far:       240 #585858 (88,88,88)
```

**near↔far RGB 차이 (olive 예시): (215,215,0) vs (108,108,108) = (107,107,108) — 매우 뚜렷!**
**효과: 가까운 원자는 선명한 색상, 먼 원자는 흐릿한 회색으로 fade out.**

#### 1-B: Chain, Rainbow 모드도 동일 원칙 적용

- Chain near: 더 밝은 변형 (기존보다 +40~80 RGB)
- Chain far: 회색조로 desaturate
- Rainbow near: 더 밝은 변형
- Rainbow far: 회색조 또는 더 어두운 변형

#### 1-C: pLDDT, conservation, interface, aligned 모드에도 fog 적용

현재 이 모드들은 fog가 적용되지 않음. depth_band 기반 색상 분기를 추가:
- 각 모드의 기존 색상을 mid band로
- near/far 변형 추가

#### 1-D: depth_band 경계 조정

현재 0.33/0.66 균등 분할 → 비선형 분할 고려:
```
near: t < 0.25  (가까운 25%)
mid:  t < 0.60  (중간 35%)
far:  t >= 0.60 (먼 40%)
```
이렇게 하면 가까운 물체의 bright 효과가 더 두드러짐.

#### 1-E: 회전 시 depth 재보정

현재 `depth_calibrated = true` 후 영구 고정. 매 프레임 또는 회전 시 `depth_calibrated = false`로 리셋하여 depth_base_min_z/max_z를 갱신:

```cpp
// 회전 입력 처리부 (do_rotation 등 호출 후):
depth_calibrated = false;
```

이렇게 하면 시점 변화에 따라 depth band가 항상 현재 뷰에 최적화됨.

### 변경 파일

| 파일 | 변경 |
|------|------|
| `Palette.hpp` | PROTEIN_NEAR_COLORS, PROTEIN_FAR_COLORS(새), CHAIN_NEAR/FAR, RAINBOW_NEAR/FAR 색상값 교체 |
| `Screen.cpp` | init_color_pairs(): 새 far 색상 pair 등록. assign_colors_to_points(): pLDDT/conservation/interface/aligned fog 추가. depth_band 경계 조정. 회전 시 depth 재보정 |
| `Camera.cpp` | 새 far pair 범위 xterm 매핑 |

---

## 문제 2: Superposition(구조 겹침)이 작동하지 않음

### 원인 분석

**핵심 버그: pan_x/pan_y가 grid offset으로 설정된 후 초기화되지 않음**

흐름:
1. `structty.cpp:54` — `normalize_proteins(utmatrix)` 호출
2. `utmatrix`가 비어있으면 `hasUT=false`
3. `Screen.cpp:312-322` — `n > 1 && !hasUT` → **grid layout**: `pan_x[i] = col_offset`, `pan_y[i] = row_offset`
4. 이후 `structty.cpp:104` — `apply_hit_transform(ti, hit)` 호출
5. `Screen.cpp:1372-1383` — `apply_foldseek_transform`: 좌표 변환 적용 + `yesUT = true`
6. **하지만 pan_x/pan_y는 그대로** grid offset 유지!

**결과**: 두 번째 protein의 좌표는 Kabsch로 회전/이동되었지만, 렌더링 시 `pan_x[1] ≠ 0`이므로 첫 번째 protein과 다른 위치에 그려짐. **Superposition은 계산되었으나 화면에서 겹치지 않음**.

동일 문제가 FoldMason에서도 발생:
1. `structty.cpp:208` — `apply_foldmason_superposition(0, 1, ...)` 호출
2. 내부에서 Kabsch + `apply_foldseek_transform()` → `yesUT = true`
3. **pan_x/pan_y는 여전히 grid offset**

### 해결 방안

#### 2-A: apply_foldseek_transform에서 pan 초기화

transform 적용 후 모든 protein의 pan을 0으로 리셋:

```cpp
void Screen::apply_foldseek_transform(...) {
    // ... 기존 transform 적용 ...
    yesUT = true;

    // pan 초기화 — overlay 모드로 전환
    for (size_t i = 0; i < pan_x.size(); i++) {
        pan_x[i] = 0.0f;
        pan_y[i] = 0.0f;
    }
}
```

#### 2-B: focal_offset 재계산

transform 적용 후 구조물의 bounding box가 변경됨. `depth_calibrated = false`는 이미 설정되지만, `focal_offset`도 재계산 필요:

```cpp
// apply_foldseek_transform 또는 apply_foldmason_superposition 끝에:
float radius = compute_scene_radius_from_render_positions(data);
focal_offset = std::clamp(2.5f * radius + 1.0f, 2.0f, 8.0f);
depth_calibrated = false;
```

#### 2-C: normalize_proteins 내 순서 재구성 (대안)

또는 `normalize_proteins`에서 `-fs`/`-fm` 플래그도 `hasUT`처럼 취급:
- 하지만 이는 normalize_proteins가 foldseek/foldmason 상태를 알아야 하므로 결합도가 높아짐
- **방안 2-A가 더 깔끔함** — transform 적용 시점에서 pan을 리셋

### 변경 파일

| 파일 | 변경 |
|------|------|
| `Screen.cpp` | `apply_foldseek_transform()`: pan_x/pan_y 전체 0으로 리셋, focal_offset 재계산 |

---

## 문제 3: 비-braille 경로 제거

### 현재 상태

- `use_braille`은 `Screen.hpp:116`에서 `true`로 하드코딩. CLI 토글 없음.
- 비-braille 경로는 실질적으로 사용되지 않으나 코드가 남아있음.

### 영향 범위

#### 제거 대상 코드

| 위치 | 내용 | 비고 |
|------|------|------|
| `Screen.hpp` | `screenPixels` 멤버, `use_braille` 멤버 | 삭제 |
| `Screen.cpp::project()` | 비-braille 경로 (lines 695-808 부근) | 삭제 |
| `Screen.cpp::print_screen()` | 비-braille 분기 (lines 1072-1095 부근) | 삭제 |
| `Screen.cpp::update_hover_info()` | 비-braille 분기 | 삭제 |
| `Screen.cpp::clear_screen()` | screenPixels 초기화 | 삭제 |
| `Screen.cpp::assign_colors_to_points()` | `use_braille` 조건 분기 (fog) | 조건 제거, braille 경로만 남김 |
| `Screen.cpp::draw_line()` | `use_braille` 조건 분기 (depth_band) | 조건 제거, 항상 계산 |
| `Parameters.hpp/cpp` | `depthcharacter`, `-d` 플래그 | 삭제 |

#### 스크린샷(Camera.cpp) 리팩토링

현재 스크린샷은 `project(screenshotPixels, w, h)` 오버로드를 사용하여 **비-braille 해상도**로 렌더링:

```cpp
// Screen.cpp, screenshot 호출부:
screenshotPixels.assign(screen_width * screen_height, RenderPoint());
project(screenshotPixels, screen_width, screen_height);
camera->screenshot(screenshotPixels);
```

비-braille 경로 제거 후, 스크린샷도 braille 해상도(2x × 4x)로 렌더링하도록 변경:

**방안 3-A: logicalPixels 기반 스크린샷**
- `Camera::screenshot()`이 `logicalPixels` (2W×4H 해상도) 수신
- Camera 내부에서 해상도 조정
- **장점**: 더 높은 해상도 스크린샷 생성
- **단점**: Camera.cpp 수정 필요

**방안 3-B: 스크린샷 전용 braille project**
- `project()` 오버로드를 제거하는 대신, 메인 `project()`의 logicalPixels를 스크린샷에도 사용
- `camera_width = screen_width * 2`, `camera_height = screen_height * 4`로 설정
- **채택**: 이 방안이 가장 깔끔

#### get_pixel_char_from_depth() 처리

이 함수는 braille 경로에서도 `RenderPoint.pixel`에 값을 할당하지만, braille 렌더링에서는 사용하지 않음. 그러나 삭제하면 `draw_line()`과 project() 내 호출부도 변경 필요. **dead code로 남겨두거나 삭제 모두 가능** — `draw_line()`에서 pixel 값을 설정하는 부분만 제거하면 깔끔.

### 변경 파일

| 파일 | 변경 |
|------|------|
| `Screen.hpp` | `screenPixels` 삭제, `use_braille` 삭제 |
| `Screen.cpp` | 비-braille project/print/hover/clear 경로 삭제, assign_colors/draw_line의 use_braille 조건 제거 |
| `Camera.hpp/cpp` | braille 해상도(2W×4H) 기반 스크린샷으로 변경 |
| `Parameters.hpp/cpp` | depthcharacter, -d 플래그 삭제 |

---

## 구현 순서

### Step 1: Superposition pan 버그 수정 (문제 2)
- `apply_foldseek_transform()`에서 pan_x/pan_y 전체 0 리셋
- focal_offset 재계산 추가
- **가장 간단하고 즉각적인 효과**

### Step 2: Fog 색상 강화 (문제 1-A, 1-B)
- `Palette.hpp`: PROTEIN_NEAR/FAR, CHAIN_NEAR/FAR, RAINBOW_NEAR/FAR 색상값 교체
- near = 매우 밝은 색, far = 회색조 desaturate
- `Screen.cpp`: init_color_pairs 업데이트

### Step 3: Fog 모드 확장 (문제 1-C)
- pLDDT, conservation, interface, aligned 모드에 depth_band 기반 fog 추가
- 각 모드별 near/far 색상 pair 추가

### Step 4: (선택적) depth_band 비선형 경계 (문제 1-D)
- 0.33/0.66 → 0.25/0.60 등으로 조정
- 가까운 물체의 bright 효과 강화

### Step 5: 비-braille 경로 제거 (문제 3)
- screenPixels 및 관련 코드 삭제
- Camera.cpp braille 해상도 기반으로 리팩토링
- use_braille 변수 제거
- Parameters에서 -d 플래그 및 depthcharacter 삭제

### Step 6: 빌드 및 검증
- 컴파일 에러 확인
- `-fs`, `-fm` 옵션으로 superposition 확인
- braille fog 원근감 확인
- 스크린샷 정상 출력 확인

---

## 수정 대상 파일 요약

| 파일 | 변경 내용 |
|------|----------|
| `src/visualization/Palette.hpp` | near/far 색상값 교체 (더 넓은 대비) |
| `src/visualization/Screen.hpp` | screenPixels, use_braille 삭제 |
| `src/visualization/Screen.cpp` | pan 리셋 버그 수정, fog 강화, 비-braille 경로 삭제 |
| `src/visualization/Camera.hpp` | braille 해상도 기반 스크린샷 |
| `src/visualization/Camera.cpp` | 동상 |
| `src/structure/Parameters.hpp` | depthcharacter 삭제 |
| `src/structure/Parameters.cpp` | -d 플래그 삭제 |
| `src/visualization/RenderPoint.hpp` | (변경 없음) |
