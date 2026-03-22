# 구현 계획: Aligned 모드 수정 + 점자 원근감 개선

## 작업 A: `-m aligned` 색상 대비 복원

### 문제 요약

`assign_colors_to_points()`에서 aligned/non-aligned 구분을 PROTEIN_COLORS(vivid) / PROTEIN_DIM_COLORS(dim)로 변경한 결과, 두 색상 간 RGB 차이가 너무 작아(특히 olive: 40,40,0) 시각적으로 구분 불가.

### 해결: 방안 D — Protein 색조 유지 + Non-aligned 통일 회색

**핵심 아이디어**: Aligned 부분은 protein 고유색의 **눈에 띄게 밝은 버전**, non-aligned 부분은 **통일된 어두운 회색**으로 처리. Protein 색은 aligned 영역에서만 보이므로 구분이 극명.

### 변경 사항

#### 1. Palette.hpp — PROTEIN_BRIGHT_COLORS 추가

```cpp
// Aligned 강조용 밝은 protein colors, pairs 101-109
inline const std::array<int, 9> PROTEIN_BRIGHT_COLORS = {
    148,  // 101  olive bright    #afd700  (기존 100 #878700 대비 +120 green)
    123,  // 102  turquoise bright #87d7ff (기존 80 #5fd7d7 대비 blue 강화)
     69,  // 103  navy bright     #5f87ff  (기존 27 #005fff 대비 밝기 상승)
    171,  // 104  purple bright   #d75fff  (기존 129 #af00ff 대비 밝기 상승)
    219,  // 105  pink bright     #ffafff  (기존 213 #ff87ff 대비 밝기 상승)
    216,  // 106  coral bright    #ffaf87  (기존 209 #ff875f 대비 밝기 상승)
    178,  // 107  brown bright    #d7af00  (기존 130 #af5f00 대비 밝기 상승)
    221,  // 108  orange bright   #ffd75f  (기존 214 #ffaf00 대비 밝기 상승)
    203,  // 109  red bright      #ff5f5f  (기존 160 #d70000 대비 밝기 상승)
};

// Non-aligned 통일 회색, pair 110
inline constexpr int ALIGNED_NONALIGNED_DIM = 240;  // xterm-240 #585858
```

**대비 검증 (protein 0, olive 기준)**:
- Aligned:     xterm-148 `#afd700` = (175, 215, 0)
- Non-aligned: xterm-240 `#585858` = (88, 88, 88)
- **RGB 차이: (87, 127, 88) → 매우 뚜렷. 색조 + 밝기 모두 다름.**

#### 2. Screen.cpp — init_color_pairs() 확장

```cpp
// Aligned bright: pairs 101-109
for (int i = 0; i < 9; ++i)
    init_pair(i + 101, Palettes::PROTEIN_BRIGHT_COLORS[i], -1);
// Aligned non-aligned dim: pair 110
init_pair(110, Palettes::ALIGNED_NONALIGNED_DIM, -1);
```

#### 3. Screen.cpp — assign_colors_to_points() aligned 분기 수정

```cpp
} else if (screen_mode == "aligned") {
    int bright_id = (protein_idx % 9) + 101;  // pairs 101-109
    for (auto& pt : points) {
        pt.color_id = pt.is_aligned ? bright_id : 110;
    }
}
```

#### 4. Camera.cpp — 스크린샷 경로에 새 pair 범위 추가

`Camera::screenshot()` 내 `cid → xterm` 매핑에 101-109, 110 범위 추가.

---

## 작업 B: 점자(Braille) 모드 원근감 개선

### 문제 분석

현재 braille 모드의 depth cue 현황:

| 기법 | 상태 | 설명 |
|------|------|------|
| Perspective projection (x/z) | O | 원근 투영으로 먼 물체가 작게 보임 |
| Z-buffer occlusion | O | 가까운 점이 먼 점을 가림 |
| Depth character (`@#%*^-.`) | **X (braille에서 무시)** | `RenderPoint.pixel`에 저장되지만 braille 합성 시 사용 안 됨 |
| Depth-based color/fog | **X** | 색상이 모드별 고정, depth 무관 |
| Structure thickness 차등 | 부분적 | Coil=3px, Helix/Sheet=5px cross. 하지만 depth 무관 고정 |
| Alpha fade (스크린샷만) | O (터미널 X) | `Camera::get_alpha_from_depth()`는 PNG 전용 |

**핵심 문제**: 비-braille 모드에서는 `@#%*^-.` 문자의 밀도 차이가 원근감을 제공하지만, braille 모드에서는 모든 점이 동일한 "dot on" 상태라 이 depth cue가 완전히 사라짐. 유일한 depth cue는 perspective projection과 occlusion뿐.

### 해결: Depth-based Fog via Color Pair Bands

**핵심 아이디어**: 각 color mode의 색상에 대해 3단계 depth band(near/mid/far)를 두고, 같은 색조의 밝기 변형을 할당. 터미널 256색 범위 내에서 fog effect를 시뮬레이션.

#### 구현 전략

##### Phase 1: Depth Band 계산 (Screen.cpp::project)

기존 `depth_base_min_z`, `depth_base_max_z` 활용:

```cpp
// depth band 경계 (near=0, mid=1, far=2)
float depth_range = depth_base_max_z - depth_base_min_z;
float band_size = depth_range / 3.0f;

// RenderPoint에 depth_band 필드 추가
int depth_band;  // 0=near, 1=mid, 2=far
float t = (z - depth_base_min_z) / depth_range;
if      (t < 0.33f) depth_band = 0;
else if (t < 0.66f) depth_band = 1;
else                 depth_band = 2;
```

##### Phase 2: RenderPoint 확장

```cpp
struct RenderPoint {
    // ... 기존 필드 ...
    int depth_band = 0;  // 0=near(vivid), 1=mid(normal), 2=far(dim)
};
```

##### Phase 3: Fog Color Palette 설계

각 모드별 색상에 대해 3단계 밝기를 xterm-256 범위 내에서 선택:

**Protein 모드 예시 (olive, pair 1)**:
| Band | xterm | RGB | 시각적 |
|------|-------|-----|--------|
| near (0) | 142 | (175,175,0) | 밝은 올리브 |
| mid (1) | 100 | (135,135,0) | 기존 올리브 (변경 없음) |
| far (2) | 58 | (95,95,0) | 어두운 올리브 |

**이렇게 하면 기존 PROTEIN_COLORS가 mid band가 되고, near는 더 밝은 변형, far는 기존 DIM_COLORS를 재사용.**

##### Phase 4: Color Pair 할당 전략

ncurses는 최대 256 color pair를 지원. 현재 사용량:

| 범위 | 용도 | 개수 |
|------|------|------|
| 1-9 | Protein vivid | 9 |
| 11-19 | Protein dim | 9 |
| 21-35 | Chain | 15 |
| 41-42 | Structure (H/S) | 2 |
| 43-44 | Interface | 2 |
| 45-46 | Aligned | 2 |
| 51-70 | Rainbow | 20 |
| 71-74 | pLDDT | 4 |
| 75-84 | Conservation | 10 |
| 101-110 | Aligned bright + dim (작업 A) | 10 |
| **합계** | | **83** |

**Fog용 추가 pair 할당** (작업 A의 101-110 이후):

| 범위 | 용도 | 개수 |
|------|------|------|
| 120-128 | Protein near (bright) | 9 |
| 1-9 | Protein mid (기존 유지) | 0 (추가 없음) |
| 11-19 | Protein far (기존 dim 재사용) | 0 (추가 없음) |
| 130-144 | Chain near/far (15색 × pair 추가) | 15 |
| 150-169 | Rainbow near (20색) | 20 |
| 170-189 | Rainbow far (20색) | 20 |

**총 추가: ~64 pairs → 합계 ~147/256. 여유 있음.**

##### Phase 5: assign_colors_to_points에 depth_band 반영

```cpp
// Protein 모드 예시
if (screen_mode == "protein") {
    int idx = protein_idx % 9;
    for (auto& pt : points) {
        if (!use_braille || pt.depth_band == 1) {
            pt.color_id = idx + 1;       // mid: 기존 색상
        } else if (pt.depth_band == 0) {
            pt.color_id = idx + 120;     // near: 밝은 색상
        } else {
            pt.color_id = idx + 11;      // far: dim 색상
        }
    }
}
```

**비-braille 경로는 변경 없음** — depth character가 이미 원근감을 제공하므로 fog 불필요.

##### Phase 6: print_screen_braille 수정

braille 셀 내 color 선택 시, 기존 "best depth" 로직을 유지하되, 선택된 RenderPoint의 color_id에 이미 depth_band가 반영되어 있으므로 추가 변경 불필요.

#### 대안 검토 (채택하지 않은 방안들)

**A_BOLD / A_DIM 속성 사용**:
- 장점: Color pair 추가 불필요
- 단점: A_BOLD는 터미널마다 동작이 다름 (밝기 증가 vs 굵은체). A_DIM은 많은 터미널에서 미지원. 일관성 보장 불가.

**Dot 패턴 간소화 (far → 점 수 감소)**:
- 장점: 구현 간단
- 단점: 구조물의 형태 정보가 손실됨. Braille의 8-dot 해상도가 이미 낮은데 더 줄이면 시각 품질 저하.

**Structure thickness 차등 (depth 기반)**:
- 장점: 형태로 원근 표현
- 단점: Helix/Sheet의 구조적 의미와 depth 정보가 혼재. 사용자가 구조 두께를 depth로 오인할 수 있음.

---

## 구현 순서

### Step 1: RenderPoint 확장
- `RenderPoint.hpp`에 `int depth_band = 0;` 추가

### Step 2: Palette 확장
- `Palette.hpp`에 `PROTEIN_BRIGHT_COLORS[9]` (작업 A), `PROTEIN_NEAR_COLORS[9]` (작업 B), chain/rainbow near/far 배열 추가
- `ALIGNED_NONALIGNED_DIM` 상수 추가

### Step 3: init_color_pairs 확장
- 작업 A: pairs 101-110
- 작업 B: pairs 120-189 (fog bands)

### Step 4: project() — depth_band 계산
- braille 경로의 RenderPoint 생성 시 depth_band 설정
- `draw_line()` 내부에서도 보간된 z값으로 depth_band 계산

### Step 5: assign_colors_to_points — Aligned 모드 수정 (작업 A)
- aligned 분기를 bright(101-109) / dim(110)으로 변경

### Step 6: assign_colors_to_points — Fog 적용 (작업 B)
- braille 모드일 때 protein, chain, rainbow 모드에서 depth_band 기반 color pair 분기
- 비-braille 경로는 변경 없음

### Step 7: Camera.cpp 스크린샷 경로 동기화
- 새 color pair 범위에 대한 xterm 매핑 추가

### Step 8: 테스트
- `-m aligned` 모드에서 aligned/non-aligned 영역의 시각적 구분 확인
- braille 모드에서 회전 시 near/mid/far 색상 변화로 원근감 확인
- 비-braille 모드에서 기존 동작 변경 없음 확인
- 스크린샷 출력에서 새 color pair가 정확히 매핑되는지 확인

---

## 수정 대상 파일 요약

| 파일 | 변경 내용 |
|------|----------|
| `src/visualization/RenderPoint.hpp` | `depth_band` 필드 추가 |
| `src/visualization/Palette.hpp` | PROTEIN_BRIGHT_COLORS, PROTEIN_NEAR_COLORS, chain/rainbow fog 배열 추가 |
| `src/visualization/Screen.cpp` | init_color_pairs() 확장, project() depth_band 계산, assign_colors_to_points() 양쪽 작업 반영 |
| `src/visualization/Camera.cpp` | 새 pair 범위 xterm 매핑 |
