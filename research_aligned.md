# -m aligned 모드 작동 불능 원인 분석

## 증상

작업 2에서 `assign_colors_to_points()`의 aligned 분기를 변경한 후, **모든 시나리오**에서 `-m aligned`가 시각적으로 구분 불가:

1. `-ut` 옵션 (Foldseek 미사용)
2. `-fs alis_BFVD.m8` (BFVD 결과 파일)
3. 기타 모든 `-m aligned` 사용 케이스

## 근본 원인: 색상 대비 부족

### 변경 전 (정상 동작)

```cpp
pt.color_id = pt.is_aligned ? 45 : 46;
```

| 상태 | color pair | xterm color | RGB | 시각적 |
|------|-----------|-------------|-----|--------|
| aligned | 45 | 46 (bright green) | (0, 255, 0) | **매우 밝은 초록** |
| non-aligned | 46 | 58 (olive dim) | (95, 95, 0) | 어두운 올리브 |

**대비: RGB 차이 (95, 160, 0) → 매우 뚜렷한 차이. 한눈에 구분 가능.**

### 변경 후 (작동 불능)

```cpp
int vivid_id = (protein_idx % 9) + 1;   // pairs 1-9 (PROTEIN_COLORS)
int dim_id   = (protein_idx % 9) + 11;  // pairs 11-19 (PROTEIN_DIM_COLORS)
pt.color_id = pt.is_aligned ? vivid_id : dim_id;
```

**protein_idx = 0 (query, 올리브):**

| 상태 | color pair | xterm color | RGB | 시각적 |
|------|-----------|-------------|-----|--------|
| aligned | 1 | 100 (olive) | (135, 135, 0) | 올리브 |
| non-aligned | 11 | 58 (olive dim) | (95, 95, 0) | 어두운 올리브 |

**대비: RGB 차이 (40, 40, 0) → 거의 구분 불가능!**

**protein_idx = 1 (target, 터콰이즈):**

| 상태 | color pair | xterm color | RGB | 시각적 |
|------|-----------|-------------|-----|--------|
| aligned | 2 | 80 (turquoise) | (95, 215, 215) | 터콰이즈 |
| non-aligned | 12 | 30 (turquoise dim) | (0, 135, 135) | 어두운 터콰이즈 |

**대비: RGB 차이 (95, 80, 80) → protein 0보다 나으나 여전히 subtle.**

## 왜 "작동 안 함"으로 보이는가

1. **PROTEIN_COLORS vs PROTEIN_DIM_COLORS는 "같은 색의 밝기 변형"** — 이들은 원래 `-m protein -s` 모드에서 coil(dim) 배경 위에 helix/sheet(yellow/cyan)을 강조하기 위한 보조 색상이다. vivid↔dim 간 대비가 주된 시각 정보가 아님.

2. **aligned 모드에서는 vivid↔dim이 유일한 시각 정보** — aligned/non-aligned 구분이 모드의 핵심인데, 그 구분을 위한 색상 대비가 너무 작다.

3. **특히 protein 0 (올리브)**: (135,135,0) vs (95,95,0)는 대부분의 터미널에서 동일하게 보임. 사실상 aligned 모드가 "없는 것"과 동일한 시각적 결과.

## 기능적(로직) 이상 여부

`is_aligned` 계산 로직 자체는 **변경 없음 — 정상 동작**:

- `compute_aligned_all()` → `compute_aligned_regions_nn()` → `sync_aligned_to_screen()` (경로 미변경)
- `compute_aligned_from_aln()` → `compute_aligned_regions_from_aln()` → `sync_aligned_to_screen()` (경로 미변경)
- `RenderPoint.is_aligned`는 `Atom.is_aligned`에서 정상 전파 (Screen.cpp:564, 578, 596, 703, 716, 830, 843)

**결론: is_aligned 값은 올바르게 계산됨. 문제는 100% 시각적 대비 문제.**

## 해결 방안

PROTEIN_COLORS vivid/dim 쌍은 aligned 모드의 주요 시각 지표로 사용하기에 대비가 부족하다. 해결 방안:

### 방안 A: 밝기(brightness) 기반 강조 추가

각 protein color의 vivid 버전 대신, 해당 색상에 **흰색 하이라이트**를 합성한 별도 밝은 색을 사용. 예: olive vivid를 xterm 100 대신 xterm 148 (175,215,0) 등 더 밝은 변형으로.

### 방안 B: aligned 부분에 bold/underline 등 속성 추가

ncurses `A_BOLD`를 aligned 부분에 추가하여 같은 색이어도 밝기 차이 생성.

### 방안 C: 단백질 색상 유지 + 밝기 단계 분리

새로운 palette 배열 `PROTEIN_BRIGHT_COLORS[9]`를 정의하여, 각 protein color의 **눈에 띄게 밝은 버전**을 할당. 기존 PROTEIN_DIM_COLORS와 충분한 대비 보장.

### 방안 D: 하이브리드 — protein 색조 + aligned 마커

aligned 부분은 protein 고유 색상의 밝은 버전, non-aligned 부분은 **통일된 회색조 dim**(예: xterm 240, 어두운 회색)으로 처리. 이렇게 하면 protein 색은 aligned 부분에서만 보이므로 구분이 극명.
