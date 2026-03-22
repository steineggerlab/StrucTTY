# 구현 계획 2: Far 색상 채도 복원 + Aligned 모드 -fs/-fm 버그 수정

---

## 문제 1: Far 색상이 완전한 회색이라 protein 간 구분 불가

### 원인 분석

plan_1에서 far 색상을 xterm grayscale(236-245)로 교체하면서, 모든 protein이 동일한 무채색 회색으로 렌더링됨.

현재 far 색상 (모두 순수 grayscale):
| Protein | Far xterm | RGB |
|---------|----------|-----|
| olive   | 242 | (108,108,108) |
| turquoise | 243 | (118,118,118) |
| navy    | 237 | (58,58,58) |
| purple  | 239 | (78,78,78) |
| pink    | 245 | (138,138,138) |
| coral   | 243 | (118,118,118) |
| brown   | 240 | (88,88,88) |
| orange  | 243 | (118,118,118) |
| red     | 240 | (88,88,88) |

**문제**: coral, orange, turquoise가 모두 243(#767676)으로 동일. brown, red도 240(#585858)으로 동일. 구분 불가.

### 해결 방안: 채도를 유지하면서 어두운 색상 사용

mid 색상(PROTEIN_COLORS)에서 2단계 아래이되, **색조(hue)를 유지**하는 xterm-256 색상 선택:

#### 1-A: PROTEIN_FAR_COLORS (pairs 200-208)

| Protein | mid xterm | mid RGB | **새 far xterm** | **새 far RGB** | 효과 |
|---------|----------|---------|-----------------|---------------|------|
| olive   | 100 (135,135,0) | | **58** | **(95,95,0)** | 어두운 올리브 (색조 유지) |
| turquoise | 73 (95,175,175) | | **30** | **(0,135,135)** | 어두운 청록 |
| navy    | 57 (95,0,255) | | **18** | **(0,0,135)** | 어두운 남색 |
| purple  | 91 (135,0,175) | | **54** | **(95,0,135)** | 어두운 보라 |
| pink    | 175 (215,135,215) | | **96** | **(135,95,135)** | 어두운 분홍 |
| coral   | 173 (215,135,95) | | **130** | **(175,95,0)** | 어두운 산호 |
| brown   | 136 (175,135,0) | | **94** | **(135,95,0)** | 어두운 갈색 |
| orange  | 179 (215,175,95) | | **136** | **(175,135,0)** | 어두운 주황 |
| red     | 160 (215,0,0) | | **88** | **(135,0,0)** | 어두운 빨강 |

**핵심**: 각 색상의 hue를 유지하되 명도/채도를 낮춤. 순수 회색 대신 "muted/dark" 색조.

#### 1-B: CHAIN_FAR_COLORS (pairs 145-159)

동일 원칙: 기존 chain mid 색상에서 2단계 아래, 색조 유지.

| Chain | **새 far xterm** | **새 far RGB** |
|-------|-----------------|---------------|
| olive | 58 | (95,95,0) |
| turquoise | 30 | (0,135,135) |
| navy | 18 | (0,0,135) |
| purple | 54 | (95,0,135) |
| pink | 96 | (135,95,135) |
| coral | 130 | (175,95,0) |
| brown | 94 | (135,95,0) |
| orange | 136 | (175,135,0) |
| red | 88 | (135,0,0) |
| teal | 23 | (0,95,95) |
| lime | 64 | (95,135,0) |
| magenta | 90 | (135,0,95) |
| gold | 136 | (175,135,0) |
| lavender | 60 | (95,95,135) |
| salmon | 131 | (175,95,95) |

#### 1-C: RAINBOW_FAR (pairs 180-199)

기존 grayscale 대신, rainbow gradient의 어두운 버전:
```
124, 130, 136, 172, 178, 142, 106, 70, 34, 28,
23, 24, 24, 18, 18, 18, 54, 54, 90, 127
```

#### 1-D: pLDDT/conservation/interface/aligned far 색상도 동일 원칙

- pLDDT far: 어두운 파랑(18), 어두운 청록(30), 어두운 노랑(58), 어두운 주황(94) → 색조 유지
- conservation far: gradient의 어두운 버전 (색조 유지)
- interface far: 어두운 마젠타(90) / 어두운 회색(236)
- aligned far: 변경 없음 (non-aligned dim은 의도적으로 회색)

### 변경 파일

| 파일 | 변경 |
|------|------|
| `Palette.hpp` | PROTEIN_FAR_COLORS, CHAIN_FAR_COLORS, RAINBOW_FAR, PLDDT_FAR, CONSERVATION_FAR, INTERFACE_FAR_COLOR 값 교체 |

---

## 문제 2: -fs/-fm에서 `-m aligned` 모드의 정렬 하이라이트가 안됨

### 원인 분석

#### -ut에서는 왜 작동하는가?

```
1. set_utmatrix() → 모든 protein의 init_atoms에 U/T 적용 (apply_ut_to_init_atoms)
2. normalize_proteins() → screen_atoms만 정규화 (do_shift/do_scale)
3. compute_aligned_all() → init_atoms 기반 거리 계산
```

**핵심**: step 1에서 **모든 protein의 init_atoms**가 동일한 좌표계(UT 변환 후 Å 공간)에 있으므로, step 3의 거리 계산이 올바르게 작동.

#### -fs에서 왜 실패하는가?

`load_next_hit()` 흐름 (Screen.cpp:1350-1608):

```
1. target protein 로드 → init_atoms = 원래 PDB Å 좌표
2. target 정규화 → screen_atoms만 변경 (do_shift/do_scale)
3. Kabsch/Foldseek transform → apply_foldseek_transform(1, U, T, T_ang)
   → target의 screen_atoms: do_naive_rotation(U) + do_shift(T)
   → target의 init_atoms: apply_ut_to_init_atoms(U, T_ang)
4. compute_aligned_from_aln(qaln, taln, 5.0f) 호출
   → query의 init_atoms (index 0): 원래 PDB Å 좌표 (변환 없음!)
   → target의 init_atoms (index 1): U/T_ang 변환 적용됨
```

**버그**: `compute_aligned_regions_from_aln()`은 `init_atoms`의 좌표를 비교하는데:
- **query init_atoms**: 원래 Å 좌표 (변환 없음)
- **target init_atoms**: Foldseek U/T가 적용된 Å 좌표 (다른 좌표계!)

좌표계가 불일치하여 `qa->x - ta->x` 등의 거리 계산이 무의미해짐.

**근본 원인**: `apply_ut_to_init_atoms`는 **target에만** 적용되고, **query에는 적용되지 않음**. Foldseek/Kabsch의 U/T는 "target을 query에 맞추는" 변환이므로, 변환 후 target이 query 좌표계로 이동해야 합니다. 하지만 query의 init_atoms는 원래 PDB 좌표 그대로이고, T_ang는 정규화 역변환을 포함하므로 좌표계가 일치하지 않음.

#### T_ang 계산의 문제

```cpp
// Screen.cpp:1578-1583
T_ang[r] = T[r] / norm_scale + q_cen[r];
for (int c = 0; c < 3; c++) T_ang[r] -= U[r*3+c] * t_cen2[c];
```

이 공식은 "정규화 공간의 T를 Å 공간으로 역변환"하려는 의도지만:
- `T[r] / norm_scale`: 정규화 스케일 역변환
- `+ q_cen[r]`: query centroid 오프셋
- `- U * t_cen2`: target centroid의 회전 보정

결과적으로 target의 init_atoms는 **query의 원래 centroid 주변**으로 이동하지만, query의 init_atoms 좌표와 정확히 일치하지 않을 수 있음.

#### -fm에서도 동일한 문제

`apply_foldmason_superposition()` (Screen.cpp:1860-1882):
- `apply_foldseek_transform(target_idx, U, T, T_ang)` → target init_atoms만 변환
- `compute_aligned_regions_from_aln()` → query와 target init_atoms 비교 → 좌표계 불일치

### 해결 방안

#### 2-A: alignment string(qaln/taln) 기반 — 좌표 비교 제거

**핵심 관찰**: `-fs`/`-fm`에서는 **alignment string이 이미 존재**. qaln/taln에서 gap이 아닌 위치는 이미 "정렬된 잔기"로 간주할 수 있음.

현재 `compute_aligned_regions_from_aln()`은 alignment string + 거리 threshold의 이중 조건을 사용:
```cpp
if (!q_gap && !t_gap) {
    if (dx*dx + dy*dy + dz*dz < thr2) {  // 이 조건이 좌표계 불일치로 실패
        qa->is_aligned = true;
```

**해결**: `-fs`/`-fm` 경로에서는 alignment string만으로 is_aligned를 설정. 거리 threshold 조건을 **선택적으로 적용**:
- alignment string이 있으면 → string 기반 (gap이 아닌 위치 = aligned)
- alignment string이 없으면 → nearest-neighbor 기반 (기존 방식)

구체적으로, `compute_aligned_regions_from_aln()`에 **좌표 비교 skip 옵션** 추가:

```cpp
void Protein::compute_aligned_regions_from_aln(Protein& other,
                                               const std::string& qaln,
                                               const std::string& taln,
                                               float threshold,
                                               bool skip_distance_check = false);
```

`skip_distance_check = true`이면:
```cpp
if (!q_gap && !t_gap) {
    if (q_idx < q_size && t_idx < t_size) {
        qa->is_aligned = true;
        ta->is_aligned = true;
        // 거리 비교 없이 바로 설정
    }
}
```

`load_next_hit()`과 `apply_foldmason_superposition()`에서 호출 시 `skip_distance_check = true` 전달.

**장점**:
- Foldseek/FoldMason의 alignment string은 이미 구조적 정렬 결과이므로, 추가 거리 검증이 불필요
- 좌표계 문제를 근본적으로 우회
- 단순하고 안전

**단점**:
- 실제로 구조적으로 먼 잔기도 aligned로 표시될 수 있으나, Foldseek/FoldMason의 정렬 품질을 신뢰하면 문제 없음

#### 2-B: (대안) init_atoms 좌표계 통일

target의 init_atoms에 적용하는 T_ang 계산을 수정하여, query의 init_atoms와 동일한 좌표계에 있도록 보장. 하지만 이는:
- 정규화(norm_scale, centroid shift) 역변환 로직이 복잡
- T_ang 공식의 정확성을 검증하기 어려움
- 다른 곳에서 init_atoms를 사용하는 코드에 부작용 가능

**방안 2-A가 더 안전하고 깔끔.**

### 변경 파일

| 파일 | 변경 |
|------|------|
| `Protein.hpp` | `compute_aligned_regions_from_aln()` 시그니처에 `skip_distance_check` 파라미터 추가 |
| `Protein.cpp` | `compute_aligned_regions_from_aln()`: `skip_distance_check` 분기 추가 |
| `Screen.cpp` | `load_next_hit()`, `apply_foldmason_superposition()`에서 `skip_distance_check = true` 전달 |

---

## 구현 순서

### Step 1: Far 색상 채도 복원 (문제 1) ✅ 완료
- `Palette.hpp`: PROTEIN_FAR_COLORS → 색조 유지 어두운 색 (58,30,18,54,96,130,94,136,88)
- CHAIN_FAR_COLORS → 동일 원칙 15색 교체
- RAINBOW_FAR → rainbow gradient의 어두운 버전으로 교체
- PLDDT_FAR → 어두운 파랑(18), 청록(30), 노랑(58), 주황(94)
- CONSERVATION_FAR → gradient 어두운 버전
- INTERFACE_FAR_COLOR → 어두운 마젠타(90)

### Step 2: Aligned 모드 -fs/-fm 버그 수정 (문제 2) ✅ 완료
- `Protein.hpp/cpp`: `compute_aligned_regions_from_aln()`에 `skip_distance_check` 파라미터 추가
- `Screen.cpp`: load_next_hit() -fs 경로 3곳에서 `skip_distance_check=true` 전달
- `Screen.cpp`: apply_foldmason_superposition() -fm 경로에서 `skip_distance_check=true` 전달

### Step 3: 빌드 및 검증 ✅ 완료
- 컴파일 에러 없이 빌드 성공

---

## 수정 대상 파일 요약

| 파일 | 변경 내용 |
|------|----------|
| `src/visualization/Palette.hpp` | far 색상을 색조 유지 어두운 색으로 교체 |
| `src/structure/Protein.hpp` | compute_aligned_regions_from_aln 시그니처 변경 |
| `src/structure/Protein.cpp` | skip_distance_check 분기 추가 |
| `src/visualization/Screen.cpp` | -fs/-fm에서 skip_distance_check=true 전달 |
