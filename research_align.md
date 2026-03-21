# `-m aligned` 작동 불량 원인 분석

## 증상

`-fs -m aligned` 또는 `--foldmason -m aligned` 실행 시 모든 잔기가 dim 색(비정렬)으로 표시되며, 정렬 구별이 전혀 이뤄지지 않음.

---

## 핵심 원인: `init_atoms`가 정규화되지 않는 상태에서 정규화 공간 T 벡터를 적용

### 좌표계 흐름 정리

**`normalize_proteins()` 이후 상태:**

| 변수 | 좌표계 |
|------|--------|
| `screen_atoms[0]` | 정규화됨: `(x_pdb - cx) * scale` |
| `init_atoms[0]` | **원래 PDB 좌표 (Å 단위, 미변환)** |

근거: `do_shift()` (Protein.cpp:424), `do_scale()` (Protein.cpp:435)는 모두 `screen_atoms`만 변환한다. `init_atoms`에는 아무것도 적용하지 않는다.

---

### `load_next_hit()` 이후 target protein 상태

`load_next_hit()` (Screen.cpp:1543-1551)의 target 정규화:

```cpp
target_protein->do_shift(t_shift);    // screen_atoms만 변환
target_protein->do_scale(norm_scale); // screen_atoms만 변환
```

그 다음 `apply_foldseek_transform(1, U_norm, T_norm)` (Screen.cpp:1626):

```cpp
void Screen::apply_foldseek_transform(int protein_idx, const float* U_flat, const float* T) {
    data[protein_idx]->do_naive_rotation(U_flat); // screen_atoms만
    data[protein_idx]->do_shift(T);               // screen_atoms만
    data[protein_idx]->apply_ut_to_init_atoms(U_flat, T); // ← 문제
    yesUT = true;
}
```

`apply_ut_to_init_atoms(U, T)` (Protein.cpp:547)는 다음을 실행:

```cpp
atom.x = U[0]*x + U[1]*y + U[2]*z + T[0];  // init_atoms에 적용
```

**이때 전달되는 U_norm, T_norm은 정규화 좌표계 기준으로 계산된 값이다.**

29컬럼 포맷의 경우 (Screen.cpp:1564-1576):
```cpp
T[0] = (Utc[0] + Tf[0] - norm_cx) * norm_scale;  // ← 정규화 공간 T
```

결과:

| 변수 | 실제 내용 |
|------|-----------|
| `screen_atoms[1]` | 정규화됨 + U_norm·x + T_norm → **정확** (표시용) |
| `init_atoms[1]` | U_norm · x_pdb_Å + T_norm → **좌표계 혼용 (오류)** |

- `init_atoms[1]`의 x, y, z는 원래 PDB Å 단위 (~100Å 규모)
- T_norm은 정규화 공간 단위 (~0.1–2.0 규모)
- 두 단위를 더하면 의미 없는 값이 됨

---

### `compute_aligned_regions_from_aln()` 에서의 거리 계산 실패

```cpp
void Protein::compute_aligned_regions_from_aln(Protein& other, ..., float threshold) {
    // init_atoms 사용 (Protein.cpp:616, 624)
    for (auto& [cid, chain] : init_atoms) ...           // query: 원래 PDB Å
    for (auto& [cid, chain] : other.init_atoms) ...     // target: 혼용 좌표

    float dx = qa->x - ta->x;
    ...
    if (dx*dx + dy*dy + dz*dz < thr2) { ... }  // thr2 = 5.0^2 = 25
}
```

- `init_atoms[0]` (query): 원래 PDB 좌표 (예: x ≈ 100Å)
- `init_atoms[1]` (target): U · x_Å + T_norm (예: x ≈ 100Å + 1.0 → 혼용)
- 두 단백질의 원점(centroid)이 다르므로 offset만 수백 Å
- 거리 ≫ threshold(5Å) → **모든 잔기 쌍이 non-aligned로 판정**
- `is_aligned = true`로 설정되는 잔기 없음 → 전체 dim 색 표시

---

## 동일 버그가 발생하는 케이스

### 1. Foldseek 29컬럼 포맷 (`has_transform == true`)

- `U_norm = hit.U`, `T_norm = (U·t_centroid + T_fs - q_centroid) * scale`
- `apply_ut_to_init_atoms(hit.U, T_norm)` → T가 정규화 공간값이므로 init_atoms에 오적용

### 2. Foldseek alis/kabsch 포맷 (`is_alis_format`)

```cpp
P_norm = (hit.alns[...] - norm_cx) * norm_scale;  // 정규화된 query coords
Q_norm = target_cas[...];                          // screen_atoms (정규화됨)
kabsch(P_norm, Q_norm, ..., U, T);                 // → U, T 모두 정규화 공간
apply_foldseek_transform(1, U, T);                 // init_atoms에 정규화 T 오적용
```

### 3. FoldMason `apply_foldmason_superposition()`

```cpp
// Screen.cpp:1668-1676
for (auto& [cid, chain] : data[query_protein_idx]->get_atoms()) { // screen_atoms (정규화됨)
    query_cas.push_back({a.x, a.y, a.z});
}
...
kabsch(P_flat, Q_flat, N, U, T);          // U, T는 정규화 공간
apply_foldseek_transform(target_idx, U, T); // init_atoms에 정규화 T 오적용
```

이후 `compute_aligned_regions_from_aln(query, target, fm_aa, fm_aa, 5.0f)` → 동일하게 거리 계산 실패.

---

## 부차적 문제

### query `init_atoms[0]`도 정규화되지 않음

`normalize_proteins()` 후 `init_atoms[0]`은 원래 PDB Å 좌표다. 만약 init_atoms에 정규화가 올바르게 적용되어 `init_atoms[1]`이 normalized query frame에 놓이더라도, `init_atoms[0]`도 normalized 공간에 있어야 비교가 가능하다. 현재는 둘 다 원래 PDB 공간이지만 서로 다른 단백질의 원점 기준이라 그것만으로도 거리가 수백 Å.

실질적 의미: 두 단백질의 Foldseek transform이 완벽히 같은 Å 좌표계(query 원점 기준)로 target을 이동시키는 올바른 U/T를 init_atoms에 적용한다면 작동하지만, 현재 코드는 정규화된 T를 사용하므로 실패한다.

---

## 올바른 수정 방향 (우선순위 순)

### Option A — 최소 수정: `apply_ut_to_init_atoms`에 올바른 (Å 공간) U/T 전달

**Foldseek 29컬럼 (`has_transform`):**
```
현재: apply_ut_to_init_atoms(hit.U, T_norm)   ← T_norm = 정규화 공간
수정: apply_ut_to_init_atoms(hit.U, hit.T)    ← hit.T = 원래 Å 공간
```
init_atoms[0]이 query 원래 PDB 좌표, init_atoms[1]이 `U·x_pdb + hit.T` (query PDB 공간으로 이동한 target) → 거리 계산 가능.

**Foldseek kabsch / FoldMason:**
kabsch 결과 U, T를 Å 공간으로 역변환:
```
T_Å = T_norm / norm_scale + q_centroid - U · t_centroid
```
이 T_Å를 `apply_ut_to_init_atoms(U, T_Å)`에 사용.

### Option B — 근본 수정: `init_atoms`에도 정규화 적용

`Protein`에 메서드 추가:
```cpp
void apply_normalize_to_init_atoms(float cx, float cy, float cz, float scale);
```

`normalize_proteins()` 및 `load_next_hit()` 에서 screen_atoms와 동일하게 init_atoms를 centroid shift + scale 적용.
그러면 `apply_ut_to_init_atoms(U_norm, T_norm)` 이 올바르게 동작하며 threshold는 정규화 단위로 변경 필요 (~`5.0f * norm_scale`).

### Option C — 간단하나 구조 변경 필요: `compute_aligned_regions_from_aln/nn`에서 screen_atoms 사용

screen_atoms는 모든 변환(normalize + U/T rotation)이 올바르게 적용된 상태.
`compute_aligned_regions_from_aln/nn`이 `init_atoms` 대신 `screen_atoms`로 거리 계산하게 변경.
threshold = 정규화 단위로 변환 필요.

단점: screen_atoms에는 helix/sheet geometry atom이 섞여있어 1:1 잔기 매핑이 어려움.

---

## 결론

모든 경우의 핵심 버그는 **`apply_foldseek_transform()` 이 `apply_ut_to_init_atoms(U_norm, T_norm)` 을 호출할 때 정규화 좌표계의 T 값을 Å 단위의 init_atoms에 그대로 적용하는 것**이다.

- query `init_atoms[0]`: 원래 PDB Å 좌표
- target `init_atoms[1]` (변환 후): U · x_Å + T_norm → 수백 Å 오차
- `compute_aligned_regions_from_aln(threshold=5.0Å)`: 항상 실패 → 전체 non-aligned

**권장 수정**: Option A (최소 침습). Foldseek 29컬럼은 hit.T를 직접 사용, kabsch 경로는 T를 Å 공간으로 역변환하여 `apply_ut_to_init_atoms` 에 전달.
