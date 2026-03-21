# StrucTTY 기능 수정 및 추가 구현 계획

## 개요

4가지 기능 수정/추가 작업에 대한 상세 구현 계획이다.
각 항목별로 현재 상태 분석, 문제점, 구현 방안, 수정 파일 목록을 기술한다.

---

## 작업 1: -fs / -fm 모드에서 query-target 구조 정렬(superposition) 보장

### 현재 상태

`Screen::load_next_hit()` (Screen.cpp:1484)에서 Foldseek hit 로드 시:

- **29컬럼 포맷** (`has_transform == true`): Foldseek 제공 U(3×3)/T(3) 를 직접 사용하여 target을 query 좌표계로 변환. **정상 동작.**
- **21컬럼 alis 포맷** (`is_alis_format && !alns.empty()`): `alns` CA 좌표(query frame)와 target CA 좌표 간 Kabsch SVD로 U/T 역산. **정상 동작.**
- **17컬럼 포맷** (`has_aln == true, has_transform == false, is_alis_format == false`): qaln/taln 정보는 있지만 U/T도 alns도 없음. → `computed_transform = false` → **superposition 없이 target이 원래 좌표에 그냥 표시됨.**
- **12컬럼 포맷** (`has_aln == false, has_transform == false`): 아무 정렬 정보 없음. → **superposition 없이 표시.**

`Screen::apply_foldmason_superposition()` (Screen.cpp:1676)에서 FoldMason hit 로드 시:

- FoldMason JSON은 `ca` 좌표가 있어 Kabsch SVD 가능. **정상 동작.**
- FoldMason FASTA는 `ca` 좌표 없음 → `build_aligned_pairs()`가 MSA 컬럼 기준 잔기 쌍을 반환하고, **screen_atoms의 정규화된 좌표**로 Kabsch 수행. **정상 동작 (정규화 이후 좌표 사용).**

### 문제점

1. **17컬럼**: qaln/taln이 있으므로 aligned 잔기 쌍을 알 수 있으나, 좌표 변환을 계산하지 않고 있음.
2. **12컬럼**: alignment 정보 자체가 없으므로, 최소한 전체 CA 좌표 매칭으로라도 superposition을 시도해야 함.

### 구현 방안

#### 1-A. 17컬럼 포맷: qaln/taln 기반 Kabsch

`Screen::load_next_hit()` 내부, `hit.has_transform` 블록 이후에 **새 분기 추가**:

```
else if (hit.has_aln && !hit.is_alis_format) {
    // 17컬럼 포맷: qaln/taln로 aligned CA 쌍 추출 → Kabsch SVD
}
```

**로직:**
1. query protein (data[0]) 의 init_atoms에서 CA flat 리스트 추출 (chain 순서대로)
2. target protein (data[1]) 의 init_atoms에서 CA flat 리스트 추출
3. qaln/taln 순회: 양쪽 비-갭 위치에서 `(qstart-1 + q_idx)`, `(tstart-1 + t_idx)` 인덱스로 CA 좌표 쌍 수집
4. 수집된 쌍을 정규화 공간으로 변환 후 `kabsch()` 호출
5. `T_ang` 계산 (기존 alis 포맷과 동일한 공식)

**수정 파일:** `src/visualization/Screen.cpp` — `load_next_hit()` 함수

#### 1-B. 12컬럼 포맷: 전체 CA 기반 Kabsch (fallback)

17컬럼 분기 이후에 **최종 fallback 분기 추가**:

```
else if (!computed_transform) {
    // 12컬럼 또는 기타: 전체 CA 순서 매칭으로 Kabsch (min(qlen, tlen) 쌍 사용)
}
```

**로직:**
1. query/target init_atoms에서 CA flat 리스트 추출
2. `min(qlen, tlen)` 개수만큼 순서대로 1:1 매칭
3. 정규화 공간 변환 → `kabsch()` → `T_ang` 계산
4. 이 fallback은 서열 순서 매칭이므로 정확도는 낮지만 "아무렇게나 띄우는" 것보다 낫다

**수정 파일:** `src/visualization/Screen.cpp` — `load_next_hit()` 함수

#### 1-C. -fm FASTA에서도 superposition 보장 확인

현재 `apply_foldmason_superposition()`은 screen_atoms 좌표로 Kabsch를 수행하므로 JSON/FASTA 모두 동작한다. 단, FASTA에서 `build_aligned_pairs()`가 빈 결과를 반환하면 superposition이 생략됨. 이 경우도 fallback (전체 CA 순서 매칭) 추가.

**수정 파일:** `src/visualization/Screen.cpp` — `apply_foldmason_superposition()` 함수

### 수정 파일 요약

| 파일 | 수정 내용 |
|------|-----------|
| `src/visualization/Screen.cpp` | `load_next_hit()` 에 17컬럼/12컬럼 Kabsch 분기 추가, `apply_foldmason_superposition()` fallback 추가 |

---

## 작업 2: -m aligned 모드에서 단백질별 고유 색상 유지

### 현재 상태

`Screen::assign_colors_to_points()` (Screen.cpp:449):

```cpp
} else if (screen_mode == "aligned") {
    for (auto& pt : points) {
        pt.color_id = pt.is_aligned ? 45 : 46;  // 정렬됨(초록) or dim(olive)
    }
}
```

모든 단백질이 동일한 초록/olive dim 색으로 표시되어 어떤 구조가 query이고 어떤 구조가 target인지 구별이 불가능하다.

### 문제점

- `-m default` (protein 모드)는 단백질마다 다른 색(PROTEIN_COLORS[0..8])을 사용.
- `-m aligned`는 단백질 구분 없이 aligned=초록, non-aligned=olive dim으로 통일.
- 사용자 요구: **단백질별 고유 색을 유지하면서** aligned 부분만 시각적으로 구분.

### 구현 방안

**색상 전략:** 각 단백질의 protein color를 기본으로 하되, aligned 부분은 **밝은(vivid) 버전**, non-aligned 부분은 **어두운(dim) 버전**으로 표시한다. 이미 `Palette.hpp`에 `PROTEIN_COLORS[9]` (vivid)와 `PROTEIN_DIM_COLORS[9]` (dim)이 정의되어 있으므로 이를 활용한다.

**수정 내용 — `assign_colors_to_points()`:**

기존:
```cpp
} else if (screen_mode == "aligned") {
    for (auto& pt : points) {
        pt.color_id = pt.is_aligned ? 45 : 46;
    }
}
```

변경:
```cpp
} else if (screen_mode == "aligned") {
    int vivid_id = (protein_idx % 9) + 1;   // pairs 1-9 (PROTEIN_COLORS)
    int dim_id   = (protein_idx % 9) + 11;  // pairs 11-19 (PROTEIN_DIM_COLORS)
    for (auto& pt : points) {
        pt.color_id = pt.is_aligned ? vivid_id : dim_id;
    }
}
```

이렇게 하면:
- protein 0 (olive): aligned = 밝은 olive(100), non-aligned = 어두운 olive(58)
- protein 1 (turquoise): aligned = 밝은 turquoise(80), non-aligned = 어두운 turquoise(30)
- 최대 9종류 구분 가능 (기존 PROTEIN_COLORS 팔레트 재활용)

**`-s` (show_structure) 플래그 연동:**

`-m aligned -s` 조합 시에도 이차 구조(H/S) 색상을 적용할지 결정해야 한다. 현재 `protein+-s` 모드에서는 H=yellow(41), S=cyan(42), coil=dim으로 오버라이드한다. aligned 모드에서는 **이차 구조 오버라이드를 적용하지 않고** 단백질별 vivid/dim만 사용하는 것이 명확하다. (이미 현재 코드에서 `screen_mode == "protein"` 일 때만 구조 오버라이드를 적용하므로 추가 수정 불필요.)

### 수정 파일 요약

| 파일 | 수정 내용 |
|------|-----------|
| `src/visualization/Screen.cpp` | `assign_colors_to_points()` 의 aligned 분기를 단백질별 vivid/dim 색상으로 변경 |

---

## 작업 3: -fs 다중 타겟 파일 시 m8 필터링

### 현재 상태

`structty.cpp` (line 71-81):

```cpp
if (!params.get_foldseek_file().empty()) {
    FoldseekParser fs_nav_parser;
    if (fs_nav_parser.load(params.get_foldseek_file()) && fs_nav_parser.hit_count() > 0) {
        screen.set_foldseek_hits(fs_nav_parser.get_hits());
        screen.set_fs_db_path(params.get_db_path());
        if ((int)params.get_in_file().size() <= 1) {
            screen.load_next_hit(+1);
        }
    }
}
```

- **단일 입력 (query.pdb만)**: m8의 모든 hit을 로드하고 N/P 키로 순회.
- **다중 입력 (query.pdb target1.pdb target2.pdb ...)**: hit을 설정하지만 `load_next_hit`을 호출하지 않음. m8의 모든 hit이 그대로 저장됨.

### 문제점

사용자가 `structty query.pdb target1.pdb target2.pdb -fs result.m8` 형식으로 실행할 때:
- 이미 CLI로 target 파일을 지정했으므로, m8에서 **해당 target들에 대한 정보만** 가져와야 함.
- 현재는 m8의 모든 hit이 설정되어 의미 없는 hit도 포함됨.
- 또한, 이미 로드된 target protein들에 대해 m8의 U/T를 적용하여 superposition을 수행해야 함.

### 구현 방안

#### 3-A. m8 hit 필터링

`structty.cpp`에서 m8 로드 후, CLI로 지정된 target 파일명과 m8의 target 컬럼을 매칭하여 필터링.

**매칭 로직:**
1. CLI target 파일명에서 stem(확장자 제외 파일명) 추출: `target1.pdb` → `target1`
2. m8 hit의 `target` 필드와 비교:
   - 완전 일치 (`hit.target == stem`)
   - 부분 포함 (`hit.target.find(stem) != npos` 또는 `stem.find(hit.target) != npos`)
   - PDB ID 매칭: stem에서 PDB 4자리 추출 → hit.target의 PDB ID와 비교
3. 매칭되는 hit만 필터링하여 `screen.set_foldseek_hits()`에 전달

**코드 위치:** `structty.cpp` — m8 로드 블록

```cpp
if ((int)params.get_in_file().size() > 1) {
    // CLI target 파일들의 stem 리스트 구성
    std::vector<std::string> target_stems;
    for (int i = 1; i < (int)params.get_in_file().size(); i++) {
        fs::path p(params.get_in_file(i));
        target_stems.push_back(p.stem().string());
    }

    // m8 hit 중 target_stems에 매칭되는 것만 필터링
    const auto& all_hits = fs_nav_parser.get_hits();
    std::vector<FoldseekHit> filtered_hits;
    for (const auto& hit : all_hits) {
        for (const auto& stem : target_stems) {
            if (hit.target.find(stem) != std::string::npos ||
                stem.find(hit.target) != std::string::npos) {
                filtered_hits.push_back(hit);
                break;
            }
        }
    }
    screen.set_foldseek_hits(filtered_hits);
}
```

#### 3-B. 이미 로드된 target들에 m8의 U/T 적용

다중 입력 모드에서는 target protein이 이미 `data[1], data[2], ...`에 로드되어 있다. m8에서 필터링된 hit의 U/T 정보를 각 target에 적용하여 query와 정렬한다.

**로직:**
1. 필터링된 hit 리스트를 순회
2. 각 hit에 대해, 매칭되는 target protein index를 찾음
3. hit의 U/T(또는 Kabsch)를 해당 protein에 `apply_foldseek_transform()` 으로 적용
4. aligned 모드일 경우 `compute_aligned_regions_from_aln()` 호출

**코드 위치:** `structty.cpp` — 새로운 함수 또는 인라인 블록

```cpp
// 각 target에 대해 매칭되는 hit의 transform 적용
for (int ti = 1; ti < (int)params.get_in_file().size(); ti++) {
    fs::path tp(params.get_in_file(ti));
    std::string tstem = tp.stem().string();

    for (const auto& hit : filtered_hits) {
        if (hit.target.find(tstem) != std::string::npos ||
            tstem.find(hit.target) != std::string::npos) {
            // hit의 U/T를 data[ti]에 적용 (load_next_hit의 transform 로직 재사용)
            // → 별도 함수로 추출 필요
            break;
        }
    }
}
```

#### 3-C. transform 적용 로직 추출

현재 `load_next_hit()`에 있는 U/T 계산 및 적용 로직(29컬럼/21컬럼/17컬럼/fallback)을 **별도 함수**로 추출하여 재사용.

```cpp
// Screen.hpp에 선언
void apply_hit_transform(int target_protein_idx, const FoldseekHit& hit);
```

이 함수는 `load_next_hit()`과 다중 타겟 모드에서 공통으로 사용.

### 수정 파일 요약

| 파일 | 수정 내용 |
|------|-----------|
| `src/structty.cpp` | m8 hit 필터링 로직 추가, 다중 target에 transform 적용 |
| `src/visualization/Screen.hpp` | `apply_hit_transform()` 선언 추가 |
| `src/visualization/Screen.cpp` | `load_next_hit()`에서 transform 로직을 `apply_hit_transform()`으로 추출 |

---

## 작업 4: PDB100 데이터베이스 target ID 파싱 수정

### 현재 상태

`PDBDownloader::detect_db_type()` (PDBDownloader.cpp:90) 의 PDB 패턴:

```cpp
// 패턴 1: PDB — [0-9][a-z0-9]{3}_[A-Za-z0-9]+
if (target_id.size() >= 6 &&
    target_id[0] >= '0' && target_id[0] <= '9' &&
    is_alnum_lower_or_digit(target_id[1]) &&
    is_alnum_lower_or_digit(target_id[2]) &&
    is_alnum_lower_or_digit(target_id[3]) &&
    target_id[4] == '_' &&
    target_id.size() >= 6) {
    return DBType::PDB;
}
```

`FoldseekParser::extract_target_id()` (FoldseekParser.cpp:39):

```cpp
static std::string extract_target_id(const std::string& raw) {
    size_t pos = raw.find(' ');
    if (pos == std::string::npos) return raw;
    return raw.substr(0, pos);
}
```

### 문제점

사용자가 제공한 PDB100 검색 결과 예시:

```
job_A	3a0d-assembly1.cif.gz_A	Crystal Structure of Polygonatum ...
```

tab으로 분리하면 col[1] = `"3a0d-assembly1.cif.gz_A Crystal Structure of Polygonatum ..."`.

1. **`extract_target_id()`에서 첫 공백 분리**: `"3a0d-assembly1.cif.gz_A"` 추출. 여기까지는 정상.
2. **`detect_db_type("3a0d-assembly1.cif.gz_A")`**: PDB 패턴 `[0-9][a-z0-9]{3}_[A-Za-z0-9]+`와 **불일치**.
   - `target_id[4]` = `'-'` (하이픈) ≠ `'_'` (언더스코어)
   - → `DBType::Unknown` 반환 → 다운로드 실패

PDB100 포맷: `{pdbid}-assembly{N}.cif.gz_{chain}` (하이픈 구분, `.cif.gz` 확장자 포함)

3. **`extract_pdb_id("3a0d-assembly1.cif.gz_A")`**: 현재 `find('_')` → 위치 24 → `substr(0, min(24, 4))` = `"3a0d"`. **PDB ID 추출은 우연히 정상.** (첫 4글자가 PDB ID이므로)
4. **`extract_chain("3a0d-assembly1.cif.gz_A")`**: `find('_')` → `"assembly1.cif.gz_A"` → `rfind('_')` → `"A"`. **체인 추출은 우연히 정상.**

결론: 핵심 문제는 `detect_db_type()`이 PDB100 포맷을 인식하지 못하는 것이다.

### 구현 방안

#### 4-A. PDB100 패턴을 detect_db_type()에 추가

기존 PDB 패턴 검사를 확장하여 PDB100 포맷도 인식:

```cpp
// 패턴 1a: PDB100 — [0-9][a-z0-9]{3}-assembly.*\.cif\.gz_[A-Za-z0-9]+
// 예: 3a0d-assembly1.cif.gz_A
if (target_id.size() >= 6 &&
    target_id[0] >= '0' && target_id[0] <= '9' &&
    is_alnum_lower_or_digit(target_id[1]) &&
    is_alnum_lower_or_digit(target_id[2]) &&
    is_alnum_lower_or_digit(target_id[3]) &&
    target_id[4] == '-' &&
    contains(target_id, ".cif.gz_")) {
    return DBType::PDB;
}
```

이 패턴을 기존 PDB 패턴(`target_id[4] == '_'`) 바로 앞에 추가한다. `extract_pdb_id()`와 `extract_chain()`은 PDB100에서도 이미 정상 동작하므로 수정 불필요.

#### 4-B. extract_chain() 보강 (안전장치)

현재 `extract_chain()`이 PDB100에서 우연히 동작하지만, `rfind('_')`에 의존하는 것이 fragile하다. 명시적으로 `.cif.gz_` 패턴을 처리하도록 보강:

```cpp
std::string PDBDownloader::extract_chain(const std::string& target_id, DBType db_type) {
    if (db_type != DBType::PDB) return "-";

    // PDB100 형식: 3a0d-assembly1.cif.gz_A
    size_t cifgz_pos = target_id.find(".cif.gz_");
    if (cifgz_pos != std::string::npos) {
        return target_id.substr(cifgz_pos + 8);  // ".cif.gz_" = 8글자
    }

    // 기존 로직: 2xyz_B, 1a0n_MODEL_1_B
    size_t pos = target_id.find('_');
    if (pos == std::string::npos) return "-";
    std::string after = target_id.substr(pos + 1);
    size_t last_us = after.rfind('_');
    if (last_us != std::string::npos) after = after.substr(last_us + 1);
    return after.empty() ? "-" : after;
}
```

#### 4-C. extract_pdb_id() 보강 (안전장치)

PDB100 형식에서 하이픈 앞 4글자가 PDB ID임을 명시적으로 처리:

```cpp
std::string PDBDownloader::extract_pdb_id(const std::string& target_id) {
    // PDB100 형식: 3a0d-assembly1.cif.gz_A
    size_t dash_pos = target_id.find('-');
    if (dash_pos == 4) return target_id.substr(0, 4);

    // 기존: 2xyz_B → 2xyz
    size_t pos = target_id.find('_');
    if (pos == std::string::npos) return target_id.substr(0, 4);
    return target_id.substr(0, std::min(pos, (size_t)4));
}
```

### 수정 파일 요약

| 파일 | 수정 내용 |
|------|-----------|
| `src/structure/PDBDownloader.cpp` | `detect_db_type()`에 PDB100 패턴 추가, `extract_chain()`/`extract_pdb_id()` 보강 |

---

## 전체 수정 파일 요약

| # | 작업 | 수정 파일 |
|---|------|-----------|
| 1 | -fs/-fm superposition 보장 | `Screen.cpp` |
| 2 | -m aligned 단백질별 색상 | `Screen.cpp` |
| 3 | -fs 다중 타겟 필터링 | `structty.cpp`, `Screen.hpp`, `Screen.cpp` |
| 4 | PDB100 파싱 수정 | `PDBDownloader.cpp` |

## 구현 순서 (의존성 고려)

1. **작업 4** (PDB100 파싱) — 독립적, 가장 단순
2. **작업 2** (aligned 색상) — 독립적, 단순 색상 로직 변경
3. **작업 1** (superposition 보장) — `load_next_hit()` transform 로직 보강
4. **작업 3** (다중 타겟 필터링) — 작업 1의 transform 추출 함수에 의존
