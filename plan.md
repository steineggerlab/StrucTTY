# StrucTTY 기능 추가 구현 계획

## 개요

7가지 신규 기능을 순서대로 구현한다. 각 기능은 독립적으로 동작 가능하도록 설계하되, 공통 인프라(Atom 필드 확장, Palette 확장)는 먼저 처리한다.
이때, 기본적으로 제공되지 않는 다른 패키지나 프로젝트를 활용하지 않도록 한다.
또한 코드의 병목 및 비효율적 병목 구조를 피하여 빠른 실행시간을 보장하도록 한다.

---

## 구현 순서

1. **공통 인프라** (Atom.hpp, Palette.hpp 확장) - 모든 기능의 전제 조건
2. **기능 7** - 코일 렌더링 버그 수정 (가장 빠르게 눈에 띄는 개선)
3. **기능 2** - pLDDT 색상 표시
4. **기능 1** - Interface region 색상 표시
5. **기능 4** - UTMatrix 정렬 구조 색상 표시
6. **기능 5** - MSA Conservation Score 색상 표시
7. **기능 6** - 커서(마우스) 기반 잔기 정보 패널 표시
8. **기능 3** - Foldseek 결과 파일 실시간 Hit 탐색

---

## 단계 0: 공통 인프라

### 0-1. `src/structure/Atom.hpp` 필드 확장

현재:
```cpp
struct Atom {
    float x, y, z;
    char structure = 'x';
};
```

추가할 필드:
```cpp
struct Atom {
    float x, y, z;
    char structure = 'x';

    // 기능 2: pLDDT
    float bfactor = 0.0f;

    // 기능 1: interface region
    bool is_interface = false;

    // 기능 4: UTMatrix 정렬 구조
    bool is_aligned = false;

    // 기능 5: MSA conservation
    float conservation_score = -1.0f;  // -1 = 미설정

    // 기능 6: 잔기 정보 hover
    int residue_number = -1;
    std::string residue_name = "";
};
```

### 0-2. `src/visualization/RenderPoint.hpp` 필드 확장

커서 hover 기능을 위해 잔기 번호와 이름도 저장:
```cpp
struct RenderPoint {
    // 기존 필드 유지
    int x, y;
    float depth;
    char pixel;
    int color_id;
    std::string chainID;
    char structure;

    // 기능 6: 잔기 정보
    int residue_number = -1;
    std::string residue_name = "";

    // 기능 2: pLDDT
    float bfactor = 0.0f;

    // 기능 5: conservation
    float conservation_score = -1.0f;
};
```

### 0-3. `src/visualization/Palette.hpp` 색상 확장

신규 color pair 할당:
```
43: interface residue 강조색 (밝은 마젠타, xterm-201)
44: non-interface dim (기존 dim 색상 재활용 가능)
45: aligned region 강조색 (밝은 초록, xterm-46)
46: non-aligned dim
71: pLDDT Very High (>90) - 파란색 (xterm-21)
72: pLDDT Confident (70-90) - 청록색 (xterm-51)
73: pLDDT Low (50-70) - 노란색 (xterm-226)
74: pLDDT Very Low (<50) - 주황색 (xterm-208)
75-84: conservation 10단계 그래디언트 (파랑→빨강)
```

`init_colors()` 함수에 위 pair들을 추가 초기화한다.

### 0-4. `src/structure/Parameters.hpp` CLI 인수 확장

신규 인수:
```
-m plddt          기능 2: pLDDT 색상 모드
-m interface      기능 1: interface region 색상 모드
-m conservation   기능 5: MSA conservation 색상 모드
-m aligned        기능 4: 정렬 구조 색상 모드
--msa <file>      기능 5: MSA 파일 경로
-fs <file>        기능 3: Foldseek 결과 파일 경로
--db-path <dir>   기능 3: target PDB 파일 디렉토리
```

기존 `-m` 옵션의 유효값 목록에 위 4가지 추가.

---

## 기능 7: 헬릭스/시트 영역의 Coil 제거 + Coil 굵기 감소

### 문제 원인
`StructureMaker::calculate_ss_points()`에서 helix/sheet 지오메트리를 생성하지만, 해당 잔기들의 원본 CA 원자들도 coil 리스트에 그대로 남아 있어 이중으로 그려진다.

### 수정 대상: `src/structure/StructureMaker.cpp`

`calculate_ss_points()` 함수 내에서 segment 분류 로직 수정:

```
현재 흐름:
  전체 chain atoms 순회
    → Helix segment → 나선 지오메트리 생성 (ss_atoms에 추가)
    → Sheet segment → 리본 지오메트리 생성 (ss_atoms에 추가)
    → Coil segment  → CA 원자 그대로 통과 (screen_atoms에 추가)
  * 문제: Helix/Sheet segment의 CA 원자도 coil로 처리되어 screen_atoms에 추가됨

수정 흐름:
  전체 chain atoms 순회
    → Helix segment → 나선 지오메트리만 ss_atoms에 추가
                       해당 CA 원자는 screen_atoms에 추가 안 함
    → Sheet segment → 리본 지오메트리만 ss_atoms에 추가
                       해당 CA 원자는 screen_atoms에 추가 안 함
    → Coil segment  → CA 원자 그대로 screen_atoms에 추가
```

구체적으로: segment를 먼저 분류(helix/sheet/coil 구간 목록 생성)한 뒤, coil 구간 원자만 coil 처리 대상에 포함시킨다. helix/sheet 구간 원자는 지오메트리 계산에만 사용하고 coil 파이프라인에서 제외한다.

### 수정 대상: `src/visualization/Screen.cpp`

`project()` 또는 `draw_line()` 함수에서 coil 포인트 렌더링 방식 변경:

```
현재: coil 포인트마다 크로스(+) 패턴으로 5픽셀 이상 그림
수정: coil 포인트는 단일 픽셀만 그림 (또는 2픽셀 크로스 최소화)
```

조건: `if (point.structure == 'x')` → 픽셀 크기 축소

---

## 기능 2: pLDDT 색상 표시

### 데이터 획득: `src/structure/Protein.cpp`

`load_init_atoms()` 내 CA 원자 로딩 부분에서 Gemmi B-factor 읽기:
```cpp
a.bfactor = atom.b_iso;  // AlphaFold 구조에서 pLDDT 값
a.residue_number = residue.seqid.num.value;
a.residue_name = residue.name;
```

`residue_number`와 `residue_name`은 기능 6에도 공통 사용.

### 색상 적용: `src/visualization/Screen.cpp`

`assign_colors_to_points()` 함수에 `"plddt"` 모드 분기 추가:

```cpp
if (mode == "plddt") {
    float plddt = point.bfactor;
    if (plddt >= 90)      color_id = 71;  // Very high: 파란색
    else if (plddt >= 70) color_id = 72;  // Confident: 청록색
    else if (plddt >= 50) color_id = 73;  // Low: 노란색
    else                  color_id = 74;  // Very low: 주황색
}
```

### CLI

```
structty protein.pdb -m plddt
```

---

## 기능 1: Interface Region 색상 표시

### Interface 정의
두 체인(또는 두 구조) 간 CA-CA 거리 < 8Å인 잔기들을 interface residue로 표시.

### 계산: `src/structure/Protein.cpp`

신규 멤버 함수 `compute_interface(chain_A, chain_B, threshold=8.0)` 추가:

```
입력: 두 체인 ID (또는 두 Protein 객체)
처리:
  chain_A의 각 CA 원자에 대해
    chain_B의 모든 CA 원자와 거리 계산
    최소 거리 < threshold이면 Atom.is_interface = true
  chain_B도 동일하게 처리
```

`Protein.hpp`에 함수 선언 추가.

### CLI

```
structty complex.pdb -m interface
```

single PDB 파일 내 두 체인 사이의 interface를 자동 계산. 체인이 2개 이상일 경우 인접한 체인쌍 전체 계산.

체인이 1개일 경우: 두 번째 PDB 파일을 `-p2 <file>` 인수로 받아 계산 (향후 확장, 현재 구현에서는 단일 PDB 내 다체인만 지원).

### 색상 적용: `src/visualization/Screen.cpp`

```cpp
if (mode == "interface") {
    color_id = point.is_interface ? 43 : 44;  // 강조 or dim
}
```

---

## 기능 4: UTMatrix 정렬 구조 색상 표시

### 전제 조건
`-ut` 인수로 변환행렬이 적용된 상태에서만 동작.

### 계산: `src/structure/Protein.cpp`

신규 멤버 함수 `compute_aligned_regions(other_protein, threshold=5.0)` 추가:

```
입력: UT 변환이 이미 적용된 두 Protein의 screen_atoms
처리:
  protein1의 각 CA 원자에 대해
    protein2의 대응 잔기 CA 원자와 거리 계산 (인덱스 기반 대응)
    거리 < threshold이면 양쪽 모두 Atom.is_aligned = true
```

거리 임계값 5Å는 TM-score 계산 기준과 동일.

### 색상 적용: `src/visualization/Screen.cpp`

```cpp
if (mode == "aligned") {
    color_id = point.is_aligned ? 45 : 46;  // 정렬됨: 밝은 초록 / 비정렬: dim
}
```

### CLI

```
structty query.pdb target.pdb -ut matrix.tsv -m aligned
```

---

## 기능 5: MSA Conservation Score 색상 표시

### 신규 파일: `src/structure/MSAParser.hpp` / `MSAParser.cpp`

**MSAParser 클래스 설계:**

```cpp
class MSAParser {
public:
    // 파일 로드: FASTA 다중정렬 형식 (.aln, .fasta, .fa) 지원
    // A3M 형식도 지원 (소문자 처리)
    bool load(const std::string& filepath);

    // 위치별 Shannon entropy 계산
    // 반환: 잔기 인덱스 → conservation score (0.0=가변, 1.0=보존)
    std::vector<float> compute_conservation();

private:
    std::vector<std::string> sequence_names;
    std::vector<std::string> sequences;  // 정렬된 시퀀스 (갭 포함)
    std::string query_sequence;          // 첫 번째 시퀀스 = 쿼리

    // MSA 위치 → 쿼리 잔기 인덱스 매핑
    std::vector<int> msa_pos_to_query_idx;
};
```

**파일 형식 판별:**
- 첫 문자 `>` → FASTA/A3M
- `# STOCKHOLM` 헤더 → Stockholm (향후 지원, 현재는 에러 메시지)

**A3M 처리:**
- 소문자 잔기(삽입)는 쿼리 기준 갭으로 처리
- 카운팅에는 포함하되 쿼리 인덱스 매핑에서 제외

**Shannon Entropy 계산:**
```
각 MSA 위치 i에 대해:
  20가지 아미노산 + 갭('-')의 빈도 계산
  갭은 엔트로피 계산에서 제외 (effective sequences 기준)
  H(i) = -Σ_aa [f(aa,i) * log2(f(aa,i))]
  conservation(i) = 1 - H(i) / log2(20)
```

### 적용: `src/structure/Protein.cpp`

`apply_conservation_scores(const std::vector<float>& scores)` 함수 추가:
- MSA 쿼리 잔기 인덱스와 PDB 잔기 번호 매핑 후 `Atom.conservation_score` 설정
- 매핑 방법: 0-based 인덱스 순서 기준 (잔기 번호가 불연속적일 수 있으므로 순서 기반)

### 색상 적용: `src/visualization/Screen.cpp`

```cpp
if (mode == "conservation") {
    float score = point.conservation_score;
    if (score < 0) {
        color_id = 11;  // 미설정: dim
    } else {
        int idx = (int)(score * 9.0f);  // 0-9
        idx = std::max(0, std::min(9, idx));
        color_id = 75 + idx;  // 75(가변)~84(보존)
    }
}
```

색상 스펙트럼: 파란색(75, 가변) → 시안 → 녹색 → 노란색 → 빨간색(84, 보존)

### CLI

```
structty protein.pdb --msa alignment.fasta -m conservation
```

---

## 기능 6: 커서(마우스) 입력에 따라 Panel에 잔기 정보 표시

### ncurses 마우스 활성화: `src/visualization/Screen.cpp`

`init()` 또는 생성자에서:
```cpp
mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, nullptr);
mouseinterval(0);
```

### 입력 처리: `Screen::handle_input()`

```cpp
case KEY_MOUSE: {
    MEVENT event;
    if (getmouse(&event) == OK) {
        int mx = event.x;
        int my = event.y;
        update_hover_info(mx, my);
    }
    break;
}
```

### hover 정보 검색: `Screen::update_hover_info(int mx, int my)`

```
screenPixels에서 좌표 (mx, my)에 해당하는 RenderPoint 검색
발견 시: Panel에 잔기 정보 전달
미발견 시: Panel 잔기 정보 섹션 클리어
```

현재 `screenPixels`는 z-buffer 정렬된 2D 배열이므로 O(1) 접근 가능.

### Panel 확장: `src/visualization/Panel.hpp` / `Panel.cpp`

패널 하단에 "Residue Info" 섹션 추가:

```
┌─────────────────────┐
│ [기존 패널 내용]     │
├─────────────────────┤
│ Residue Info        │
│ Chain: A            │
│ Res:   GLU 42       │
│ SS:    Helix        │
│ pLDDT: 87.3         │  (plddt 모드일 때만)
│ Cons:  0.82         │  (conservation 모드일 때만)
└─────────────────────┘
```

`Panel::set_hover_residue(chain, res_num, res_name, ss_type, bfactor, conservation)` 함수 추가.

현재 활성화된 색상 모드에 따라 표시 항목 조절.

### RenderPoint에 잔기 정보 전파

`Screen::project()` 함수에서 RenderPoint 생성 시 `Atom.residue_number`, `Atom.residue_name`, `Atom.bfactor`, `Atom.conservation_score`를 RenderPoint로 복사.

---

## 기능 3: Foldseek 결과 파일 실시간 Hit 탐색

### 데이터 구조: `src/structure/FoldseekParser.hpp` (신규)

```cpp
struct FoldseekHit {
    std::string query;
    std::string target;
    float fident;
    int alnlen;
    float evalue;
    float bitscore;
    float lddt;          // --format-output에 lddt 포함 시
    float qtmscore;
    float ttmscore;
    // U matrix (3x3) 및 T vector: 있으면 구조 정렬 가능
    float U[3][3];
    float T[3];
    bool has_transform;  // U,T 컬럼 존재 여부
    std::string qaln;    // 쿼리 정렬 문자열 (기능 4 확장용)
    std::string taln;    // 타겟 정렬 문자열
};

class FoldseekParser {
public:
    bool load(const std::string& filepath);
    const std::vector<FoldseekHit>& get_hits() const;
    int hit_count() const;

private:
    std::vector<FoldseekHit> hits;
    void parse_line(const std::string& line);
    // 컬럼 헤더 자동 감지 (첫 줄이 '#query ...' 형식인지 확인)
    std::vector<std::string> column_headers;
};
```

### Foldseek 파일 형식 지원

두 가지 형식 자동 감지:
1. **기본 m8 형식** (12컬럼): `query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits`
2. **확장 형식**: `--format-output`으로 `lddt,qtmscore,ttmscore,u,t` 등이 추가된 경우

첫 줄이 `#`으로 시작하면 헤더로 처리하여 컬럼 이름 파싱. 헤더 없으면 기본 m8 형식으로 가정.

U,T 컬럼이 있으면 `has_transform = true`로 설정, 구조 정렬 적용 가능.

### Screen 통합: `src/visualization/Screen.hpp/cpp`

**멤버 변수 추가:**
```cpp
std::vector<FoldseekHit> foldseek_hits;
int current_hit_idx = -1;
std::string db_path;  // target PDB 파일 경로
```

**키 처리 (`handle_input()`):**
```cpp
case 'n': case 'N':
    if (!foldseek_hits.empty())
        load_next_hit(+1);
    break;
case 'p': case 'P':
    if (!foldseek_hits.empty())
        load_next_hit(-1);
    break;
```

**`Screen::load_next_hit(int delta)` 함수:**
```
1. current_hit_idx += delta (범위 클램핑)
2. FoldseekHit& hit = foldseek_hits[current_hit_idx]
3. target PDB 파일 경로 결정:
   - db_path + "/" + hit.target + ".pdb" 시도
   - db_path + "/" + hit.target + ".cif" 시도
   - 파일 없으면 패널에 "File not found" 표시 후 skip
4. 두 번째 Protein 객체 로드 (target)
5. hit.has_transform이면 set_utmatrix(hit.U, hit.T) 적용
6. 재렌더링 트리거
7. 패널 hit 정보 업데이트
```

**Panel hit 정보 표시:**
```
┌─────────────────────┐
│ Foldseek Hits       │
│ [3 / 47]            │
│ Target: 2xyz_B      │
│ E-val:  1.2e-15     │
│ TM:     0.87        │
│ lDDT:   0.74        │
│ [N]ext  [P]rev      │
└─────────────────────┘
```

### target PDB 파일 경로 결정 로직

Foldseek target ID 형식: `2xyz_B` (PDB ID + chain) 또는 `AF-P12345-F1` (AlphaFold) 등 다양.

지원할 경로 패턴 (순서대로 시도):
```
{db_path}/{target}.pdb
{db_path}/{target}.cif
{db_path}/{target_lower}.pdb
{db_path}/{pdb_id}/{pdb_id}.pdb  (일부 DB 구조)
```

파일 미발견 시: 패널에 경고 표시, 해당 hit 건너뜀.

---

## 파일별 수정 상세 목록

| 파일 | 수정 유형 | 주요 변경 내용 |
|---|---|---|
| `src/structure/Atom.hpp` | 수정 | `bfactor`, `is_interface`, `is_aligned`, `conservation_score`, `residue_number`, `residue_name` 필드 추가 |
| `src/visualization/RenderPoint.hpp` | 수정 | `residue_number`, `residue_name`, `bfactor`, `conservation_score` 필드 추가 |
| `src/visualization/Palette.hpp` | 수정 | color pair 43-46, 71-74, 75-84 추가 |
| `src/structure/Parameters.hpp` | 수정 | `-m plddt/interface/conservation/aligned`, `--msa`, `-fs`, `--db-path` 인수 추가 |
| `src/structure/Protein.hpp` | 수정 | `compute_interface()`, `compute_aligned_regions()`, `apply_conservation_scores()` 선언 |
| `src/structure/Protein.cpp` | 수정 | B-factor/residue 정보 로딩, interface/aligned 계산 구현 |
| `src/structure/StructureMaker.cpp` | 수정 | Helix/Sheet 구간 CA 원자 coil 파이프라인 제외, coil 굵기 축소 |
| `src/visualization/Screen.hpp` | 수정 | foldseek_hits, current_hit_idx, db_path, hover 관련 멤버 추가 |
| `src/visualization/Screen.cpp` | 수정 | `assign_colors_to_points()` 신규 모드 추가, 마우스 입력 처리, `load_next_hit()`, `update_hover_info()` 추가, `project()`에서 잔기 정보 전파 |
| `src/visualization/Panel.hpp` | 수정 | `set_hover_residue()`, `set_foldseek_hit_info()` 선언 |
| `src/visualization/Panel.cpp` | 수정 | 잔기 정보 섹션, hit 정보 섹션 렌더링 추가 |
| `src/structure/MSAParser.hpp` | **신규** | MSAParser 클래스 선언 |
| `src/structure/MSAParser.cpp` | **신규** | FASTA/A3M 파싱, Shannon entropy 계산 구현 |
| `src/structure/FoldseekParser.hpp` | **신규** | FoldseekHit 구조체, FoldseekParser 클래스 선언 |
| `src/structure/FoldseekParser.cpp` | **신규** | m8/확장 형식 파싱 구현 |
| `src/structty.cpp` | 수정 | 신규 인수 처리, MSAParser/FoldseekParser 초기화 연결 |

---

## CLI 사용 예시 (완성 후)

```bash
# 기능 2: pLDDT 색상
structty af_protein.pdb -m plddt

# 기능 1: Interface region (다체인 복합체)
structty complex.pdb -m interface

# 기능 4: 구조 정렬 색상
structty query.pdb target.pdb -ut matrix.tsv -m aligned

# 기능 5: Conservation 색상
structty protein.pdb --msa alignment.fasta -m conservation

# 기능 3: Foldseek hit 탐색 (N/P 키로 이동)
structty query.pdb -fs search_result.tsv --db-path /data/pdb

# 기능 6: 마우스 hover (항상 활성화, 패널 필요)
structty protein.pdb  # 마우스 이동 시 자동으로 잔기 정보 표시

# 복합 사용
structty query.pdb -fs result.tsv --db-path /data/pdb -m aligned -s
```

---

## 키 매핑 정리 (추가분)

| 키 | 기능 |
|---|---|
| `N` / `n` | 다음 Foldseek hit 로드 (기능 3) |
| `P` / `p` | 이전 Foldseek hit 로드 (기능 3) |
| 마우스 이동 | 커서 위치 잔기 정보 패널 표시 (기능 6) |

---

## 구현 시 주의사항

1. **기존 기능 무결성**: 기존 `-m protein/chain/rainbow` 모드와 `-s` (이차구조 표시)는 변경 없이 동작해야 함
2. **Foldseek hit 전환 시 메모리**: 매 hit 전환마다 target Protein 객체를 새로 생성하므로 이전 객체 명시적 해제 필요
3. **마우스 지원 여부**: 일부 터미널에서 마우스 이벤트 미지원. `mousemask()` 반환값으로 지원 여부 확인 후, 미지원 시 패널에 "Mouse not supported" 표시
4. **MSA 잔기 인덱스 매핑**: PDB 잔기 번호가 불연속적(insertion code 등)일 수 있으므로 순서(0-based index) 기준 매핑 사용
5. **pLDDT vs 일반 B-factor**: `-m plddt`는 B-factor 컬럼 값을 pLDDT로 해석하므로, 일반 X-ray 구조에 사용 시 의미 없음을 패널에 표시할 것
6. **coil 굵기 수정**: StructureMaker 수정과 Screen.cpp draw_line 수정 두 곳 모두 필요. 한쪽만 수정하면 불완전함
