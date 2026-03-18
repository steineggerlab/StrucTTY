# StrucTTY 기능 추가 구현 계획

## 개요

7가지 신규 기능을 순서대로 구현한다. 각 기능은 독립적으로 동작 가능하도록 설계하되, 공통 인프라(Atom 필드 확장, Palette 확장)는 먼저 처리한다.
이때, 기본적으로 제공되지 않는 다른 패키지나 프로젝트를 활용하지 않도록 한다. 단, `curl`이나 `wget` 등 시스템에 이미 설치된 CLI 툴을 `popen()`으로 호출하는 것은 허용한다.
또한 코드의 병목 및 비효율적 병목 구조를 피하여 빠른 실행시간을 보장하도록 한다.

---

## 구현 순서

1. **공통 인프라** (Atom.hpp, RenderPoint.hpp, Palette.hpp 확장) - 모든 기능의 전제 조건
2. **기능 7** - 코일 렌더링 버그 수정
3. **기능 2** - pLDDT 색상 표시
4. **기능 1** - Interface region 색상 표시
5. **기능 4** - UTMatrix 정렬 구조 색상 표시
6. **기능 5** - MSA Conservation Score 색상 표시
7. **기능 6** - 커서(마우스) 기반 잔기 정보 패널 표시
8. **기능 3** - Foldseek 결과 파일 실시간 Hit 탐색

---

## 단계 0: 공통 인프라 ✅

### 0-1. `src/structure/Atom.hpp` 필드 확장 ✅

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

### 0-2. `src/visualization/RenderPoint.hpp` 필드 확장 ✅

**색상 적용 코드(`assign_colors_to_points()`)는 Atom이 아닌 RenderPoint를 참조한다.
따라서 색상 판단에 필요한 모든 플래그와 값은 반드시 RenderPoint에도 존재해야 한다.**

```cpp
struct RenderPoint {
    // 기존 필드 유지
    int x, y;
    float depth;
    char pixel;
    int color_id;
    std::string chainID;
    char structure;

    // 기능 1: interface region
    bool is_interface = false;

    // 기능 4: 정렬 구조
    bool is_aligned = false;

    // 기능 2: pLDDT
    float bfactor = 0.0f;

    // 기능 5: conservation
    float conservation_score = -1.0f;

    // 기능 6: 잔기 정보
    int residue_number = -1;
    std::string residue_name = "";
};
```

### 0-3. `src/visualization/Screen.cpp` — `project()` 함수 내 RenderPoint 전파 목록 ✅

`project()` 함수에서 Atom → RenderPoint로 값을 복사할 때, **아래 필드 전체**를 빠짐없이 복사해야 한다:

```cpp
rp.bfactor            = atom.bfactor;
rp.is_interface       = atom.is_interface;
rp.is_aligned         = atom.is_aligned;
rp.conservation_score = atom.conservation_score;
rp.residue_number     = atom.residue_number;
rp.residue_name       = atom.residue_name;
```

이 복사가 누락되면 해당 기능의 색상 분기가 RenderPoint 기본값(false, 0.0f, -1)만 읽게 되어 전체가 단일 색상으로 표시된다. **구현 후 반드시 각 필드가 복사되는지 코드에서 확인한다.**

### 0-4. `src/visualization/Palette.hpp` 색상 확장 ✅

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

### 0-5. `src/structure/Parameters.hpp` CLI 인수 확장 ✅

신규 인수:
```
-m plddt          기능 2: pLDDT 색상 모드
-m interface      기능 1: interface region 색상 모드
-m conservation   기능 5: MSA conservation 색상 모드
-m aligned        기능 4: 정렬 구조 색상 모드
--msa <file>      기능 5: MSA 파일 경로
-fs <file>        기능 3: Foldseek 결과 파일 경로
--db-path <dir>   기능 3: target 구조 파일 로컬 디렉토리 (선택)
```

기존 `-m` 옵션의 유효값 목록에 위 4가지 추가.

---

## 기능 7: 헬릭스/시트 영역의 Coil 제거 + Coil 굵기 조정 ✅

### 문제 원인
`StructureMaker::calculate_ss_points()`에서 helix/sheet 지오메트리를 생성하지만, 해당 잔기들의 원본 CA 원자들도 coil 리스트에 그대로 남아 있어 이중으로 그려진다.

### 수정 대상: `src/structure/StructureMaker.cpp` ✅

**Helix/Sheet CA 원자를 coil 파이프라인에서 제외하되, 접합부(junction) 처리를 반드시 포함한다.**

```
수정 흐름:
  1. 전체 chain atoms를 순회하여 각 잔기의 SS 유형(helix/sheet/coil)을 먼저 분류한다.
  2. coil 구간 원자만 screen_atoms에 추가한다.
  3. helix/sheet 구간 원자는 지오메트리 계산에만 사용하고 coil 파이프라인에서 제외한다.

접합부 처리 (끊김 방지):
  helix/sheet segment의 첫 번째 잔기와 마지막 잔기의 CA 원자는
  coil 파이프라인에도 동시에 추가한다.
  이렇게 하면 coil ↔ helix/sheet 경계에서 선이 끊기지 않는다.
```

즉, helix/sheet segment 경계 잔기(첫 번째, 마지막)는 이중 등록(ss_atoms + screen_atoms)이 허용된다.

### 수정 대상: `src/visualization/Screen.cpp` — coil 굵기 ✅

`project()` 또는 `draw_line()` 함수에서 coil 포인트 렌더링 방식 변경:

```
현재(이전): coil 포인트마다 크로스(+) 패턴으로 5픽셀 이상 그림  ← 너무 두꺼움
이번 목표:  coil 포인트를 3픽셀 크로스(중심 + 상하 또는 좌우 1픽셀씩)로 그림
            단일 픽셀은 너무 가늘어 구조가 보이지 않으므로 사용하지 않는다.
```

구체적으로: `if (point.structure == 'x')` 분기에서 중심 픽셀 + 인접 2픽셀(총 3픽셀) 패턴을 사용한다. 5픽셀 크로스와 단일 픽셀의 중간값이다.

---

## 기능 2: pLDDT 색상 표시 ✅

### 데이터 획득: `src/structure/Protein.cpp` ✅

`load_init_atoms()` 내 CA 원자 로딩 부분에서 Gemmi B-factor 읽기:
```cpp
a.bfactor = atom.b_iso;  // AlphaFold 구조에서 pLDDT 값
a.residue_number = residue.seqid.num.value;
a.residue_name = residue.name;
```

`residue_number`와 `residue_name`은 기능 6에도 공통 사용.

### 색상 적용: `src/visualization/Screen.cpp` ✅

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

### [추가] `src/structure/StructureMaker.cpp` — geometry 포인트 메타데이터 보간 ✅

plan.md에 명시되지 않았으나 구현 중 발견된 암묵적 요구사항:
`calculate_ss_points()`에서 H/S geometry 포인트를 `Atom(x, y, z, 'H'/'S')`로 생성하면 `bfactor=0.0f`가 됩니다.
`project()`의 `before_draw` 패턴은 `screen_atoms[i].bfactor`를 그대로 RenderPoint에 복사하므로,
geometry 포인트의 bfactor가 0이면 pLDDT 색상이 항상 Very Low(74, 주황색)로 표시됩니다.

**수정 내용**: geometry Atom 생성 시 세그먼트 양 끝 잔기의 값을 `t`(0→1)로 선형 보간하여 설정.
- helix spiral: `segment.front()` ~ `segment.back()` bfactor/conservation_score 보간
- sheet ribbon: 각 Ca-Ca pair의 `pa` ~ `pb` bfactor/conservation_score 보간
- is_interface, is_aligned: pa(또는 front) 값 사용 (보간 없이)

### CLI

```
structty protein.pdb -m plddt
```

---

## 기능 1: Interface Region 색상 표시 ✅

### Interface 정의
두 체인(또는 두 구조) 간 CA-CA 거리 < 8Å인 잔기들을 interface residue로 표시.

### 계산: `src/structure/Protein.cpp` ✅

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

**주의: `Atom.is_interface`를 설정한 뒤, 해당 값이 `project()` 과정에서 RenderPoint로 복사되는지 반드시 확인한다 (0-3절 참조).**

### CLI ✅

```
structty complex.pdb -m interface
```

single PDB 파일 내 두 체인 사이의 interface를 자동 계산. 체인이 2개 이상일 경우 인접한 체인쌍 전체 계산.

체인이 1개일 경우: 두 번째 PDB 파일을 `-p2 <file>` 인수로 받아 계산 (향후 확장, 현재 구현에서는 단일 PDB 내 다체인만 지원).

### [추가] 구현 세부사항

- `compute_interface_pair(chain_A, chain_B, threshold)` private 헬퍼: 두 체인 쌍에 대한 CA-CA 거리 계산
- `compute_interface(threshold=8.0f)` public 진입점: 모든 체인 쌍 자동 계산 후 `sync_interface_to_screen()` 호출
- `sync_interface_to_screen()` private: init_atoms → screen_atoms 전파
  - residue_number >= 0 (coil/junction): residue_number로 직접 매핑
  - residue_number == -1 (H/S geometry): 3D 최근접 init_atom 값 사용
- `Screen::compute_interface_all()` 추가: 로드된 모든 Protein에 대해 호출
- `structty.cpp`: `-m interface` 시 `screen.compute_interface_all()` 호출

### 색상 적용: `src/visualization/Screen.cpp` ✅

```cpp
if (mode == "interface") {
    color_id = point.is_interface ? 43 : 44;  // 강조 or dim
}
```

---

## 기능 4: UTMatrix 정렬 구조 색상 표시 ✅

### 전제 조건 및 `-ut` 인수 처리 방침

**`-fs`가 있는 경우: `-ut` 인수를 받지 않는다.**
Foldseek hit에 이미 U matrix와 T vector가 포함되어 있다(`has_transform`). hit을 로드할 때 해당 hit의 transform을 자동으로 적용한다. 사용자가 별도로 `-ut` 파일을 제공할 필요가 없으며, 제공하더라도 무시한다.

**`-fs`가 없는 경우: `-ut` 인수로 변환행렬 파일을 받는다.**
이 경우 두 번째 PDB 파일과 함께 사용한다.

### 잔기 대응 방식

**인덱스 기반 대응(i번째 잔기끼리 비교)은 사용하지 않는다.** 두 단백질의 서열 길이가 다르거나 alignment gap이 있으면 같은 인덱스끼리의 거리 비교는 의미가 없다.

**`-fs`가 있는 경우 (alignment string 기반, 권장):**

```
FoldseekHit.qaln과 FoldseekHit.taln을 동시에 순회:
  qaln[i] != '-' and taln[i] != '-' 인 위치만 대응 쌍으로 수집
  qaln 기준 잔기 인덱스 → protein1의 CA 원자
  taln 기준 잔기 인덱스 → protein2의 CA 원자
  hit의 U/T 변환 후 3D 거리 < 10.0이면 양쪽 모두 Atom.is_aligned = true
```

**`-fs`가 없는 경우 (nearest-neighbor fallback):**

```
UT 변환 후 nearest-neighbor 방식으로 대응:
  protein1의 각 CA 원자에 대해
    protein2의 모든 CA 원자 중 최소 거리 탐색
    최소 거리 < 10.0Å이면 양쪽 모두 is_aligned = true
```

패널에 현재 사용 중인 대응 방식을 표시한다:
- `-fs` 있음: `Align: aln-string`
- `-fs` 없음: `Align: nearest-nbr`

### 계산: `src/structure/Protein.cpp`

신규 멤버 함수 두 가지:

```cpp
// Foldseek alignment string 기반
void compute_aligned_regions_from_aln(
    Protein& other,
    const std::string& qaln,
    const std::string& taln,
    float threshold = 5.0f
);

// nearest-neighbor 기반 (alignment string 없을 때 fallback)
void compute_aligned_regions_nn(
    Protein& other,
    float threshold = 10.0f
);
```

**주의: `Atom.is_aligned`를 설정한 뒤, 해당 값이 `project()` 과정에서 RenderPoint로 복사되는지 반드시 확인한다 (0-3절 참조).**

### 색상 적용: `src/visualization/Screen.cpp` ✅

```cpp
if (mode == "aligned") {
    color_id = point.is_aligned ? 45 : 46;  // 정렬됨: 밝은 초록 / 비정렬: dim
}
```

### CLI ✅

```bash
# -fs와 함께 사용 (alignment string 기반, -ut 불필요) — 기능 3 구현 시 연결
structty query.pdb -fs result.m8 -m aligned

# -ut 단독 사용 (nearest-neighbor fallback, threshold=10Å)
structty query.pdb target.pdb -ut matrix.tsv -m aligned
```

### [추가] 구현 세부사항

- `apply_ut_to_init_atoms(const float* U, const float* T)` public: UT transform을 `init_atoms`에도 적용.
  `Screen::set_utmatrix()` 에서 `do_naive_rotation()` / `do_shift()` 직후 호출하여 `screen_atoms`와 `init_atoms`의 좌표 공간을 일치시킨다.
- `compute_aligned_regions_nn(Protein& other, threshold=10.0f)` public: nearest-neighbor 기반 정렬 잔기 표시.
  UT transform 적용 후 `init_atoms` 좌표 공간에서 CA-CA 거리를 계산한다. `-ut` 단독 사용 시 진입점.
- `compute_aligned_regions_from_aln(Protein& other, qaln, taln, threshold=10.0f)` public: alignment string 기반 정렬 잔기 표시.
  함수 자체는 기능 4에서 구현 완료. **실제 호출은 기능 3 (`-fs`) 구현 시 `Screen::load_next_hit()` 에서 연결한다.**
  (`hit.has_aln == true` 이면 `compute_aligned_regions_from_aln()` 호출, 아니면 `compute_aligned_regions_nn()` 호출 — plan.md 기능 3, `load_next_hit()` 7번 항목 참조)
- `sync_aligned_to_screen()` private: `init_atoms.is_aligned` → `screen_atoms` 전파.
  `sync_bool_field_to_screen()` 공통 헬퍼(pointer-to-member)로 구현. H/S geometry 원자는 index 기준 최근접 coil 원자 값 사용.
- `Screen::compute_aligned_all(threshold=10.0f)` 추가: 로드된 모든 Protein 쌍에 대해 NN 계산 후 패널에 `"nearest-nbr"` 표시.
- `Panel::set_align_method(method)` 추가: `panel_mode == "aligned"` 일 때 구분선 아래에 `Align: nearest-nbr` 또는 `Align: aln-string` 표시.
- `structty.cpp`: `-m aligned` 시 `screen.compute_aligned_all()` 호출.

---

## 기능 5: MSA Conservation Score 색상 표시

### 신규 파일: `src/structure/MSAParser.hpp` / `MSAParser.cpp`

**MSAParser 클래스 설계:**

```cpp
class MSAParser {
public:
    bool load(const std::string& filepath);
    std::vector<float> compute_conservation();

private:
    std::vector<std::string> sequence_names;
    std::vector<std::string> sequences;
    std::string query_sequence;
    std::vector<int> msa_pos_to_query_idx;
};
```

**파일 형식 판별:**
- 첫 문자 `>` → FASTA/A3M
- `# STOCKHOLM` 헤더 → Stockholm (향후 지원, 현재는 에러 메시지)

**A3M 처리:**
- 소문자 잔기(삽입)는 쿼리 기준 갭으로 처리

**Shannon Entropy 계산:**
```
각 MSA 위치 i에 대해:
  갭 제외 빈도 계산
  H(i) = -Σ_aa [f(aa,i) * log2(f(aa,i))]
  conservation(i) = 1 - H(i) / log2(20)
```

### 적용: `src/structure/Protein.cpp`

`apply_conservation_scores(const std::vector<float>& scores)` 함수 추가.
매핑: 0-based 인덱스 순서 기준.

### 색상 적용: `src/visualization/Screen.cpp`

```cpp
if (mode == "conservation") {
    float score = point.conservation_score;
    if (score < 0) {
        color_id = 11;
    } else {
        int idx = std::max(0, std::min(9, (int)(score * 9.0f)));
        color_id = 75 + idx;
    }
}
```

### CLI

```
structty protein.pdb --msa alignment.fasta -m conservation
```

---

## 기능 6: 커서(마우스) 기반 잔기 정보 패널 표시

### ncurses 마우스 활성화: `src/visualization/Screen.cpp`

`init()` 또는 생성자에서:
```cpp
mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, nullptr);
mouseinterval(0);
```

### 이벤트 처리 방식

**마우스 이동(hover) 이벤트만 사용한다. 클릭 이벤트는 잔기 정보 표시에 사용하지 않는다.**

ncurses에서 `REPORT_MOUSE_POSITION`을 활성화하면 마우스 이동 시 `KEY_MOUSE` 이벤트가 발생한다. 단, 터미널에 따라 이동 이벤트 지원 여부가 다르다.

```cpp
case KEY_MOUSE: {
    MEVENT event;
    if (getmouse(&event) == OK) {
        // 버튼 이벤트와 이동 이벤트 모두 처리
        update_hover_info(event.x, event.y);
    }
    break;
}
```

`mousemask()` 반환값이 0이면 마우스 미지원 터미널이므로 패널의 Residue Info 섹션에 "Mouse not supported"를 고정 표시한다.

### 패널 Residue Info 섹션 — 고정 크기 유지

**패널의 크기는 잔기 정보 유무에 관계없이 항상 동일해야 한다.**

Residue Info 섹션은 항상 고정된 줄 수를 차지한다. 정보가 없을 때는 해당 줄을 빈 칸으로 채운다.

```
┌─────────────────────┐
│ [기존 패널 내용]     │
├─────────────────────┤
│ Residue Info        │   ← 항상 표시
│ Chain: A            │   ← 없으면 "Chain: -"
│ Res:   GLU 42       │   ← 없으면 "Res:   -"
│ SS:    Helix        │   ← 없으면 "SS:    -"
│ pLDDT: 87.3         │   ← plddt 모드일 때만, 없으면 "pLDDT: -"
│ Cons:  0.82         │   ← conservation 모드일 때만, 없으면 "Cons:  -"
└─────────────────────┘
```

섹션 높이는 활성화된 모드에 따라 컴파일 시점에 결정되며, 런타임에 줄 수가 바뀌지 않는다.

`Panel::set_hover_residue()` 호출 시 항상 모든 줄을 다시 그린다(빈 칸 포함). 패널 전체를 다시 그리지 않고 Residue Info 섹션만 갱신하여 화면 깜빡임을 최소화한다.

---

## 기능 3: Foldseek 결과 파일 실시간 Hit 탐색

### 설계 방침

**`--db-path`는 선택 인수다.** 로컬 DB 디렉토리를 제공하지 않으면, 프로그램이 Foldseek 결과 파일의 상위 N개 hit을 자동으로 다운로드한다.

자동 다운로드는 시스템에 설치된 `curl` 또는 `wget`을 `popen()`으로 호출하여 수행한다. 별도 라이브러리를 추가하지 않는다.

### 입력 파일 형식 명세

**Foldseek 실행 명령어 (사용자가 미리 실행해야 함):**

```bash
foldseek easy-search query.pdb /path/to/db result.m8 tmp \
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddt,qtmscore,ttmscore,u,t,qaln,taln"
```

**컬럼 구성 (탭 구분, 헤더 없음):**

| 인덱스 | 컬럼명 | 설명 |
|--------|--------|------|
| 0 | query | 쿼리 파일명 또는 ID |
| 1 | target | 타겟 ID (예: `2xyz_B`, `AF-P12345-F1-model_v4`) |
| 2 | fident | 서열 동일성 (0.0~1.0) |
| 3 | alnlen | alignment 길이 |
| 4 | mismatch | mismatch 수 |
| 5 | gapopen | gap open 수 |
| 6 | qstart | 쿼리 시작 잔기 |
| 7 | qend | 쿼리 끝 잔기 |
| 8 | tstart | 타겟 시작 잔기 |
| 9 | tend | 타겟 끝 잔기 |
| 10 | evalue | E-value |
| 11 | bits | bitscore |
| 12 | lddt | lDDT 점수 |
| 13 | qtmscore | query 기준 TM-score |
| 14 | ttmscore | target 기준 TM-score |
| 15~23 | u | 3x3 rotation matrix (9개 값, 탭 구분) |
| 24~26 | t | translation vector (3개 값, 탭 구분) |
| 27 | qaln | 쿼리 alignment 문자열 |
| 28 | taln | 타겟 alignment 문자열 |

**U matrix와 T vector 파싱 주의사항:**
- `u` 필드는 9개의 탭 구분 값으로 확장된다 (단일 컬럼이 아님)
- `t` 필드는 3개의 탭 구분 값으로 확장된다
- 전체 컬럼 수: 29개

**컬럼 수 기반 자동 포맷 감지:**

```
컬럼 수 == 29: 전체 형식 (U, T, qaln, taln 포함)
컬럼 수 == 17: 축소 형식 (qaln, taln 포함, U/T 없음)
컬럼 수 == 12: 기본 m8 형식 (qaln, taln, U, T 모두 없음)
그 외: 파싱 실패, 에러 메시지 출력
```

헤더 행은 지원하지 않는다. Foldseek의 기본 출력에는 헤더가 없다.

### 자동 다운로드 로직

`--db-path`가 제공되지 않으면, 프로그램 시작 시 결과 파일의 상위 10개 hit을 자동 다운로드한다.

다운로드 디렉토리: `~/.cache/structty/pdb/` (없으면 생성)

**target ID 형식 분류 및 다운로드 URL:**

```
패턴 1: 숫자+문자 4글자 + 언더스코어 + 체인 (예: 2xyz_B)
  → PDB ID = "2xyz", chain = "B"
  → URL: https://files.rcsb.org/download/2xyz.pdb
  → 저장: ~/.cache/structty/pdb/2xyz_B.pdb
  → 로드 시 체인 B만 추출

패턴 2: AF- 로 시작 (예: AF-P12345-F1-model_v4)
  → UniProt ID = "P12345"
  → URL: https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.pdb
  → 저장: ~/.cache/structty/pdb/AF-P12345-F1-model_v4.pdb

패턴 3: 그 외
  → 다운로드 시도하지 않음, 패널에 "Unknown target format" 표시
```

**다운로드 실행:**

```cpp
// curl 우선, 없으면 wget
std::string download_cmd =
    "curl -sf -o " + dest_path + " " + url +
    " || wget -q -O " + dest_path + " " + url;
popen(download_cmd.c_str(), "r");
```

다운로드 성공 여부는 파일 존재 및 크기로 확인한다. 실패 시 패널에 "Download failed: {target}" 표시.

**캐시 활용:** `~/.cache/structty/pdb/{target}.pdb`가 이미 존재하면 다운로드를 건너뛴다.

**`--db-path`가 제공된 경우:** 로컬 디렉토리에서만 탐색하고 다운로드는 수행하지 않는다.

```
탐색 순서:
  {db_path}/{target}.pdb
  {db_path}/{target}.cif
  {db_path}/{target_lowercase}.pdb
  {db_path}/{target_lowercase}.cif
```

### 프로그램 시작 시 동작

```
1. -fs 파일을 파싱하여 hit 목록 로드
2. --db-path 없으면: 상위 10개 hit을 백그라운드에서 순차 다운로드
   (다운로드 중에는 패널에 "Downloading [3/10]..." 표시)
3. 첫 번째 hit을 자동으로 로드하여 표시
```

### 데이터 구조: `src/structure/FoldseekParser.hpp` (신규)

```cpp
struct FoldseekHit {
    std::string query;
    std::string target;
    float fident;
    int alnlen;
    float evalue;
    float bits;
    float lddt = -1.0f;
    float qtmscore = -1.0f;
    float ttmscore = -1.0f;
    float U[3][3] = {};
    float T[3] = {};
    bool has_transform = false;
    std::string qaln;
    std::string taln;
    bool has_aln = false;
};

class FoldseekParser {
public:
    bool load(const std::string& filepath);
    const std::vector<FoldseekHit>& get_hits() const;
    int hit_count() const;

private:
    std::vector<FoldseekHit> hits;
    int detected_format = 0;  // 12, 17, 29
    bool parse_line(const std::string& line, FoldseekHit& hit);
};
```

### Screen 통합: `src/visualization/Screen.hpp/cpp`

**멤버 변수 추가:**
```cpp
std::vector<FoldseekHit> foldseek_hits;
int current_hit_idx = -1;
std::string db_path;  // 비어 있으면 자동 다운로드 모드
```

**키 처리:**
```cpp
case 'n': case 'N':
    if (!foldseek_hits.empty()) load_next_hit(+1);
    break;
case 'p': case 'P':
    if (!foldseek_hits.empty()) load_next_hit(-1);
    break;
```

**`Screen::load_next_hit(int delta)` 함수:**
```
1. current_hit_idx += delta (범위 클램핑)
2. FoldseekHit& hit = foldseek_hits[current_hit_idx]
3. 구조 파일 경로 결정:
   - db_path 있음: 로컬 탐색
   - db_path 없음: ~/.cache/structty/pdb/{target}.pdb 확인
4. 파일 없으면 패널에 "File not found: {target}" 표시 후 return
5. 두 번째 Protein 객체 로드
6. hit.has_transform이면 해당 hit의 U, T를 자동 적용 (-ut 파일 불필요)
7. -m aligned이고 hit.has_aln이면 compute_aligned_regions_from_aln() 호출
   hit.has_aln이 false이면 compute_aligned_regions_nn(threshold=10.0f) 호출
8. 재렌더링 트리거
9. 패널 hit 정보 업데이트
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
│ Align:  aln-string  │
│ [N]ext  [P]rev      │
└─────────────────────┘
```

---

## 파일별 수정 상세 목록

| 파일 | 수정 유형 | 주요 변경 내용 |
|---|---|---|
| `src/structure/Atom.hpp` | 수정 | `bfactor`, `is_interface`, `is_aligned`, `conservation_score`, `residue_number`, `residue_name` 필드 추가 |
| `src/visualization/RenderPoint.hpp` | 수정 | `is_interface`, `is_aligned`, `bfactor`, `conservation_score`, `residue_number`, `residue_name` 필드 추가 |
| `src/visualization/Palette.hpp` | 수정 | color pair 43-46, 71-74, 75-84 추가 |
| `src/structure/Parameters.hpp` | 수정 | `-m plddt/interface/conservation/aligned`, `--msa`, `-fs`, `--db-path`(선택) 인수 추가 |
| `src/structure/Protein.hpp` | 수정 | `compute_interface()`, `compute_aligned_regions_from_aln()`, `compute_aligned_regions_nn()`, `apply_conservation_scores()` 선언 |
| `src/structure/Protein.cpp` | 수정 | B-factor/residue 정보 로딩, interface/aligned 계산 구현 |
| `src/structure/StructureMaker.cpp` | 수정 | Helix/Sheet 구간 CA 원자 coil 파이프라인 제외 + **경계 잔기 이중 등록으로 끊김 방지**, coil 3픽셀 크로스로 조정 |
| `src/visualization/Screen.hpp` | 수정 | foldseek_hits, current_hit_idx, db_path, hover 관련 멤버 추가 |
| `src/visualization/Screen.cpp` | 수정 | `project()`에서 모든 Atom 필드를 RenderPoint로 복사, `assign_colors_to_points()` 신규 모드 추가, 마우스 이동 이벤트 처리(hover), `load_next_hit()`, `update_hover_info()` 추가 |
| `src/visualization/Panel.hpp` | 수정 | `set_hover_residue()`, `set_foldseek_hit_info()` 선언, Residue Info 섹션 고정 높이 상수 정의 |
| `src/visualization/Panel.cpp` | 수정 | Residue Info 섹션 **고정 크기** 렌더링, hit 정보 섹션 렌더링 추가 |
| `src/structure/MSAParser.hpp` | **신규** | MSAParser 클래스 선언 |
| `src/structure/MSAParser.cpp` | **신규** | FASTA/A3M 파싱, Shannon entropy 계산 구현 |
| `src/structure/FoldseekParser.hpp` | **신규** | FoldseekHit 구조체, FoldseekParser 클래스 선언 |
| `src/structure/FoldseekParser.cpp` | **신규** | 컬럼 수 기반 포맷 자동 감지, 파싱 구현 |
| `src/structure/PDBDownloader.hpp` | **신규** | target ID 형식 분류, 다운로드 URL 생성, 캐시 관리 |
| `src/structure/PDBDownloader.cpp` | **신규** | curl/wget popen 기반 다운로드 구현 |
| `src/structty.cpp` | 수정 | 신규 인수 처리, MSAParser/FoldseekParser/PDBDownloader 초기화 연결 |

---

## CLI 사용 예시 (완성 후)

```bash
# 기능 2: pLDDT 색상
structty af_protein.pdb -m plddt

# 기능 1: Interface region (다체인 복합체)
structty complex.pdb -m interface

# 기능 4: 구조 정렬 색상 (-fs와 함께, -ut 불필요)
structty query.pdb -fs result.m8 -m aligned

# 기능 4: 구조 정렬 색상 (nearest-neighbor fallback, threshold=10Å)
structty query.pdb target.pdb -ut matrix.tsv -m aligned

# 기능 5: Conservation 색상
structty protein.pdb --msa alignment.fasta -m conservation

# 기능 3: Foldseek hit 탐색 — 자동 다운로드 (상위 10개)
structty query.pdb -fs result.m8

# 기능 3: Foldseek hit 탐색 — 로컬 DB 사용
structty query.pdb -fs result.m8 --db-path /data/pdb

# 기능 6: 마우스 hover (항상 활성화)
structty protein.pdb

# 복합 사용
structty query.pdb -fs result.m8 -m aligned
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

1. **RenderPoint 전파 필수 확인**: `assign_colors_to_points()`는 RenderPoint를 읽는다. Atom에 값을 설정한 뒤 `project()`에서 RenderPoint로 복사하지 않으면 색상이 전혀 반영되지 않는다. 특히 `is_interface`, `is_aligned`는 기본값이 `false`이므로 복사가 누락되면 전체가 dim 색상으로만 표시된다.
2. **Coil 굵기**: 단일 픽셀은 너무 가늘고 5픽셀 크로스는 너무 두껍다. 3픽셀(중심 + 2픽셀) 패턴을 사용한다.
3. **Coil-Helix/Sheet 접합부**: 경계 잔기를 ss_atoms와 screen_atoms 양쪽에 등록하여 끊김을 방지한다. 이중 등록이 의도된 동작이다.
4. **마우스 hover**: 클릭이 아닌 이동 이벤트(`REPORT_MOUSE_POSITION`)로 잔기 정보를 갱신한다. 패널 Residue Info 섹션은 항상 고정 높이를 유지하여 화면 레이아웃이 변하지 않도록 한다.
5. **패널 고정 크기**: `Panel::set_hover_residue()` 호출 시 빈 칸으로 패딩하여 항상 동일한 줄 수를 출력한다.
6. **`-fs`와 `-ut` 상호 배제**: `-fs`가 있으면 hit의 U/T를 자동 사용하고 `-ut` 인수는 무시한다. `-fs` 없이 `-m aligned`를 사용하려면 반드시 `-ut`가 필요하다.
7. **nearest-neighbor threshold**: `-fs` 없이 `-ut`만 사용하는 경우 10.0Å을 기본값으로 한다.
8. **자동 다운로드**: `curl` → `wget` 순으로 시도한다. 둘 다 없으면 패널에 "curl/wget not found" 표시. 네트워크 오류는 파일 크기 0 또는 파일 미생성으로 감지한다.
9. **기존 기능 무결성**: 기존 `-m protein/chain/rainbow` 모드와 `-s` (이차구조 표시)는 변경 없이 동작해야 함.
10. **Foldseek hit 전환 시 메모리**: 매 hit 전환마다 target Protein 객체를 새로 생성하므로 이전 객체 명시적 해제 필요.
