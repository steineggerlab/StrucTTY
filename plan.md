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
    char residue_name[4] = {};  // 최대 3글자(GLU, ALA 등) + null terminator
                                // std::string 금지: RenderPoint는 매 프레임 대량
                                // 생성·소멸하므로 heap allocation이 발생하면 안 됨
};
```

**`RenderPoint`에 `std::string`을 사용해서는 안 된다.** `RenderPoint`는 렌더링 루프의 핫 패스에서 매 프레임 화면 전체 픽셀 수만큼 생성·소멸한다. `std::string` 멤버가 있으면 픽셀마다 heap allocation이 발생하여 렌더링 성능이 크게 저하된다. 잔기 이름은 최대 3글자이므로 `char[4]`로 충분하다. 복사 시에는 `strncpy(rp.residue_name, atom.residue_name.c_str(), 3)`을 사용한다.

### 0-3. `src/visualization/Screen.cpp` — `project()` 함수 내 RenderPoint 전파 목록

`project()` 함수에서 Atom → RenderPoint로 값을 복사할 때, **아래 필드 전체**를 빠짐없이 복사해야 한다:

```cpp
rp.bfactor            = atom.bfactor;
rp.is_interface       = atom.is_interface;
rp.is_aligned         = atom.is_aligned;
rp.conservation_score = atom.conservation_score;
rp.residue_number     = atom.residue_number;
strncpy(rp.residue_name, atom.residue_name.c_str(), 3);
rp.residue_name[3]    = '\0';
```

이 복사가 누락되면 해당 기능의 색상 분기가 RenderPoint 기본값(false, 0.0f, -1)만 읽게 되어 전체가 단일 색상으로 표시된다. **구현 후 반드시 각 필드가 복사되는지 코드에서 확인한다.**

**helix/sheet geometry atom의 `residue_number` 처리:**
`calculate_ss_points()`에서 생성되는 geometry atom은 실제 CA atom이 아니므로 기본값 `residue_number = -1`을 가진다. 이 atom들도 화면에 렌더링되는 픽셀을 생성하므로, hover 기능이 이 픽셀들에서도 동작하려면 geometry atom 생성 시 해당 segment의 대표 잔기 정보를 설정해야 한다.

```
helix spiral geometry atom:  segment 내 가장 가까운 CA atom의 residue_number/residue_name 사용
sheet ribbon geometry atom:  pair의 pa atom의 residue_number/residue_name 사용
```

즉 `bfactor`/`conservation_score` 보간과 동일한 시점에, `residue_number`와 `residue_name`도 함께 설정한다. 보간하지 않고 구간 내 가장 가까운 원자의 값을 그대로 사용한다.

### 0-4. `src/visualization/Palette.hpp` 색상 확장

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

### 0-5. `src/structure/Parameters.hpp` CLI 인수 확장

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

## 기능 7: 헬릭스/시트 영역의 Coil 제거 + Coil 굵기 조정

### 문제 원인
`StructureMaker::calculate_ss_points()`에서 helix/sheet 지오메트리를 생성하지만, 해당 잔기들의 원본 CA 원자들도 coil 리스트에 그대로 남아 있어 이중으로 그려진다.

### 수정 대상: `src/structure/StructureMaker.cpp`

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

### 수정 대상: `src/visualization/Screen.cpp` — coil 굵기

`project()` 또는 `draw_line()` 함수에서 coil 포인트 렌더링 방식 변경:

```
현재(이전): coil 포인트마다 크로스(+) 패턴으로 5픽셀 이상 그림  ← 너무 두꺼움
이번 목표:  coil 포인트를 3픽셀 크로스(중심 + 상하 또는 좌우 1픽셀씩)로 그림
            단일 픽셀은 너무 가늘어 구조가 보이지 않으므로 사용하지 않는다.
```

구체적으로: `if (point.structure == 'x')` 분기에서 중심 픽셀 + 인접 2픽셀(총 3픽셀) 패턴을 사용한다. 5픽셀 크로스와 단일 픽셀의 중간값이다.

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

**주의: `Atom.is_interface`를 설정한 뒤, 해당 값이 `project()` 과정에서 RenderPoint로 복사되는지 반드시 확인한다 (0-3절 참조).**

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
  hit의 U/T 변환 후 3D 거리 < 5.0Å이면 양쪽 모두 Atom.is_aligned = true
```

**`-fs`가 없는 경우 (nearest-neighbor fallback):**

```
UT 변환 후 nearest-neighbor 방식으로 대응:
  protein1의 각 CA 원자에 대해
    protein2의 모든 CA 원자 중 최소 거리 탐색
    최소 거리 < 10.0Å이면 양쪽 모두 is_aligned = true
```

threshold를 10.0Å으로 설정한다. alignment string 없이 nearest-neighbor만으로 대응하는 경우 5.0Å은 너무 엄격하여 거의 모든 잔기가 non-aligned로 표시된다. 10.0Å은 구조적으로 유사한 영역을 포착하는 실용적인 기준이다.

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

### 색상 적용: `src/visualization/Screen.cpp`

```cpp
if (mode == "aligned") {
    color_id = point.is_aligned ? 45 : 46;  // 정렬됨: 밝은 초록 / 비정렬: dim
}
```

### CLI

```bash
# -fs와 함께 사용 (alignment string 기반, -ut 불필요)
structty query.pdb -fs result.m8 -m aligned

# -ut 단독 사용 (nearest-neighbor fallback, threshold=10Å)
structty query.pdb target.pdb -ut matrix.tsv -m aligned
```

---

## 기능 5: MSA Conservation Score 색상 표시

### 입력 파일 형식 명세

**`--msa`에 전달하는 파일은 반드시 a3m 또는 FASTA 형식의 다중 서열 정렬(MSA) 파일이어야 한다.**

Foldseek/MMseqs2의 `.m8` 파일은 pairwise 검색 결과 테이블이므로 MSA가 아니다. conservation 계산에 사용할 수 없다.

**MMseqs2로 MSA 파일을 생성하는 방법 (사용자가 미리 실행해야 함):**

```bash
# easy-search 실행 시 tmp 파일 보존
mmseqs easy-search query.fasta targetDB result.m8 tmp \
  --remove-tmp-files 0

# tmp 안의 내부 DB로 a3m 생성
mmseqs result2msa tmp/query tmp/targetDB tmp/result tmp/msaDB \
  --msa-format-mode 5

mmseqs convert2fasta tmp/msaDB alignment.a3m
```

또는 DB를 직접 생성하는 방법:

```bash
mmseqs createdb query.fasta queryDB
mmseqs search queryDB targetDB resultDB tmp -a
mmseqs result2msa queryDB targetDB resultDB msaDB --msa-format-mode 5
mmseqs convert2fasta msaDB alignment.a3m
```

**`--msa-format-mode` 옵션:**

| 값 | 출력 형식 |
|---|---|
| 3 | FASTA |
| 5 | A3M (권장) |
| 6 | A3M (alignment 정보 포함) |

a3m 형식(mode 5 또는 6)을 권장한다. 소문자로 표기된 삽입(insertion) 잔기를 포함하므로 쿼리 기준 매핑이 명확하다.

### 신규 파일: `src/structure/MSAParser.hpp` / `MSAParser.cpp`

**MSAParser 클래스 설계:**

```cpp
class MSAParser {
public:
    // a3m 또는 FASTA 형식 파일 로드
    // .m8 파일은 지원하지 않으며, 로드 시도 시 에러 메시지 출력
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
- 첫 문자 `>` → FASTA 또는 A3M으로 처리
- 첫 줄이 탭 구분 숫자열 → `.m8` 파일로 판단, 에러 메시지 출력 후 `false` 반환
  - 에러 메시지: `"Error: --msa expects an a3m/fasta MSA file. .m8 search result files cannot be used for conservation scoring. Use 'mmseqs result2msa' to generate an MSA first."`
- `# STOCKHOLM` 헤더 → Stockholm (향후 지원, 현재는 에러 메시지)

**A3M 처리:**
- 소문자 잔기(삽입)는 다른 서열의 갭으로 처리
- 쿼리(첫 번째 서열) 기준 갭이 없는 위치만 `msa_pos_to_query_idx`에 등록

**Shannon Entropy 계산:**
```
각 MSA 위치 i에 대해:
  갭 제외 빈도 계산
  H(i) = -Σ_aa [f(aa,i) * log2(f(aa,i))]
  conservation(i) = 1 - H(i) / log2(20)
```

### 적용: `src/structure/Protein.cpp`

`apply_conservation_scores(const std::vector<float>& scores)` 함수 추가.
매핑: 0-based 인덱스 순서 기준 (PDB 잔기 번호가 불연속적일 수 있으므로).

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

```bash
# MMseqs2로 생성한 a3m 파일 사용
structty protein.pdb --msa alignment.a3m -m conservation
```

---

## 기능 6: 커서(마우스) 기반 잔기 정보 패널 표시

### ncurses 마우스 활성화: `src/visualization/Screen.cpp`

`init()` 또는 생성자에서:
```cpp
mousemask(ALL_MOUSE_EVENTS | REPORT_MOUSE_POSITION, nullptr);
mouseinterval(0);

// ※ 중요: escape sequence는 mousemask() 반환값과 무관하게 항상 전송한다.
// VSCode 내장 터미널, Windows Terminal(SuperShell) 등 많은 현대 터미널은
// 실제로 xterm-1003 마우스 모드를 지원하지만 terminfo에 올바르게 등록되지 않아
// mousemask()가 0을 반환한다. 반환값에 따라 escape sequence를 조건부로 전송하면
// 이런 터미널에서 마우스가 전혀 동작하지 않는다.
printf("\033[?1003h");  // 모든 마우스 이동 추적 활성화 (xterm 계열)
fflush(stdout);
mouse_supported = true;  // escape sequence를 보냈으므로 지원으로 간주
```

프로그램 종료 시(소멸자 또는 cleanup 함수)에는 반드시 추적을 해제한다:
```cpp
printf("\033[?1003l");  // 마우스 이동 추적 비활성화
fflush(stdout);
```

이 escape sequence 없이는 `KEY_MOUSE` 이벤트가 클릭 시에만 발생하고, 이동 시에는 발생하지 않는다.

### 이벤트 처리 방식

**마우스 이동(hover) 이벤트만 사용한다. 클릭 이벤트는 잔기 정보 표시에 사용하지 않는다.**

ncurses에서 `REPORT_MOUSE_POSITION`을 활성화하면 마우스 이동 시 `KEY_MOUSE` 이벤트가 발생한다. 단, 터미널에 따라 이동 이벤트 지원 여부가 다르다.

**⚠️ `flushinp()` 호출 순서 주의:**

`handle_input(int key)` 내부에서 `flushinp()`를 switch 이전에 일괄 호출하는 패턴은 `KEY_MOUSE` 처리 시 문제를 일으킨다. 일부 ncurses 구현에서 `flushinp()`가 내부 mouse event 큐까지 비워 `getmouse()`가 `ERR`를 반환한다.

따라서 반드시 **`getmouse()`를 먼저 호출하고 그 이후에 `flushinp()`를 호출**하거나, `KEY_MOUSE` 케이스에서는 `flushinp()`를 생략한다:

```cpp
case KEY_MOUSE: {
    MEVENT event;
    if (getmouse(&event) == OK) {  // ← flushinp() 이전에 호출
        update_hover_info(event.x, event.y);
    }
    // KEY_MOUSE에서는 flushinp() 호출하지 않음
    // (또는 getmouse() 이후에 호출)
    break;
}
```

`mousemask()` 반환값이 0이어도 escape sequence를 이미 전송했으므로 KEY_MOUSE는 수신될 수 있다. 패널의 Residue Info 섹션은 항상 유효 정보 또는 "-"를 표시한다.

### main loop 설계: 마우스 이벤트 시 전체 재렌더링 생략

**현재 main loop 구조의 문제점:**

```cpp
while(run) {
    screen.draw_screen();   // clear() + 전체 재렌더링 (50~100ms)
    run = screen.handle_input();
}
```

`REPORT_MOUSE_POSITION` 모드에서는 마우스 이동 1픽셀마다 `KEY_MOUSE` 이벤트가 발생한다. 위 구조에서는 마우스 이동 시 매번 `draw_screen()` → `clear()` → 전체 구조 재렌더링이 실행되어 극심한 화면 깜빡임이 발생하고, 패널의 잔기 정보가 지속적으로 덮어써진다.

**개선 설계: `handle_input()`에 재렌더링 필요 여부 반환 추가**

`handle_input()`이 "전체 재렌더링이 필요한가"를 caller에게 알릴 수 있어야 한다. 마우스 이벤트는 패널 부분 갱신만 필요하고, 키보드 이벤트(회전·이동·줌)는 전체 재렌더링이 필요하다.

구현 옵션:

```cpp
// 옵션 A: out 파라미터
bool Screen::handle_input(bool& needs_redraw);

// 옵션 B: 반환 구조체
struct InputResult { bool keep_running; bool needs_redraw; };
InputResult Screen::handle_input();
```

main loop:
```cpp
bool needs_redraw = true;
while(run) {
    if (needs_redraw) {
        screen.draw_screen();
    }
    run = screen.handle_input(needs_redraw);
    // KEY_MOUSE → needs_redraw = false (패널 부분 갱신만 수행)
    // 키보드   → needs_redraw = true  (전체 재렌더링 필요)
}
```

이 구조로 마우스 이동 시 구조가 재렌더링되지 않으므로 깜빡임이 없어지고, `hover_info` 갱신 → 패널 부분 `refresh()`만 수행된다.

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

### hover 정보 검색: `Screen::update_hover_info(int mx, int my)`

```
screenPixels에서 좌표 (mx, my)에 해당하는 RenderPoint 검색
발견 시: Panel에 잔기 정보 전달
미발견 시: Panel 잔기 정보 섹션 클리어
```

현재 `screenPixels`는 z-buffer 정렬된 2D 배열이므로 O(1) 접근 가능.

**유효 픽셀 판단 기준:** `residue_number >= 0`인 픽셀만 유효로 처리한다.

helix/sheet geometry atom은 실제 CA atom이 아니므로 기본적으로 `residue_number = -1`이지만, 0-3절에 명시된 대로 `calculate_ss_points()`에서 geometry atom 생성 시 반드시 인접 CA atom의 `residue_number`와 `residue_name`을 설정한다. 이 설정이 올바르게 이루어지면 coil뿐 아니라 helix, sheet 위에서도 hover가 동작한다.

결과적으로 구조가 렌더링된 모든 픽셀에서 `residue_number >= 0`이 보장되어야 한다.

**⚠️ non-braille z-buffer 전체 복사 필수:**

non-braille 모드의 z-buffer resolve 코드에서 `depth`, `pixel`, `color_id`만 복사하고 나머지 필드를 누락하면 `screenPixels`의 `residue_number`가 항상 -1로 남는다:

```cpp
// ❌ 잘못된 예: residue_number 등이 복사되지 않음
if (pt.depth < screenPixels[idx].depth) {
    screenPixels[idx].depth    = pt.depth;
    screenPixels[idx].pixel    = pt.pixel;
    screenPixels[idx].color_id = pt.color_id;
}

// ✅ 올바른 예: 전체 복사
if (pt.depth < screenPixels[idx].depth) {
    screenPixels[idx] = pt;  // RenderPoint 전체 복사
}
```

braille 경로는 `logicalPixels[idx] = pt` 전체 복사이므로 문제 없으나, non-braille 경로도 반드시 동일하게 처리해야 한다.

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
11. **`\033[?1003h`는 무조건 전송**: `mousemask()` 반환값이 0이어도 escape sequence를 반드시 전송한다. VSCode, Windows Terminal 등 현대 터미널은 terminfo와 무관하게 실제로 xterm-1003 마우스를 지원한다. 반환값으로 조건부 처리하면 주요 환경에서 마우스가 전혀 동작하지 않는다.
12. **`flushinp()` 위치**: `KEY_MOUSE` 케이스에서는 반드시 `getmouse()`를 `flushinp()` 보다 먼저 호출한다. 일부 ncurses 구현에서 `flushinp()`가 내부 mouse event 큐까지 flush하여 `getmouse()`가 실패한다. `flushinp()`를 switch 진입 전에 일괄 호출하는 패턴은 KEY_MOUSE와 함께 사용 시 반드시 피한다.
13. **마우스 이벤트 시 전체 재렌더링 금지**: `REPORT_MOUSE_POSITION`은 이동마다 이벤트를 생성한다. main loop에서 모든 이벤트 후 `draw_screen()`(= `clear()` + 전체 재렌더링)을 실행하면 극심한 깜빡임이 발생한다. `handle_input()`이 "재렌더링 필요 여부"를 반환하도록 설계하고, KEY_MOUSE 이벤트 후에는 패널 부분 갱신(`refresh()`)만 수행한다.
14. **non-braille z-buffer 전체 복사**: non-braille 경로의 z-buffer resolve에서 `depth/pixel/color_id`만 복사하면 `residue_number`가 항상 -1로 남아 hover가 동작하지 않는다. `screenPixels[idx] = pt`로 RenderPoint 전체를 복사해야 한다.
6. **`-fs`와 `-ut` 상호 배제**: `-fs`가 있으면 hit의 U/T를 자동 사용하고 `-ut` 인수는 무시한다. `-fs` 없이 `-m aligned`를 사용하려면 반드시 `-ut`가 필요하다.
7. **nearest-neighbor threshold**: `-fs` 없이 `-ut`만 사용하는 경우 10.0Å을 기본값으로 한다. 5.0Å은 alignment 정보 없는 상황에서 지나치게 엄격하다.
8. **자동 다운로드**: `curl` → `wget` 순으로 시도한다. 둘 다 없으면 패널에 "curl/wget not found" 표시. 네트워크 오류는 파일 크기 0 또는 파일 미생성으로 감지한다.
9. **기존 기능 무결성**: 기존 `-m protein/chain/rainbow` 모드와 `-s` (이차구조 표시)는 변경 없이 동작해야 함.
10. **Foldseek hit 전환 시 메모리**: 매 hit 전환마다 target Protein 객체를 새로 생성하므로 이전 객체 명시적 해제 필요.
