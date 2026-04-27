# StrucTTY 프로젝트 연구 문서

## 1. 프로젝트 목적 및 전체 구조 요약

**StrucTTY**는 터미널에서 단백질 3D 구조를 인터랙티브하게 시각화하는 C++17 CLI 도구다.  
유니코드 점자(Braille) 문자를 이용한 서브픽셀 렌더링(일반 문자 대비 8× 해상도)과 ncurses 기반 TUI를 통해 최대 9개 단백질 구조를 동시에 표시하며, Foldseek·FoldMason과 연동하여 구조 비교 및 MSA 분석을 지원한다.

- **언어:** C++17
- **라이선스:** GPLv3
- **지원 플랫폼:** Linux, macOS (Windows 빌드 현재 비활성화됨)
- **전체 구현 코드:** 약 5,778 lines (`.cpp` 파일 합계)

### 디렉터리 구조

```
StrucTTY/
├── CMakeLists.txt           # 루트 빌드 설정
├── lib/
│   ├── CMakeLists.txt       # 써드파티 빌드 설정
│   ├── gemmi/               # git submodule: PDB/mmCIF 파서
│   └── lodepng/             # git submodule: PNG 인코더
├── src/
│   ├── structty.cpp         # main() 진입점 (261 lines)
│   ├── structure/           # 구조 데이터 및 파싱
│   │   ├── Atom.hpp
│   │   ├── Protein.cpp/hpp  # 핵심 단백질 데이터 구조 (851 lines)
│   │   ├── Parameters.cpp/hpp
│   │   ├── FoldseekParser.cpp/hpp
│   │   ├── FoldseekDBReader.cpp/hpp
│   │   ├── FoldMasonParser.cpp/hpp
│   │   ├── MSAParser.cpp/hpp
│   │   ├── PDBDownloader.cpp/hpp
│   │   ├── SSPredictor.cpp/hpp
│   │   └── StructureMaker.cpp/hpp
│   ├── visualization/       # 렌더링 및 UI
│   │   ├── Screen.cpp/hpp   # 메인 렌더 엔진 (1942 lines)
│   │   ├── Camera.cpp/hpp
│   │   ├── Panel.cpp/hpp    # 정보 패널 (696 lines)
│   │   ├── Palette.hpp      # xterm-256 컬러 팔레트 (463 lines)
│   │   └── RenderPoint.hpp
│   └── utils/
│       ├── Benchmark.hpp
│       ├── Common.hpp
│       └── Curses.hpp
├── example/                 # 예제 구조 파일 (.cif), Foldseek/FoldMason 결과
├── plan_replace_ncurses.md  # ncurses 제거 계획 문서
└── THIRD_PARTY_NOTICES.md
```

---

## 2. 핵심 모듈/컴포넌트와 각각의 역할

### 2.1 structure/

#### `Atom.hpp`
단일 원자(Cα)를 표현하는 간단한 struct.

| 필드 | 타입 | 설명 |
|------|------|------|
| `x, y, z` | `float` | 3D 좌표 |
| `structure` | `char` | 2차 구조 ('x'=coil, 'H'=helix, 'S'=sheet) |
| `bfactor` | `float` | AlphaFold pLDDT 신뢰도 점수 |
| `is_interface` | `bool` | 체인 간 인터페이스 영역 여부 |
| `is_aligned` | `bool` | 정렬 영역 시각화용 플래그 |
| `conservation_score` | `float` | Shannon 엔트로피 기반 보존도 (-1=미설정) |
| `residue_number` | `int` | 호버 정보용 잔기 번호 |
| `residue_name[4]` | `char[]` | 잔기 3-letter code |

#### `Protein.cpp` (851 lines)
핵심 단백질 데이터 구조. Gemmi를 통해 PDB/mmCIF 파일을 읽고, 두 가지 원자 집합을 관리한다.
- `init_atoms`: 원본 상태 Cα 좌표 (체인별 `std::map<char, std::vector<Atom>>`)
- `screen_atoms`: 변환 적용 후 디스플레이용 좌표 (동일 구조)
- 주요 메서드:
  - `load_data()`: Gemmi로 PDB/mmCIF 로드, HELIX/SHEET 레코드 파싱
  - `compute_interface()`: 체인 간 Cα-Cα 거리 8.0 Å 이하 접촉 계산
  - `compute_aligned_regions_from_aln()`: Foldseek 정렬 문자열로 정렬 잔기 표시
  - `compute_aligned_regions_nn()`: 최근접이웃 fallback 정렬
  - `apply_conservation_scores()`: MSA Shannon 엔트로피 점수 적용
  - `apply_ut_to_init_atoms()`: UT 회전-이동 행렬 적용
  - `compute_bounding_box()`: 3D 장면 경계 계산

#### `Parameters.cpp` (137 lines)
CLI 인수 파싱. 최대 9개 입력 파일 지원.

| 옵션 | 설명 |
|------|------|
| `-m/--mode` | 색상 모드 (protein/chain/rainbow/plddt/interface/conservation/aligned) |
| `-c/--chains` | 체인 선택 TSV 파일 |
| `-s/--structure` | 2차 구조 기하학 표시 여부 |
| `-ut/--utmatrix` | U/T 행렬 TSV 파일 경로 |
| `--msa` | MSA 파일 (FASTA/A3M) |
| `-fs/--foldseek` | Foldseek easy-search 결과 .m8 파일 |
| `--db/--db-path` | Foldseek DB 직접 읽기 |
| `-fm/--foldmason` | FoldMason JSON/FASTA 결과 파일 |
| `-n/--nopanel` | 정보 패널 비표시 |
| `-b/--benchmark` | 벤치마크 모드 |

#### `FoldseekParser.cpp` (181 lines)
Foldseek `easy-search` 출력 `.m8` 포맷 파싱. 4종 컬럼 포맷 지원:
- **12 cols:** 기본 (정렬 없음)
- **17 cols:** qaln/taln 포함
- **21 cols:** `alis` 포맷 (qaln, taln, alns, tseq, taxid, taxname) — Foldseek 기본
- **29 cols:** 전체 포맷 (U/T 행렬, lddt, qtmscore, ttmscore 포함)

`FoldseekHit` struct: query, target, fident, alnlen, qstart/qend/tstart/tend, evalue, prob, U[9], T[3], 정렬 문자열

#### `FoldseekDBReader.cpp` (279 lines)
Foldseek 구조 DB(`_ca`, `_seq` 파일) 직접 읽기. 히트 기반 선택적 스캔으로 메모리 효율화 (AFDB50 1000 히트에 약 152KB).

#### `FoldMasonParser.cpp` (295 lines)
FoldMason JSON/FASTA MSA 결과 파싱.
- JSON 구조: entries(name, aa, ss, ca coords), scores, tree, statistics
- `build_query_col_map()`: 잔기 인덱스 → MSA 열 인덱스 매핑
- `build_aligned_pairs()`: gap-free 잔기 쌍 추출 (Kabsch 중첩용)
- `compute_column_entropy()`: 열별 보존도 계산 (aa 또는 ss 기준)

#### `MSAParser.cpp` (127 lines)
FASTA/A3M 포맷 MSA 파싱 및 Shannon 엔트로피 기반 보존도 계산.
- A3M 정규화: 소문자(삽입) 제거, 대문자 + 갭 유지
- Shannon 엔트로피: `H(i) = -Σ f·log₂(f)`
- 보존도: `conservation = 1 - H(i) / log₂(20)` (0.0=가변, 1.0=보존)

#### `PDBDownloader.cpp` (419 lines)
멀티-DB 타겟 탐지 및 다운로드 관리. 지원 DB: PDB, AlphaFold, ESMAtlas30, CATH50, BFVD(로컬/공식), GMGCL, TED, BFMD. 캐시 경로: `~/.cache/structty/pdb/`

#### `SSPredictor.cpp` (156 lines)
거리 및 토션 각도 기반 2차 구조 예측 (3-잔기 윈도우 투표).
- 나선 기준: d₁₃∈[4.8,6.0]Å, d₁₄∈[4.8,6.0]Å, torsion∈[35,75]°
- 시트 기준: d₁₃≥6.4Å, d₁₄≥8.5Å, torsion∈[110,180]∪[-180,-110]°
- 노이즈 제거 스무딩 (최소 길이: 나선=4, 시트=3)

#### `StructureMaker.cpp` (323 lines)
2차 구조 시각화 기하학 생성.
- 나선 렌더링: PCA로 축 계산(거듭제곱 반복법), Cα로부터 실린더 생성 (반지름=2.5Å)
- 베타 시트 렌더링: Cα 백본에 수직 리본 (폭=4, 스텝=0.28)

---

### 2.2 visualization/

#### `Screen.cpp` (1942 lines)
메인 렌더 엔진. 모든 단백질, 카메라, 패널, 컬러 쌍을 관리.

**핵심 점자 렌더링 원리:**  
터미널 문자 1개 = 논리 픽셀 2×4 = 8배 해상도. `logicalPixels` 버퍼에 RenderPoint 저장.

**주요 메서드:**
- `set_protein()`: Protein 인스턴스 로드
- `load_next_hit()`: Foldseek 히트 탐색 + Kabsch 중첩
- `apply_foldmason_superposition()`: FoldMason 데이터로 Kabsch 정렬
- `draw_screen()`: 메인 렌더 함수
- `handle_input()`: 키보드/마우스 입력 처리
- `update_hover_info()`: 마우스 잔기 탐지 → 패널 부분 갱신

#### `Camera.cpp` (111 lines)
3D → 2D 원근 투영 관리.
- FOV = 90° (하드코딩)
- 초점 오프셋: 일반 뷰=3.0, 스크린샷=10.0
- lodepng를 통한 PNG 스크린샷 저장, RGBA 이미지 렌더링

#### `Panel.cpp` (696 lines)
정보 패널 표시 (터미널 우측). 구성 섹션:
- **파일 정보**: 단백질별 체인 및 원자 수
- **Foldseek 히트**: 현재 히트 #/전체, 타겟명, evalue, prob, lddt, qtmscore, 정렬 방법, 다운로드 상태
- **호버 잔기**: 체인, 잔기명/번호, 구조 타입, bfactor, 보존도 점수 (4줄 고정)
- **FoldMason**: 엔트리 수, 정렬 방법

#### `Palette.hpp` (463 lines)
xterm-256 컬러 팔레트 전체 정의. 250개 이상 ncurses 컬러 쌍 사전 할당.

| 팔레트 | 쌍 번호 | 용도 |
|--------|---------|------|
| 단백질 색 (9색) | 1-9, 11-19, 101-109, 120-128, 200-208 | near/mid/far 깊이 변형 포함 |
| 체인 색 (15색) | 21-35 + 깊이 변형 | |
| 무지개 그라디언트 (20색조) | 51-70 + 깊이 변형 | |
| pLDDT 신뢰도 (4밴드) | 71-74 + 깊이 변형 | |
| 보존도 (10단계 그라디언트) | 75-84 + 깊이 변형 | |
| 인터페이스 | 43-44, 237-240, 165-166 | |
| 정렬 영역 | 45-46, 101-109, 110, 241-250 | |
| xterm-256 RGBA LUT | 256 엔트리 | PNG 내보내기용 |

#### `RenderPoint.hpp` (34 lines)
논리 픽셀 데이터 struct.

| 필드 | 설명 |
|------|------|
| `x, y` | 픽셀 위치 |
| `z` | Z-버퍼용 깊이 |
| `pixel` | 점자 문자 |
| `color_id` | ncurses 컬러 쌍 번호 |
| `chainID` | 체인 식별자 |
| `structure` | 2차 구조 타입 |
| `is_interface`, `is_aligned` | 기능 플래그 |
| `is_conserved` | 보존도 점수 |
| `residue_number`, `residue_name[4]` | 호버 정보 |
| `depth_band` | 깊이 밴드 (0=near, 1=mid, 2=far) |

---

### 2.3 utils/

| 파일 | 역할 |
|------|------|
| `Common.hpp` | `current_timestamp()` — 밀리초 정밀도 타임스탬프 |
| `Curses.hpp` | OS별 ncurses 헤더 감지 (다중 include 경로 대응) |
| `Benchmark.hpp` | 프레임 타이밍 CSV 출력, FPS/레이턴시/렌더시간 추적 |

---

## 3. 빌드 시스템 및 의존성

### 3.1 CMake 구성

**루트 `CMakeLists.txt`:**
- 최소 버전: CMake 3.15
- C++17 표준 (required)
- Wide character ncurses 필수 (`CURSES_NEED_WIDE TRUE`)
- `file(GLOB_RECURSE ...)` 로 `src/structure/*.cpp`, `src/visualization/*.cpp` 수집
- 실행 파일 타겟: `StrucTTY`
- 링크: `gemmi::gemmi_cpp`, `lodepng`, `${CURSES_LIBRARIES}`

**`lib/CMakeLists.txt`:**
- gemmi: `add_subdirectory(gemmi)` (CMake 타겟 `gemmi::gemmi_cpp` 제공)
- lodepng: `lodepng.cpp` 단일 파일을 정적 라이브러리로 빌드 (`POSITION_INDEPENDENT_CODE ON`)

**`src/CMakeLists.txt`:**
- `structty_core` 정적 라이브러리 생성 (최종 바이너리에서 미사용 — 루트 CMakeLists가 직접 소스 glob)

### 3.2 외부 의존성

| 라이브러리 | 버전/형태 | 라이선스 | 용도 |
|-----------|---------|---------|------|
| **Gemmi** | git submodule (github.com/project-gemmi/gemmi) | MPL-2.0 | PDB/mmCIF 파싱, 구조 읽기 |
| **LodePNG** | git submodule (github.com/lvandeve/lodepng) | zlib | PNG 스크린샷 인코딩 |
| **ncurses (wide)** | 시스템 패키지 | Public Domain | 터미널 렌더링, 컬러 쌍, 마우스 입력 |
| **C++ stdlib** | C++17 | — | chrono, fstream, iostream, algorithm 등 |

**플랫폼별 설치:**
- Linux: `libncursesw-dev`
- macOS: `brew install ncurses` + CMake homebrew 경로 설정
- Windows: MSYS2 `mingw-w64-x86_64-ncurses` (현재 빌드 비활성화)

---

## 4. 주요 데이터 흐름

### 4.1 초기화 및 로딩 단계 (`structty.cpp`)

```
main()
  ├─ Parameters::parse() — CLI 인수 파싱
  ├─ initscr() + ncurses 설정
  ├─ Screen 초기화 (width, height, show_structure, mode)
  ├─ 각 입력 파일에 대해:
  │   └─ Screen::set_protein() → Protein() → Gemmi로 PDB/mmCIF 로드
  ├─ Screen::normalize_proteins() — 스케일 정규화 (최대 차원 3.0 Å 최소)
  ├─ Screen::update_total_len_ca()
  │
  ├─ [-ut] → Screen::set_utmatrix() → Protein::apply_ut_to_init_atoms()
  ├─ [-m interface] → Screen::compute_interface_all()
  │                   → Protein::compute_interface_pair() (8.0 Å 임계값)
  ├─ [-m aligned, -ut 만 있을 때] → Screen::compute_aligned_all()
  │                                 → Protein::compute_aligned_regions_nn()
  ├─ [--db] → Screen::open_foldseek_db() → FoldseekDBReader::prepare()
  │
  ├─ [-fs] → FoldseekParser::load(.m8)
  │           ├─ 다중 입력: CLI 타겟 파일명 기준 히트 필터링
  │           └─ 단일 입력: Screen::load_next_hit(+1)
  │                         → UT 변환 (Kabsch, alis 포맷일 때)
  │                         → 정렬 영역 계산 (aln-string 또는 nn)
  │
  ├─ [-m conservation, --msa] → MSAParser::load() + compute_conservation()
  │                             → Screen::apply_msa_conservation()
  │
  └─ [-fm] → FoldMasonParser::load_json() or load_fasta()
              ├─ [-m conservation] → compute_column_entropy() → apply_msa_conservation()
              └─ [2+ 구조] → Screen::apply_foldmason_superposition() (Kabsch)
```

### 4.2 렌더링 루프 (`Screen::draw_screen()`)

```
매 프레임:
  ├─ 모든 원자 투영: 3D 월드 → 2D 논리 픽셀 좌표
  │   └─ 원근 투영: z_screen = (z_world - camera_z) / focal_offset
  │   └─ 논리 픽셀 격자에 클램핑 [0, 2×width) × [0, 4×height)
  ├─ 깊이 정렬 + Z-버퍼 (논리 픽셀당)
  ├─ 색상 모드별 컬러 할당:
  │   ├─ protein: 9색 순환 (coil은 dim)
  │   ├─ chain: 체인당 15색
  │   ├─ rainbow: N→C 20단계 색조 그라디언트
  │   ├─ plddt: 4밴드 [≥90, 70-90, 50-70, <50]
  │   ├─ interface: 마젠타(인터페이스), dim(비인터페이스)
  │   ├─ conservation: 10단계 그라디언트 (가변→보존)
  │   └─ aligned: bright(정렬), dim(비정렬)
  ├─ 3-밴드 깊이 안개 적용:
  │   ├─ near (z < base_min): 밝은 변형
  │   ├─ mid: 기본 컬러
  │   └─ far (z > base_max): 색조 유지 어두운 변형
  ├─ 점자 문자 렌더링: 논리 픽셀 2×4 → 터미널 문자 1개
  └─ draw_panel() (패널 활성 시)
```

### 4.3 입력 처리

```
handle_input():
  ├─ 0: 전체 단백질 제어
  ├─ 1-9: 특정 단백질 제어
  ├─ W/A/S/D: 팬 (x/y 이동)
  ├─ X/Y/Z: 해당 축 회전
  ├─ R/F: 줌 인/아웃
  ├─ Q: 종료
  └─ KEY_MOUSE: update_hover_info() → 패널 부분 갱신 (전체 재렌더링 없음)
```

### 4.4 입력/출력 포맷 요약

| 방향 | 포맷 | 설명 |
|------|------|------|
| 입력 | PDB / mmCIF | Gemmi로 파싱, Cα 원자만 추출 |
| 입력 | Foldseek `.m8` (21 col alis) | qaln/taln/alns/tseq + taxid/taxname |
| 입력 | FoldMason JSON | entries(name, aa, ss, ca coords), scores, tree |
| 입력 | MSA FASTA / A3M | A3M 소문자=삽입(제거), Shannon 엔트로피 계산 |
| 입력 | UT matrix TSV | U[9] + T[3] 회전-이동 행렬 |
| 출력 | 터미널 | ncurses + Unicode 점자 + ANSI 컬러 |
| 출력 | PNG | lodepng RGBA 렌더링 |
| 출력 | CSV | 벤치마크 타이밍 (t_ms, tag, key, dt_ms, num_ca, num_file) |

---

## 5. 현재 구현된 기능 목록

### 렌더링 및 시각화
- [x] Unicode 점자 서브픽셀 렌더링 (8× 해상도)
- [x] 최대 9개 단백질 동시 표시
- [x] 7가지 색상 모드 (protein, chain, rainbow, plddt, interface, conservation, aligned)
- [x] 3-밴드 깊이 안개 (near/mid/far, 색조 보존)
- [x] 2차 구조 시각화 (나선 실린더, 시트 리본)
- [x] 팬/회전/줌 키보드 제어
- [x] 마우스 호버 잔기 정보 표시
- [x] PNG 스크린샷 저장

### 입력 및 구조 로딩
- [x] Gemmi 기반 PDB/mmCIF 파싱
- [x] TSV 파일을 통한 체인 선택
- [x] Cα 원자만 추출
- [x] 바운딩 박스 계산 및 장면 정규화

### 정렬 및 비교
- [x] UT 행렬(회전-이동) 적용
- [x] Foldseek 히트 키보드 탐색
- [x] 정렬 문자열(qaln/taln) 시각화
- [x] 최근접이웃 정렬 fallback
- [x] FoldMason MSA Kabsch 중첩
- [x] 구조 정렬 영역 하이라이팅

### Foldseek 연동
- [x] M8 포맷 파싱 (12/17/21/29 컬럼 변형)
- [x] 멀티 DB 지원 (PDB, AlphaFold, ESMAtlas, CATH, BFVD, TED, BFMD, GMGCL)
- [x] 구조 자동 다운로드 (curl/wget)
- [x] 로컬 캐시 (`~/.cache/structty/pdb/`)
- [x] Foldseek DB 직접 읽기 (`_ca`, `_seq` 파일)
- [x] 히트 기반 선택적 스캔 (메모리 효율)

### MSA 분석
- [x] FASTA/A3M 파싱
- [x] Shannon 엔트로피 보존도 점수 계산
- [x] 잔기별 점수 시각화
- [x] FoldMason JSON 지원
- [x] FoldMason FASTA 중첩

### 2차 구조
- [x] PDB/mmCIF HELIX/SHEET 레코드 탐지
- [x] 거리-토션 각도 기반 예측 (대안)
- [x] 나선 실린더 렌더링
- [x] 베타 시트 리본 렌더링

### 사용자 인터페이스
- [x] ncurses 터미널 렌더링
- [x] 키보드 제어 (팬, 회전, 줌, 단백질 선택)
- [x] 마우스 추적 (호버 잔기 정보)
- [x] 정보 패널 (파일, 체인, 히트 정보, 호버, FoldMason 데이터)
- [x] Wide 문자 지원 (Unicode 점자)

### 벤치마킹
- [x] 프레임 타이밍 CSV 출력
- [x] FPS 측정
- [x] 입력→프레임 레이턴시 추적
- [x] 워밍업(200 프레임) + 측정(2000 프레임) 분리

---

## 6. 미구현 또는 TODO로 남겨진 항목

### `plan_replace_ncurses.md`에 명시된 계획 (미실행)

- [ ] **`libstructty-render` 라이브러리 추출** — ncurses 의존성 없는 독립 렌더링 라이브러리
  - **배경:** Martin이 Foldseek 연동 시 ncurses 링크(너무 비대) 및 subprocess 호출(비수용) 방식 모두 거부. 렌더링 엔진을 standalone 라이브러리로 분리 지시.
  - **대상 범위:** 3D→2D 투영, 레이캐스팅/Z-버퍼, RenderPoint 버퍼 생성, ANSI 이스케이프 코드 stdout 출력
  - **비대상:** 패널 레이아웃, 입력 처리, 멀티 구조 TUI 라이프사이클
  - **오픈 질문:** 공개 API 표면(`render_structure(atoms, width, height) -> ansi_string`?), 단일/다중 구조 지원 여부, Foldseek 호출 시점
  - **다음 단계 (미실행):**
    1. 렌더링 파이프라인 내 ncurses 호출 감사
    2. 라이브러리 공개 API 정의
    3. ncurses 링크 없는 별도 CMake 타겟으로 렌더링 로직 추출
    4. ANSI 출력 경로 구현 (비-ncurses 디스플레이 백엔드)
    5. Foldseek 링크 및 호출 검증

### 코드에서 발견된 미구현/제한 사항

- [ ] **Windows 빌드 비활성화** — `.github/workflows/` 내 Windows 빌드 주석 처리됨 (커밋 `40a105f` "Temporarily disable Windows build")
- [ ] **`src/CMakeLists.txt`의 `structty_core` 정적 라이브러리** — 정의되어 있으나 루트 CMakeLists의 최종 바이너리에서 실제로 사용되지 않음 (잠재적 빌드 불일치)
- [ ] **다운로드 진행바** — 상태 메시지 표시는 되나 진행률 표시줄 미구현
- [ ] **멀티스레드 구조 다운로드** — 단일 스레드 I/O만 지원
- [ ] **히트 필터링/정렬 UI** — 점수 표시는 되나 인터랙티브 필터/정렬 미구현
- [ ] **9개 초과 단백질** — `MAX_STRUCT_NUM` 하드코딩

---

## 7. 코드 컨벤션 및 네이밍 규칙

### 네이밍 규칙

| 카테고리 | 규칙 | 예시 |
|---------|------|------|
| **클래스** | PascalCase | `Screen`, `Protein`, `FoldseekParser`, `SSPredictor` |
| **함수/메서드** | snake_case | `set_protein()`, `compute_interface()`, `draw_screen()` |
| **멤버 변수** | snake_case | `screen_width`, `screen_atoms`, `foldseek_hits` |
| **상수** | UPPER_CASE | `FOV`, `PI`, `MAX_STRUCT_NUM`, `MAX_DEPTH` |
| **열거형 타입** | PascalCase | `DBType` |
| **열거형 값** | PascalCase | `DBType::PDB`, `DBType::AlphaFold` |
| **struct 멤버** | snake_case | `fident`, `alnlen`, `conservation_score` |
| **불리언 멤버** | `is_*` / `has_*` 접두사 | `is_interface`, `is_aligned`, `has_transform`, `has_aln` |

### 파일 구성

- **헤더 가드:** `#pragma once` (이식성 주석 없이 일관 사용)
- **include 순서:** 시스템 헤더(`<gemmi/...>`, `<vector>`) → 로컬 헤더(`"Protein.hpp"`)
- **소스/헤더 쌍:** 각 클래스마다 `.cpp` + `.hpp` 쌍 (단, 간단한 struct/헤더전용은 `.hpp`만)

### 메모리 관리

- **장수명 객체:** 수동 `new`/`delete` (Screen, Protein, Camera, Panel)
- **스택 할당:** 임시 객체
- **스마트 포인터:** `std::unique_ptr`는 `structty.cpp`의 `FoldMasonParser`에만 사용
- **참조 전달:** 복사 회피를 위해 `const std::string&`, `const std::vector<...>&`
- **출력용 비const 참조:** `std::vector<...>&` (in-place 채움)

### 데이터 타입 선택

| 용도 | 타입 |
|------|------|
| 3D 좌표 | `float` |
| 인덱스, 카운트 | `int` |
| 2차 구조 타입 | `char` ('x', 'H', 'S') |
| 이름, 경로, 시퀀스 | `std::string` |
| 동적 배열 | `std::vector<>` |
| 체인→원자 룩업 | `std::map<>` |

### 에러 처리

- **최소 예외 사용** — 대부분 불리언 반환
- **실패 시 stderr** — `std::cerr` 출력 후 `false` 반환
- **Assert 미사용** — 방어적 체크만
- **우아한 저하** — 누락 데이터는 -1 값 또는 빈 벡터로 표현

### 아키텍처 패턴

1. **이중 원자 표현** — `init_atoms`(원본) vs `screen_atoms`(변환 적용): 다중 변환 합성 허용, 재계산 불필요
2. **원자별 기능 플래그** — 각 기능(인터페이스, 정렬, 보존도)을 bool/float로 직접 저장: 다중 모드 동시 시각화 가능
3. **논리 픽셀 버퍼** — `logicalPixels` (2×width × 4×height): 서브문자 해상도 구현
4. **컬러 쌍 사전 할당** — `init_color_pairs()`로 250개 이상 ncurses 쌍 초기화: 즉시 컬러 룩업, 깊이 안개 변형 사전 베이킹
5. **선택적 DB 스캔** — FoldseekDBReader가 히트 타겟만 인덱싱: GB 단위 AFDB50에서도 최소 메모리 사용

### 주석 스타일

- Doxygen/JavaDoc 스타일 미사용
- 복잡한 로직 인라인 주석
- **한국어와 영어 혼용** — `structty.cpp`의 기능 설명 주석은 한국어로 작성됨 (예: `// 기능 1: interface 모드일 때 inter-chain interface 계산`)

---

---

## 8. libstructty-render 추출 작업 진행 현황

### Phase 0 완료 (2026-04-24)

`plan_replace_ncurses.md`의 Phase 0 (설계 확정 및 감사)이 완료됐다. 코드 변경 없음.

#### 주요 발견

**ncurses 교체 대상은 4개뿐이다 (`print_screen_braille()` 내)**

| 라인 | 현재 (ncurses) | 교체 대상 (ANSI) |
|------|----------------|-----------------|
| 826 | `getmaxyx(stdscr, rows, cols)` | `width`/`height` 파라미터로 대체 |
| 867 | `attron(COLOR_PAIR(id))` | `\033[38;5;{n}m` |
| 868 | `mvaddstr(row, tx, utf8)` | 위치 기반 string 빌딩 |
| 869 | `attroff(COLOR_PAIR(id))` | `\033[0m` |

나머지 ncurses 호출(init_pair, getch, refresh 등)은 모두 TUI 레이어에 속하며 StrucTTY에 그대로 유지된다.

**color_id → xterm-256 변환 테이블이 이미 존재한다**

`Camera::renderPoint2image()` (`Camera.cpp`:41-71)에 color_id(1–250) → xterm-256 번호의 완전한 변환 로직이 이미 구현되어 있다. Phase 1의 `palette_to_xterm256_fg()` 함수는 이 로직을 Palette.hpp로 추출만 하면 된다.

**depth_band 처리 단순화**

Palette.hpp의 near/mid/far 색상 변형이 이미 별도 xterm-256 번호로 정의되어 있으므로 (`PROTEIN_NEAR_COLORS`, `PROTEIN_FAR_COLORS` 등), ANSI 출력 시 `\033[2m` dim escape 없이 color_id만 변환해도 깊이 안개가 정확히 재현된다.

### Phase 1 완료 (2026-04-24)

`plan_replace_ncurses.md`의 Phase 1 (컬러 시스템 이식)이 완료됐다.

#### 변경 파일

**`src/visualization/Palette.hpp`** — `Palettes` 네임스페이스 끝에 인라인 함수 3개 추가:
- `palette_to_xterm256_fg(int cid) -> int`: color_id(1–250) → xterm-256 전경색 번호. `Camera::renderPoint2image()` 변환 테이블을 정확히 이식.
- `palette_to_ansi_fg_str(int color_id) -> std::string`: `"\033[38;5;{n}m"` 포맷 문자열 반환.
- `palette_to_ansi_reset() -> const char*`: `"\033[0m"` 반환.

**`CMakeLists.txt`** — 빌드 환경 수정 2건:
- `target_include_directories`에 `${CURSES_INCLUDE_DIRS}` 추가 (누락되어 있었음).
- `find_package(ZLIB REQUIRED)` 및 `ZLIB::ZLIB` 링크 추가 (gemmi가 요구하는 zlib가 링커 경로에 없어 빌드 실패하던 문제 수정).

#### 설계 결정

- `RenderPoint.color_id` 값 자체는 변경하지 않았다. 기존 ncurses pair 번호 범위(1–250)를 그대로 사용하며, ANSI 변환 경로만 추가했다.
- 깊이 안개(near/mid/far)는 이미 color_id에 인코딩되어 있으므로 (`PROTEIN_NEAR_COLORS`/`PROTEIN_FAR_COLORS` 등) 별도 `\033[2m` escape 불필요.
- 26개 테스트 케이스(Python spot-check)로 변환 함수 정확성 사전 검증 완료.

#### 빌드 결과

Zero new warnings. 기존 `Protein.cpp` narrowing 경고들은 Phase 1 이전부터 존재하던 것으로, Phase 1과 무관하다.

*최종 업데이트: 2026-04-24*

---

## Phase 2 구현 요약: 렌더링 코어 추출

**브랜치:** `replace-ncurses`  
**완료일:** 2026-04-25

### 목표

`Screen.cpp`의 렌더링 파이프라인(투영·Z-버퍼·색상 할당)을 ncurses와 무관한 `Renderer` 클래스로 분리하여 `src/render/` 아래에 배치.

### 변경 파일

**`src/render/Renderer.hpp`** (신규):
- `RenderAtom` 구조체: Protein/Atom/gemmi/ncurses와 완전 분리된 순수 데이터. `chain_id`는 `std::string` (plan의 `char` 대신 — `Protein::get_atoms()` 키가 `std::string`이므로). per-protein 화면 오프셋 `pan_x`, `pan_y` 필드 추가 (plan 미명시).
- `Renderer` 클래스: 공개 API는 `render(atoms)`, `get_pixels()`, `get_logical_width/height()`, `set_depth_params()`.

**`src/render/Renderer.cpp`** (신규):
- `project_and_fill()`: 원근 투영 + Catmull-Rom 코일 보간 (평탄한 `RenderAtom` 벡터를 protein → chain 순으로 그루핑하여 처리).
- `zbuffer_resolve()`: Z-버퍼 합성, `logical_pixels_` 소유.
- `assign_colors_impl()`: 7개 색상 모드(protein/chain/rainbow/plddt/interface/aligned/conservation) × 3 깊이 밴드 — `Screen::assign_colors_to_points()`와 동일 논리.
- `draw_line_impl()`: `compute_depth_band()` 내부 호출로 min_z/max_z 파라미터 제거.

**`src/visualization/Screen.hpp`** (수정):
- `#include "Renderer.hpp"` 추가, `Renderer renderer_` 멤버 추가.
- `logicalPixels` 벡터 제거, `project()`/`clear_screen()` 선언 제거.
- `to_render_atoms()` 비공개 메서드 선언 추가.

**`src/visualization/Screen.cpp`** (수정):
- 생성자에 `renderer_(width, height, mode, show_structure)` 이니셜라이저 추가.
- `Screen::project()`(~155줄) 및 `clear_screen()` 구현 제거.
- `to_render_atoms()` 구현 추가: `data[]` → `RenderAtom` 벡터 변환 어댑터.
- `draw_screen()`: `calibrate_depth_baseline_first_view()` → `renderer_.set_depth_params(...)` → `renderer_.render(to_render_atoms())` 순으로 위임.
- `print_screen_braille()`, `update_hover_info()`, `camera->screenshot()`: `logicalPixels` → `renderer_.get_pixels()` 교체.

**`CMakeLists.txt`** (수정):
- `src/render/*.cpp`, `src/render/*.hpp` 글로브 추가.
- `target_include_directories`에 `src/render` 경로 추가.

### 설계 결정 및 plan 이탈 사항

| 항목 | plan 명세 | 실제 구현 | 이유 |
|------|-----------|-----------|------|
| `RenderAtom.chain_id` | `char` | `std::string` | `Protein::get_atoms()` 키가 `std::string` |
| `pan_x`, `pan_y` | 미명시 | `RenderAtom`에 추가 | 투영 시 per-protein 화면 오프셋 필요 |
| `generate_braille()` 추출 | Phase 2-6 | Screen에 유지 | 브레일 생성은 ncurses 출력과 결합되어 있어 Phase 3에서 분리하는 것이 적절 |
| `calibrate_depth_baseline` | `project()` 내부 호출 | `draw_screen()`으로 이동 | Renderer 생성 후 `set_depth_params()` 전에 호출해야 하는 순서 보장 |

### 빌드 결과

Zero new warnings (기존 `Protein.cpp` narrowing 경고는 Phase 2 이전부터 존재). 빌드 성공, `1_1CRN.cif -m chain -s` 실행 결과 동일.

*최종 업데이트: 2026-04-25*

---

## Phase 3 구현 요약: ANSI 출력 백엔드

**브랜치:** `replace-ncurses`  
**완료일:** 2026-04-26

### 목표

`Renderer::get_pixels()`가 반환하는 논리 픽셀 버퍼(`logical_pixels_`)를 ANSI escape code 문자열로 변환하는 `AnsiOutput` 클래스 구현. ncurses 미포함.

### 변경 파일

**`src/render/AnsiOutput.hpp`** (신규):
- `AnsiOutput` 클래스: 두 개의 static 메서드로 구성.
  - `to_ansi_string(pixels, logical_width, logical_height) -> std::string`: 픽셀 버퍼 → ANSI string. 터미널 크기 = (logical_width/2) × (logical_height/4).
  - `print_to_stdout(...)`: `setlocale(LC_ALL, "")` 후 `to_ansi_string()` 결과를 `fwrite`로 stdout 출력.

**`src/render/AnsiOutput.cpp`** (신규):
- `Screen::print_screen_braille()`와 동일한 2×4 점자 bitmask 계산 로직 사용.
- `DOT_BITS[2][4]` 상수: Screen.cpp와 동일한 dot bit 배열.
- 각 터미널 셀에서 최소 depth(가장 앞면)의 `color_id`를 선택.
- `Palettes::palette_to_ansi_fg_str(color_id)` + 브레일 UTF-8 3바이트 인코딩(`0xE2, 0xA0|(mask>>6), 0x80|(mask&0x3F)`) + `Palettes::palette_to_ansi_reset()` 조합.
- 빈 셀(bitmask==0 또는 color_id==0)은 공백(' ')으로 채워 열 정렬 유지.
- 각 행 끝에 `'\n'` 추가.
- `out.reserve()` 로 heap reallocation 최소화.

### 설계 결정

| 항목 | 결정 | 이유 |
|------|------|------|
| 커서 위치 초기화 | 미포함 | Foldseek UX 요구사항 Phase 5에서 결정; 호출부에서 필요 시 추가 |
| depth_band far dim | 미적용 | color_id에 이미 어두운 xterm 번호 인코딩됨 (`PROTEIN_FAR_COLORS` 등) |
| 빈 셀 처리 | 공백(' ') | ANSI 출력에서 열 정렬 유지 필요 (ncurses는 커서 이동으로 스킵 가능했으나 ANSI는 행 연속 출력) |

### 빌드 결과

Zero new warnings. 빌드 성공, `1_1CRN.cif -m chain -s` 실행 결과 동일.

## Phase 4 구현 요약: 공개 API 헤더 및 CMake 타겟

**브랜치:** `replace-ncurses`  
**완료일:** 2026-04-26

### 목표

Phase 1–3에서 구현한 `Renderer`/`AnsiOutput`/`Palette`를 ncurses 미링크 독립 정적 라이브러리 (`structty_render`)로 패키징하고, Foldseek이 `add_subdirectory` 방식으로 링크할 수 있는 공개 C++ API 헤더를 작성.

### 변경 파일

**`include/structty_render.h`** (신규):
- `#include "Renderer.hpp"` + `#include "AnsiOutput.hpp"` 포함.
- `namespace structty {}` 내 두 개의 inline 편의 함수:
  - `render_to_ansi(atoms, width, height, mode, show_structure) -> std::string`
  - `render_to_stdout(atoms, width, height, mode, show_structure)`
- 두 함수 모두 `Renderer` 인스턴스 생성 → `render()` 호출 → `AnsiOutput` 변환의 3단계.
- 기본값: `mode="protein"`, `show_structure=false`.

**`src/render/CMakeLists.txt`** (신규):
- `add_library(structty_render STATIC Renderer.cpp AnsiOutput.cpp)`
- PUBLIC include dirs: `${PROJECT_SOURCE_DIR}/include`, `${CMAKE_CURRENT_SOURCE_DIR}` (src/render/), `${PROJECT_SOURCE_DIR}/src/visualization`.
- PRIVATE 의존: `gemmi::gemmi_cpp`. ncurses 링크 없음.

**`CMakeLists.txt`** (수정):
- `add_subdirectory(src/render)` 추가 (lib 다음).
- `APP_SOURCES` glob에서 `src/render/*.cpp` 제거 — `structty_render` 라이브러리가 컴파일하므로 중복 방지.
- `target_link_libraries(StrucTTY ...)` 에 `structty_render` 추가. ncurses는 StrucTTY에만 링크.

### 설계 결정

| 항목 | 결정 | 이유 |
|------|------|------|
| inline vs .cpp | inline 함수 | `structty_render.a`에 컴파일 단위 추가 없이 헤더만으로 완결; 얇은 래퍼이므로 코드 비용 무시 가능 |
| src/render 중복 제거 | glob에서 제외 | `structty_render`가 동일 .cpp를 컴파일하므로 StrucTTY에 직접 포함 시 duplicate symbol 링크 오류 발생 |
| PUBLIC include 전파 | structty_render PUBLIC | StrucTTY가 `src/render/*.hpp` 를 include할 때 별도 PRIVATE include 없이 transitive 경로로 해결 가능 |

### 빌드 및 검증 결과

```
[85%] Linking CXX static library libstructty_render.a
[85%] Built target structty_render
[100%] Built target StrucTTY
```

Zero new warnings. ncurses 심볼 검증:

```
nm -u build/src/render/libstructty_render.a | grep -i ncurses
# → 출력 없음 (exit 1) — ncurses 미의존 확인
```

*최종 업데이트: 2026-04-26*

---

## Phase 5 구현 요약: Foldseek 연동 검증

**브랜치:** `replace-ncurses`  
**완료일:** 2026-04-27

### 목표

Phase 1–4에서 구현한 `structty_render` 정적 라이브러리를 Foldseek가 실제로 임베딩할 수 있는지 검증하고, 두 가지 CMake 연동 방식(add_subdirectory / find_package)을 모두 지원.

### 변경 파일

**`example/render_test.cpp`** (신규):
- `structty_render.h` + `AnsiOutput.hpp`만 include하는 독립 바이너리.
- Gemmi로 `.cif`/`.pdb` 파일을 직접 파싱해 Cα 원자 추출.
- `center_and_scale()`: Screen::normalize_proteins()와 동일 로직으로 원점 중심화 + 2.0 스케일.
- z 범위 자동 계산 후 `Renderer::set_depth_params()` 호출로 깊이 안개 정확도 확보.
- 사용법: `render_test <file> [mode] [width] [height]`

**`CMakeLists.txt`** (수정):
- `option(STRUCTTY_BUILD_APP "Build StrucTTY app" ON)` 추가 — OFF 시 Curses 미요구, Foldseek의 `add_subdirectory` 임베딩에 최적.
- `option(BUILD_RENDER_TEST "Build render_test" OFF)` 추가 — ON 시 `render_test` 바이너리 빌드.
- `find_package(Curses)` 및 StrucTTY 실행 파일 타겟을 `STRUCTTY_BUILD_APP` 조건으로 감쌈.

**`src/render/CMakeLists.txt`** (수정):
- `target_include_directories`를 generator expression으로 변경 (`$<BUILD_INTERFACE:...>` / `$<INSTALL_INTERFACE:...>`) — 빌드 트리와 설치 트리 모두 정확한 include 경로 보장.
- install 규칙 추가: `structty_render.a`, 공개 헤더 4종(`structty_render.h`, `Renderer.hpp`, `AnsiOutput.hpp`, `RenderPoint.hpp`), CMake 패키지 파일.

**`cmake/structty_renderConfig.cmake`** (신규):
- `find_package(structty_render)` 지원용 CMake config 파일.
- `find_dependency(ZLIB)` + `structty_renderTargets.cmake` include.
- 소비자는 Gemmi 심볼 해결을 위해 `gemmi::gemmi_cpp`도 별도 링크 필요 (Foldseek는 이미 링크 중).

### 설계 결정

| 항목 | 결정 | 이유 |
|------|------|------|
| add_subdirectory 지원 | `STRUCTTY_BUILD_APP=OFF` 옵션 | ncurses 없는 환경에서 render lib만 빌드 가능하게 |
| find_package 지원 | install 규칙 + Config 파일 | 시스템 설치 후 `find_package` 사용 케이스 지원 |
| render_test의 depth param | `set_depth_params()` 직접 호출 | 편의 함수 `render_to_stdout()`은 depth 설정 미지원; 테스트는 정확한 깊이 안개 필요 |
| render_test의 정규화 | 원점 중심화 + 스케일 2.0 | Screen::normalize_proteins() 단일 구조 경로와 동일 로직 |

### 빌드 및 검증 결과

```
# 기본 빌드 (StrucTTY 앱)
[100%] Built target StrucTTY  ← Zero new warnings

# render_test 활성화
cmake -DBUILD_RENDER_TEST=ON ...
[100%] Built target render_test  ← Zero new warnings

# 7가지 색상 모드 ANSI 출력 확인
render_test example/1CJK-assembly1.cif protein 60 12  → ANSI 시퀀스 정상 출력
render_test example/1CJK-assembly1.cif chain 60 12   → 체인별 색상 정상
(나머지 5종 동일 확인)

# 회귀 없음 확인
StrucTTY 1_1CRN.cif -m chain -s  → 기존 동작과 동일
```

*최종 업데이트: 2026-04-27*

---

## 9. Foldseek 연동 가이드

### 9.1 전제 조건

- Foldseek은 이미 Gemmi를 빌드 시스템에 포함하고 있어야 함 (structty_render의 정적 라이브러리가 Gemmi 심볼을 참조하므로 링크 단계에서 필요)
- C++17 이상
- Gemmi, ZLIB는 Foldseek 측에서 이미 제공됨을 가정

### 9.2 방식 A: add_subdirectory (권장)

StrucTTY 저장소를 Foldseek 저장소 안에 서브모듈 또는 복사본으로 포함하는 방식.

```bash
# Foldseek 저장소 루트에서
git submodule add https://github.com/steineggerlab/StrucTTY.git lib/structty
```

Foldseek의 `CMakeLists.txt`에 추가:

```cmake
# ncurses 없이 render 라이브러리만 빌드
set(STRUCTTY_BUILD_APP OFF CACHE BOOL "" FORCE)
add_subdirectory(lib/structty)

# Foldseek 타겟에 링크
target_link_libraries(foldseek
    PRIVATE
        structty_render
        gemmi::gemmi_cpp   # structty_render가 참조하는 Gemmi 심볼 해결
)
```

`STRUCTTY_BUILD_APP OFF`를 설정하면 ncurses `find_package`가 실행되지 않으므로 ncurses가 없는 빌드 환경에서도 안전하다.

### 9.3 방식 B: find_package (설치 후)

StrucTTY를 시스템에 설치한 뒤 `find_package`로 찾는 방식.

```bash
# StrucTTY 빌드 및 설치
cmake -B build -DSTRUCTTY_BUILD_APP=OFF -DCMAKE_INSTALL_PREFIX=/opt/structty
cmake --build build --target structty_render
cmake --install build --component structty_render
```

Foldseek의 `CMakeLists.txt`에 추가:

```cmake
find_package(structty_render REQUIRED
    PATHS /opt/structty/lib/cmake/structty_render)

target_link_libraries(foldseek
    PRIVATE
        structty_render::structty_render
        gemmi::gemmi_cpp
)
```

### 9.4 Foldseek 코드에서의 사용법

```cpp
#include "structty_render.h"  // 또는 <structty_render.h>

// 1. Foldseek의 내부 구조 데이터 → RenderAtom 변환
std::vector<RenderAtom> atoms;
for (const auto& residue : foldseek_hit.residues) {
    RenderAtom ra{};
    ra.x              = residue.ca_x;
    ra.y              = residue.ca_y;
    ra.z              = residue.ca_z;
    ra.structure      = 'x';         // 'x'=coil, 'H'=helix, 'S'=sheet
    ra.bfactor        = residue.plddt;
    ra.chain_id       = residue.chain_name;
    ra.protein_index  = 0;
    ra.pan_x = ra.pan_y = 0.0f;
    ra.conservation_score = -1.0f;
    ra.is_interface = ra.is_aligned = false;
    atoms.push_back(ra);
}

// 2. 원점 중심화 + 스케일 정규화 (필수)
//    → example/render_test.cpp의 center_and_scale() 참고
center_and_scale(atoms);

// 3. 렌더링 (터미널 크기를 width/height로 전달)
// 간단한 경우: 편의 함수 사용
structty::render_to_stdout(atoms, term_width, term_height, "protein");

// 깊이 안개 정밀도가 필요한 경우: Renderer 직접 사용
Renderer r(term_width, term_height, "chain", false);
r.set_depth_params(focal, zoom, min_z, max_z);  // example/render_test.cpp 참고
r.render(atoms);
AnsiOutput::print_to_stdout(r.get_pixels(), r.get_logical_width(), r.get_logical_height());
```

### 9.5 RenderAtom 필드 설명

| 필드 | 타입 | 필수 | 설명 |
|------|------|------|------|
| `x, y, z` | `float` | ✅ | Cα 3D 좌표 (정규화 후) |
| `structure` | `char` | ✅ | `'x'`=coil, `'H'`=helix, `'S'`=sheet |
| `bfactor` | `float` | plddt 모드 | pLDDT 점수 (0–100). 다른 모드에서는 무시 |
| `chain_id` | `std::string` | chain 모드 | 체인 이름. 다른 모드에서는 무시 |
| `protein_index` | `int` | ✅ | 0-based 단백질 인덱스 (단일 구조면 항상 0) |
| `pan_x, pan_y` | `float` | ✅ | 화면 오프셋 (단일 구조면 0.0) |
| `is_interface` | `bool` | interface 모드 | 인터페이스 잔기 여부 |
| `is_aligned` | `bool` | aligned 모드 | 정렬 잔기 여부 |
| `conservation_score` | `float` | conservation 모드 | 0.0–1.0, 미설정 시 -1.0 |
| `residue_number` | `int` | 선택 | 잔기 번호 (호버 정보용) |
| `residue_name[4]` | `char[]` | 선택 | 3-letter code + null (호버 정보용) |

### 9.6 색상 모드별 필요 데이터

| 모드 | 필요 필드 | 용도 |
|------|----------|------|
| `protein` | 없음 (protein_index만) | 구조별 9색 순환 |
| `chain` | `chain_id` | 체인별 15색 |
| `rainbow` | 없음 | N→C 방향 색상 그라디언트 |
| `plddt` | `bfactor` (pLDDT 0–100) | AlphaFold 신뢰도 4밴드 |
| `interface` | `is_interface` | 인터페이스 잔기 강조 |
| `conservation` | `conservation_score` (0–1) | 보존도 10단계 |
| `aligned` | `is_aligned` | 정렬 잔기 강조 |

### 9.7 좌표 정규화 유틸리티 (center_and_scale)

`example/render_test.cpp`에 구현된 정규화 함수를 Foldseek 코드에 복사해서 사용:

```cpp
static void center_and_scale(std::vector<RenderAtom>& atoms) {
    float min_x = std::numeric_limits<float>::max(), max_x = -min_x;
    float min_y = min_x, max_y = max_x;
    float min_z = min_x, max_z = max_x;
    for (const auto& a : atoms) {
        min_x = std::min(min_x, a.x); max_x = std::max(max_x, a.x);
        min_y = std::min(min_y, a.y); max_y = std::max(max_y, a.y);
        min_z = std::min(min_z, a.z); max_z = std::max(max_z, a.z);
    }
    float cx = 0.5f*(min_x+max_x), cy = 0.5f*(min_y+max_y), cz = 0.5f*(min_z+max_z);
    float ext = std::max({max_x-min_x, max_y-min_y, max_z-min_z});
    float scale = (ext > 0.0f) ? (2.0f / ext) : 1.0f;
    for (auto& a : atoms) {
        a.x = (a.x-cx)*scale;
        a.y = (a.y-cy)*scale;
        a.z = (a.z-cz)*scale;
    }
}
```
