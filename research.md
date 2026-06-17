# Project Research

## 1. 프로젝트 목적 및 전체 구조 요약

**StrucTTY**는 C++17로 작성된 경량 터미널 기반 단백질 구조 뷰어다. 3D 단백질 구조를 **Unicode Braille 서브픽셀 렌더링**으로 터미널에 직접 출력한다 (각 문자 셀 = 2×4 논리 픽셀 → 일반 문자 대비 8× 해상도). Foldseek / FoldMason과 통합되어 구조 검색 결과 탐색과 다중 구조 정렬을 지원한다.

```
StrucTTY/
├── src/
│   ├── main.cpp                  # 진입점: Parameters → RunOptions → structty::run()
│   ├── structty.cpp              # structty::run() 구현: Scene 셋업 + 이벤트 루프
│   ├── structure/                # 데이터 파싱·처리 레이어
│   │   ├── Atom.hpp              # 단일 Cα 잔기 데이터 (좌표, SS, pLDDT, 플래그들)
│   │   ├── Protein.{hpp,cpp}     # 단백질 구조 (init_atoms + screen_atoms, gemmi 파싱)
│   │   ├── StructureMaker.{hpp,cpp}  # 2차구조 3D 기하 (헬릭스 실린더, 베타 리본)
│   │   ├── SSPredictor.{hpp,cpp}     # CA-거리+이면각 기반 SS 예측
│   │   ├── FoldseekParser.{hpp,cpp}  # .m8 파일 파싱 (12/17/21/29 컬럼)
│   │   ├── FoldseekDBReader.{hpp,cpp}# Foldseek DB 직접 읽기 (선택적 스캔)
│   │   ├── FoldMasonParser.{hpp,cpp} # FoldMason JSON/FASTA 파싱 + Kabsch superposition
│   │   ├── MSAParser.{hpp,cpp}       # FASTA/A3M MSA → Shannon entropy conservation
│   │   ├── MultimerReportParser.{hpp,cpp} # _report 14-컬럼 TSV (complex 단위)
│   │   ├── PDBDownloader.{hpp,cpp}   # target ID 패턴 감지 + curl/wget 다운로드
│   │   └── Parameters.{hpp,cpp}      # CLI arg 파싱
│   ├── visualization/            # 렌더링 코디네이터·시각화 레이어
│   │   ├── Screen.{hpp,cpp}      # 중앙 오케스트레이터 (Protein[], Renderer, Panel, Camera)
│   │   ├── Camera.{hpp,cpp}      # PNG 스크린샷 (lodepng)
│   │   ├── Panel.{hpp,cpp}       # 정보 패널 ANSI 출력 (단백질 목록, 잔기 hover, Foldseek info)
│   │   ├── Palette.hpp           # xterm-256 색상 팔레트 + depth fog 색상 세트
│   │   └── RenderPoint.hpp       # 논리 픽셀 데이터 구조체
│   ├── render/                   # ncurses-free 렌더링 라이브러리 (structty_render)
│   │   ├── Renderer.{hpp,cpp}    # RenderAtom → 2D 논리 픽셀 버퍼 (z-buffer 해소)
│   │   └── AnsiOutput.{hpp,cpp}  # 논리 픽셀 → Braille Unicode ANSI 문자열
│   └── utils/
│       ├── Terminal.{hpp,cpp}    # Raw mode, 크기 조회, 키/마우스 입력, ANSI 화면 제어
│       ├── Common.hpp            # current_timestamp() 유틸
│       └── Benchmark.hpp        # FPS/latency CSV 로그 (--benchmark 모드)
├── include/
│   ├── structty.h                # 공개 API: structty::RunOptions, structty::run()
│   └── structty_render.h         # 독립 렌더 API: render_to_ansi(), render_to_stdout()
├── lib/
│   ├── gemmi/                    # git submodule: mmCIF/PDB 파싱 (헤더-온리, MPL-2.0)
│   └── lodepng/                  # git submodule: PNG 인코더 (단일 TU, zlib)
├── cmake/
│   └── structty_renderConfig.cmake  # find_package(structty_render) 지원
├── example/                      # 샘플 CIF/PDB, Foldseek/FoldMason 결과물
└── .github/
    └── workflows/windows-mingw64.yml  # Windows MSYS2/mingw64 CI (수동 트리거)
```

---

## 2. 핵심 모듈/컴포넌트와 각각의 역할

### `structty::run()` — `src/structty.cpp:22`
최상위 진입점. Screen 생성 → 씬 셋업(단백질 로드, 파서 실행, transform 적용) → 이벤트 루프. 세 가지 분기:
- 일반 모드: CLI 파일 순차 로드 → Foldseek/MSA/FoldMason 선택 적용
- `query_from_db`: Foldseek query tmp DB에서 쿼리 구조 직접 읽기 (멀티 쿼리 내비게이션)
- `multimer_report`: `_report` 14-컬럼 TSV → complex 단위 superposition

### `Screen` — `src/visualization/Screen.{hpp,cpp}`
중앙 오케스트레이터. 소유하는 객체:
- `std::vector<Protein*> data` — 최대 9개 구조
- `Renderer renderer_` — 논리 픽셀 버퍼
- `Camera* camera` — PNG 스크린샷
- `Panel* panel` — 정보 패널
- `FoldseekDBReader fs_db_reader_` — target DB
- `FoldseekDBReader query_db_reader_` — query tmp DB (Step 4)
- `std::unique_ptr<FoldMasonParser> foldmason_parser`

주요 메서드: `draw_screen()`, `handle_input()`, `set_protein()`, `load_next_hit()`, `apply_foldseek_transform()`, `apply_foldmason_superposition()`, `normalize_proteins()`, `set_tmatrix()`, `compute_interface_all()`, `switch_query()`, `activate_multimer_query()`

### `Protein` — `src/structure/Protein.{hpp,cpp}`
단일 구조 래퍼. gemmi로 CIF/PDB 파일 파싱. 두 개의 원자 맵 유지:
- `init_atoms`: 원본 좌표 (Kabsch/UT transform 기준)
- `screen_atoms`: 회전·이동 적용 후 렌더링 좌표

메서드: `do_rotation()`, `do_shift()`, `do_scale()`, `compute_interface()`, `apply_conservation_scores()`, `compute_aligned_regions_from_aln()`, `compute_aligned_regions_nn()`, `load_from_ca()` (Foldseek DB Cα only 로드)

### `Atom` — `src/structure/Atom.hpp`
단일 Cα 잔기. 필드:
- `x, y, z: float` — 좌표
- `structure: char` — 'x'(코일), 'H'(헬릭스), 'S'(베타)
- `bfactor: float` — pLDDT (B-factor 컬럼)
- `is_interface: bool` — 계면 잔기 여부
- `is_aligned: bool` — 정렬 구역 여부
- `conservation_score: float` — MSA entropy score
- `residue_number: int`, `residue_name: string` — hover 표시용

### `Renderer` — `src/render/Renderer.{hpp,cpp}`
`RenderAtom` 배열 → 논리 픽셀 버퍼(width×2, height×4). 파이프라인:
1. `project_and_fill()`: 투영(FOV=90°) + 선 그리기(draw_line_impl) + color/depth_band 할당
2. `zbuffer_resolve()`: 깊이별 최전방 픽셀 선택 → `logical_pixels_` 에 기록

### `AnsiOutput` — `src/render/AnsiOutput.{hpp,cpp}`
논리 픽셀 버퍼 → Braille Unicode ANSI 문자열. 각 2×4 서브픽셀 블록을 U+2800+bitmask로 인코딩. xterm-256 색상 escape 적용.

### `RenderAtom` — `src/render/Renderer.hpp:6`
Screen → Renderer 전달 데이터 구조체. `Atom` + `protein_index`, `pan_x/y`, `chain_id`.

### `RenderPoint` — `src/visualization/RenderPoint.hpp`
논리 픽셀 1개. `x, y, depth, color_id, chainID, structure, depth_band`, 잔기 정보 필드.

### `FoldseekParser` — `src/structure/FoldseekParser.{hpp,cpp}`
Foldseek `.m8` 파일 파싱. 컬럼 수(12/17/21/29) 자동 감지. `FoldseekHit` 벡터 생성. 21컬럼(alis) 포맷: `alns`(정렬된 CA 좌표 리스트), `tseq`, `taxid` 포함.

### `FoldseekDBReader` — `src/structure/FoldseekDBReader.{hpp,cpp}`
Foldseek DB 직접 읽기. `prepare()` 에서 hit ID 세트만 인덱싱 (`.lookup` + `_ca.index` 순차 스캔, ~152KB/1000 hits). `read_entry()` 로 Cα 좌표 + AA 서열 조회.

### `FoldMasonParser` — `src/structure/FoldMasonParser.{hpp,cpp}`
FoldMason JSON(Cα 좌표 포함) 또는 FASTA(서열만) 파싱. `build_aligned_pairs()`: gap-free 공통 열 잔기 쌍 추출(Kabsch superposition 입력). `compute_column_entropy()`: 열별 Shannon entropy.

### `MSAParser` — `src/structure/MSAParser.{hpp,cpp}`
FASTA/A3M MSA 파싱 → Shannon entropy 기반 conservation score 계산. A3M의 소문자 삽입(insertion) 열 자동 제거.

### `MultimerReportParser` — `src/structure/MultimerReportParser.{hpp,cpp}`
Foldseek `_report` 14-컬럼 TSV 파서. `MultimerHit` 벡터: complex 이름, 체인 ID 목록, 9-float U + 3-float T.

### `PDBDownloader` — `src/structure/PDBDownloader.{hpp,cpp}`
Foldseek target ID 패턴 감지 → DB 종류 판별 (PDB/AlphaFold/ESMAtlas30/CATH50/BFVD/TED/BFMD/GMGCL/Unknown). `resolve_target_file()`: db_path 탐색 → 캐시 → curl/wget 다운로드. 캐시 경로: `~/.cache/structty/pdb/`.

### `SSPredictor` — `src/structure/SSPredictor.{hpp,cpp}`
Cα 거리(d13, d14) + 이면각 기반 투표(vote) → 헬릭스/베타 레이블. smooth + 최소 길이 필터. CA-only 구조에도 적용 가능.

### `StructureMaker` — `src/structure/StructureMaker.{hpp,cpp}`
2차구조 3D 렌더링 기하 계산. 헬릭스: 축 계산 + 실린더 점 생성. 베타: ribbon 오프셋(±4 steps, step=0.28).

### `Panel` — `src/visualization/Panel.{hpp,cpp}`
ANSI 정보 패널 (우측 컬럼). 섹션: 단백질 목록(파일명+체인+잔기수), Foldseek hit info(target/evalue/prob/lddt), Residue hover info(체인/잔기명/번호/SS/pLDDT/conservation), FoldMason info. `draw_hover_section()`: 마우스 이동 시 패널 부분 갱신만 수행.

### `Camera` — `src/visualization/Camera.{hpp,cpp}`
`screenshot()`: RenderPoint 벡터 → PNG 저장 (lodepng). `renderPoint2image()`: color_id → RGBA 변환(Palette::ID2RGBA). 저장 경로: `~/[camera_dir]/`.

### `Terminal` — `src/utils/Terminal.{hpp,cpp}`
Raw mode 관리(SIGTERM/SIGINT atexit). `read_key()`: blocking. `read_mouse()`: SGR 1006 마우스 이벤트. `get_size()`: ioctl TIOCGWINSZ.

### `Palette` — `src/visualization/Palette.hpp`
xterm-256 팔레트. depth fog를 위한 3세트(near/mid/far) × 각 모드별 색상 배열. `palette_to_xterm256_fg(int cid)`: color_id → xterm256 번호. `palette_to_ansi_fg_str()`: ANSI escape 문자열.

### `Parameters` — `src/structure/Parameters.{hpp,cpp}`
CLI arg 파싱. `-m/--mode`, `-c/--chains`, `-s/--structure`, `-ut/--utmatrix`, `--msa`, `-fs/--foldseek`, `--db`, `--db-path`, `-fm/--foldmason`, `-n/--nopanel`, `-b/--benchmark`, `--report-format`, `--query-db`. 최대 9개 입력 파일.

---

## 3. 빌드 시스템 및 의존성

| 항목 | 내용 |
|------|------|
| 빌드 도구 | CMake ≥ 3.15 |
| 언어 표준 | C++17 (required) |
| 기본 빌드 타입 | Release |
| 지원 플랫폼 | Linux, macOS (공식); Windows MSYS2/mingw64 (CI 있음, 수동 트리거) |

**빌드 타겟:**
- `structty_render` (STATIC): `Renderer.cpp` + `AnsiOutput.cpp`. ncurses 의존 없음. 설치 가능(find_package 지원).
- `structty` (STATIC): 전체 라이브러리 (`structty.cpp` + `structure/*.cpp` + `visualization/*.cpp` + `utils/*.cpp`). Foldseek가 `add_subdirectory`로 임베딩 가능.
- `StrucTTY` (EXECUTABLE): `main.cpp` → `structty` 링크. `STRUCTTY_BUILD_APP=ON`(기본값)일 때만 빌드.
- `render_test` (EXECUTABLE): `BUILD_RENDER_TEST=ON`일 때만 빌드.

**외부 의존성:**

| 라이브러리 | 버전 | 방식 | 라이선스 | 용도 |
|-----------|------|------|---------|------|
| Gemmi | 0.5.1 (git submodule) | INTERFACE (헤더-온리) | MPL-2.0 | mmCIF/PDB 파싱 |
| LodePNG | HEAD (git submodule) | STATIC | zlib | PNG 인코더 |
| ZLIB | 시스템 (선택적) | CMake find_package | zlib | Gemmi gz 지원; 없으면 Gemmi 내장 zlib 사용 |

**빌드 명령:**
```bash
git clone --recurse-submodules https://github.com/steineggerlab/StrucTTY.git
cmake -B build
cmake --build build -j$(nproc)
# 바이너리: build/StrucTTY
```

---

## 4. 주요 데이터 흐름

### 일반 모드 (단일/다중 파일)
```
CLI args
  → Parameters::Parameters()
  → structty::RunOptions
  → structty::run()
    → Terminal::enter_raw_mode()
    → Screen(width, height, mode)
    → Screen::set_protein(file, i)  [gemmi로 CIF/PDB 파싱]
      → Protein() → load_init_atoms() → (SSPredictor::run() if show_structure)
    → Screen::set_tmatrix()         [회전행렬 배열 할당]
    → Screen::normalize_proteins()  [centroid → 원점 이동, 스케일 정규화]
    → [옵션] compute_interface_all() / apply_msa_conservation() / apply_foldmason_superposition()
    → 이벤트 루프:
        Screen::draw_screen()
          → to_render_atoms()       [screen_atoms → RenderAtom[]]
          → Renderer::render()
              → project_and_fill()  [투영 + 선 그리기 + 색상/depth_band 할당]
              → zbuffer_resolve()   [깊이 해소 → logical_pixels_]
          → AnsiOutput::print_to_stdout()  [Braille ANSI 출력]
          → Panel::draw_panel()     [우측 패널]
        Screen::handle_input()      [키/마우스 → 회전·이동·줌·구조 전환]
    → Terminal::exit_raw_mode()
```

### Foldseek hit 탐색 흐름
```
.m8 파일 → FoldseekParser::load()
  → Screen::set_foldseek_hits()
  → Screen::prepare_foldseek_db()   [hit ID 세트 → FoldseekDBReader::prepare()]
  → Screen::load_next_hit(+1)
      → FoldseekDBReader::read_entry() 또는 PDBDownloader::resolve_target_file()
      → Protein::load_from_ca() 또는 Protein(file)
      → Kabsch SVD (21컬럼 alis 포맷) 또는 U/T 직접 사용 (29컬럼)
      → Screen::apply_foldseek_transform()
      → Panel::set_foldseek_hit_info()
```

### FoldMason 슈퍼포지션 흐름
```
JSON/FASTA → FoldMasonParser::load_json() / load_fasta()
  → build_aligned_pairs(query_idx, target_idx)  [gap-free 잔기 쌍]
  → Kabsch SVD → U, T
  → Screen::apply_foldmason_superposition()
      → Protein::apply_ut_to_init_atoms()
      → compute_aligned_regions_from_aln()
  → [conservation] compute_column_entropy() → apply_msa_conservation()
```

### Multimer report 흐름
```
_report TSV → MultimerReportParser::load()
  → Screen::set_multimer_report()
      → FoldseekDBReader::open(query_db) + prepare(query_complexes)
      → activate_multimer_query(0)
          → load_chain_into_data() per query chain
          → normalize_complex()
          → load_multimer_hit(0)
              → load_chain_into_data() per target chain
              → transform_target_chain() [U,T Å → query frame]
```

---

## 5. 현재 구현된 기능 목록

코드에서 확인된 기능:

1. **Braille 서브픽셀 렌더링** — `AnsiOutput.cpp:14` (DOT_BITS 레이아웃), 각 셀 2×4 논리 픽셀
2. **다중 구조 (최대 9개)** — `Parameters.cpp:69` (`in_file.size() < 9`), `Screen.hpp:10` (`MAX_STRUCT_NUM = 9`)
3. **7가지 색상 모드** — `Parameters.cpp:43` 검증: protein/chain/rainbow/plddt/interface/conservation/aligned
4. **3-band depth fog** — `Renderer.cpp:41` (`compute_depth_band`: t<0.25=near, t<0.60=mid, else far), `Palette.hpp` near/far 색상 세트
5. **2차구조 시각화** — `SSPredictor.cpp` (CA 거리+이면각 투표), `StructureMaker.cpp` (헬릭스 실린더, 베타 리본)
6. **Foldseek .m8 내비게이션** — `FoldseekParser.hpp:12` (12/17/21/29 컬럼 지원)
7. **Foldseek DB 직접 읽기** — `FoldseekDBReader.hpp:20` (선택적 스캔)
8. **멀티 쿼리 내비게이션 (]/[)** — `Screen.hpp:49` (`set_query_nav`), `Screen.cpp` (`switch_query`)
9. **Multimer _report 지원** — `MultimerReportParser.hpp:5`, `Screen.hpp:61` (`set_multimer_report`)
10. **PDB 자동 다운로드** — `PDBDownloader.hpp:20` (PDB/AlphaFold/ESMAtlas/CATH/BFVD/TED 지원)
11. **FoldMason 슈퍼포지션** — `FoldMasonParser.hpp:28` (Kabsch SVD), JSON/FASTA 모두 지원
12. **MSA conservation** — `MSAParser.hpp:18` (Shannon entropy, FASTA/A3M)
13. **계면 감지 (Interface)** — `Protein.hpp:72` (CA-CA < 8Å)
14. **정렬 구역 표시 (Aligned)** — `Protein.hpp:81` (aln-string 또는 nearest-neighbor)
15. **체인 선택** — `Parameters.cpp:52` (`-c/--chains`), TSV 형식
16. **UT 행렬 외부 정렬** — `Screen.hpp:69` (`set_utmatrix`)
17. **PNG 스크린샷** — `Camera.hpp:29` (lodepng, `~/` 하위 저장)
18. **마우스 hover** — `Terminal.hpp:14` (SGR 1006), `Panel.hpp:76` (`set_hover_residue`), `Screen.hpp:118`
19. **벤치마크 모드** — `Benchmark.hpp` (FPS/latency CSV, 200 warmup + 2000 frames)
20. **패널 없는 모드** — `-n/--nopanel`
21. **structty_render 독립 설치** — `include/structty_render.h` (`render_to_ansi`, `render_to_stdout`)
22. **Foldseek `--view` 연동** — README: `foldseek easy-search ... --view 1`
23. **Multi-input Foldseek** — `structty.cpp:130` (복수 target 파일과 .m8 hit 매칭)

---

## 6. 미구현 또는 TODO 항목

코드 전체 스캔 결과 TODO/FIXME/HACK/UNIMPLEMENTED 마커는 **1개**만 발견:

- `src/structure/PDBDownloader.cpp:127` — `// AF-XXXXX-F\d+-model_v\d+` (정규식 설명 주석, 실제 미구현 아님)

**기타 알려진 제한사항 (코드 및 주석에서 확인):**
- GMGCL DB 타입은 다운로드 URL 없음 — `PDBDownloader.hpp:16` (`// 다운로드 URL 없음`)
- BFMD DB 타입은 다운로드 URL 없음, `--db-path` 필요 — `PDBDownloader.hpp:17`
- Windows 공식 지원 안 함 (CI는 있으나 플랫폼 지원 목록에서 제외): `README.md:49`
- `src/CMakeLists.txt`의 `structty_core` 타겟은 최상위 `CMakeLists.txt`와 달리 레거시 구조로 보임 (현재 사용 안 함)
- **테스트 스위트 없음** — 공식 단위/통합 테스트 없음 (`BUILD_RENDER_TEST`는 통합 테스트 바이너리 옵션이지만 `example/render_test.cpp` 하나뿐)
- `Atom::get_position()`의 `static float coords[3]` — thread-unsafe: `Atom.hpp:39`

---

## 7. 코드 컨벤션 및 네이밍 규칙

**클래스/구조체:** PascalCase — `Screen`, `Protein`, `Renderer`, `FoldseekParser`, `RenderPoint`

**함수/메서드:** snake_case 혼용 — `set_protein()`, `draw_screen()`, `load_next_hit()`, `compute_depth_band()`

**멤버 변수:** 두 패턴 혼용
- 일부: trailing `_` — `renderer_`, `depth_base_min_z_`, `query_ids_`, `multi_query_show_structure_`
- 일부: 없음 — `screen_width`, `screen_mode`, `zoom_level`, `foldseek_hits`

**네임스페이스:** lowercase — `namespace structty`, `namespace Terminal`, `namespace Palettes`, `namespace fs (= std::filesystem)`

**상수:** 매크로 또는 `constexpr` — `MAX_STRUCT_NUM = 9` (`Screen.cpp:10`), `FOV_ = 90.0f` (`Renderer.hpp:43`)

**주석 언어:** 한국어 광범위 사용. "기능 N" 번호 체계:
- 기능 1: interface 감지
- 기능 2: pLDDT (B-factor)
- 기능 3: Foldseek hit 내비게이션
- 기능 4: UTMatrix 정렬 / alignment string
- 기능 5: MSA conservation
- 기능 6: 마우스 hover
- 기능 7: depth fog
- 기능 8: FoldMason

**파일 구성:** 기능/관심사 단위 (`structure/`, `visualization/`, `render/`, `utils/`)

**헤더 가드:** `#pragma once` (통일)

**들여쓰기:** 4 스페이스

**오류 처리:** 대부분 `std::cerr` + `return false`; 생성자에서 예외 없음; Parameters만 `std::runtime_error` 사용

---

## 구현 완료: Parameters CLI 플래그 추가 — 2026-06-17

### 구현된 기능
- `--report-format` 플래그: `--foldseek`로 전달된 파일이 `_report` 14-컬럼 TSV임을 명시
- `--query-db <PATH>` 플래그: query complex Foldseek DB 경로 (`--report-format`과 함께 필수)
- `--report-format` 사용 시 `in_file` 없이도 실행 가능 (이전에는 무조건 입력 파일 필요)
- `--report-format` 단독 사용 시 오류 메시지 출력 (`--foldseek`와 `--query-db` 모두 필요)

### 변경된 파일
- `src/structure/Parameters.hpp` — `report_format`, `foldseek_query_db` private 멤버 + getter 2개 추가
- `src/structure/Parameters.cpp` — help 텍스트, 파싱 루프, 검증 로직 추가

### 계획 대비 변경사항
- 변경 없음

---

## 구현 완료: main.cpp RunOptions 연결 + Screen 숫자키 complex 단위 선택 — 2026-06-17

### 구현된 기능
- `main.cpp`: `opts.report_format` / `opts.foldseek_query_db` 연결 완료
- `Screen.cpp`: `mm_complex_range(complex_idx)` 헬퍼 추가 — complex_idx=0→query 체인 범위, 1→target 체인 범위 반환 (경계 방어 포함)
- `Screen.cpp`: 숫자키 핸들러 multimer 분기 추가 — 키 1→query complex(structNum=0), 키 2→target complex(structNum=1), 키 3+→무시
- `Screen.cpp`: A/D/W/S pan 핸들러 — multimer+structNum≥0 시 complex 전체 체인에 pan 적용
- `Screen.cpp`: X/Y/Z 회전 핸들러 — multimer+structNum≥0 시 complex 전체 체인에 rotation 적용
- `Screen.hpp`: `mm_complex_range` private 헬퍼 선언 추가

### 변경된 파일
- `src/main.cpp` — `opts.report_format`, `opts.foldseek_query_db` 연결 2줄 추가
- `src/visualization/Screen.hpp` — `mm_complex_range` 선언 추가
- `src/visualization/Screen.cpp` — `mm_complex_range` 구현, 숫자키/pan/rotation 핸들러 multimer 분기

### 계획 대비 변경사항
- 변경 없음
