# StrucTTY 프로젝트 리서치

## 1. 프로젝트 개요

**StrucTTY**는 C++로 작성된 터미널 기반 단백질 구조 시각화 도구이다. 최대 6개의 단백질 구조를 터미널에서 동시에 비교·시각화할 수 있으며, 회전/이동/체인 선택/정렬 시각화 등의 인터랙티브 기능을 지원한다. Foldseek 통합 및 확장성을 염두에 두고 설계되었다.

**핵심 특징:**
- 경량·고속 터미널 네이티브 GUI
- 최대 6개 단백질 동시 시각화
- 독립/동기화 구조 제어
- 인터페이스 검출, 보존성 스코어링, 정렬 시각화 등 고급 분자생물학 기능
- 크로스 플랫폼: Linux, macOS, Windows (Git Bash/MinTTY)

---

## 2. 프로젝트 구조

```
/home/user/StrucTTY/
├── src/                         # 메인 소스 코드 (C++)
│   ├── structty.cpp            # 메인 진입점 (221줄)
│   ├── structure/              # 단백질 파싱 및 처리
│   │   ├── Atom.hpp            # 원자 좌표 데이터 구조
│   │   ├── Protein.hpp/.cpp    # Protein 클래스 (761줄)
│   │   ├── Parameters.hpp/.cpp  # CLI 인자 파싱 (155줄)
│   │   ├── StructureMaker.hpp/.cpp # 이차 구조 시각화 (323줄)
│   │   ├── SSPredictor.hpp/.cpp # 이차 구조 예측
│   │   ├── FoldseekParser.hpp/.cpp # Foldseek 결과 파일 파싱 (181줄)
│   │   ├── FoldMasonParser.hpp/.cpp # FoldMason MSA 파싱 (295줄)
│   │   ├── MSAParser.hpp/.cpp  # MSA 파일 파싱 (127줄)
│   │   └── PDBDownloader.hpp/.cpp # Foldseek 히트 자동 다운로드 (397줄)
│   ├── visualization/          # 렌더링 및 UI
│   │   ├── Screen.hpp/.cpp     # 메인 렌더/입력 루프 (1745줄)
│   │   ├── Camera.hpp/.cpp     # 스크린샷 기능 (87줄)
│   │   ├── Panel.hpp/.cpp      # 정보 패널 표시 (654줄)
│   │   ├── RenderPoint.hpp     # 화면 픽셀 데이터 구조
│   │   ├── Palette.hpp         # 색상 체계 및 xterm256 색상
│   │   └── Curses.hpp          # ncurses 래퍼
│   └── utils/
│       ├── Common.hpp          # 타임스탬프 유틸리티
│       └── Benchmark.hpp       # 성능 벤치마킹
├── lib/                         # 서드파티 라이브러리
│   ├── gemmi/                  # Gemmi PDB/mmCIF 파서 (서브모듈)
│   └── lodepng/                # PNG 이미지 라이브러리 (서브모듈)
├── example/                     # 예제 PDB/mmCIF 파일 및 체인파일
├── benchmark/                   # 벤치마크 스크립트
├── CMakeLists.txt              # 빌드 설정
├── README.md                    # 사용자 문서
├── plan.md                      # 구현 로드맵 (68KB)
├── research_*.md               # Foldseek/FoldMason 관련 리서치 노트
└── structty-improvements.patch # 최근 구조 개선 패치
```

**총 소스 코드량:** 약 5,101줄의 C++ (헤더, 테스트, 라이브러리 제외)

---

## 3. 기술 스택

| 항목 | 상세 |
|------|------|
| **언어** | C++17 (GCC ≥ 7.1 또는 Clang ≥ 5.0) |
| **빌드 시스템** | CMake ≥ 3.15 |
| **PDB/mmCIF 파싱** | Gemmi (MPL-2.0) |
| **PNG 출력** | LodePNG (zlib) |
| **터미널 UI** | ncurses |
| **기타** | STL 컨테이너, filesystem, threading |

외부 패키지 매니저 없이 시스템 `curl`/`wget`을 `popen()`으로 호출하여 PDB 다운로드를 처리한다.

---

## 4. 아키텍처 및 설계 패턴

### 핵심 설계 원칙

1. **관심사 분리:** 구조 파싱(Protein) vs 시각화(Screen/Panel/Camera)
2. **이중 좌표 시스템:**
   - `init_atoms`: 원본 PDB Å 좌표 (정렬 계산용)
   - `screen_atoms`: 정규화된 뷰포트 좌표 (렌더링용)
3. **다형적 원자 속성:** 단일 Atom 구조체에 선택적 필드를 통한 기능 확장
   - `structure` (H/S/x: 나선/시트/코일)
   - `bfactor` (pLDDT 신뢰도)
   - `is_interface`, `is_aligned` (bool 플래그)
   - `conservation_score` (float)
   - `residue_number`, `residue_name` (호버 정보)

### 데이터 흐름

```
PDB/mmCIF 파일
    ↓ (Gemmi 파서)
Protein::init_atoms (Å 좌표)
    ↓ (normalize_proteins)
Protein::screen_atoms (뷰포트 좌표)
    ↓ (apply_foldseek_transform / apply_msa_conservation)
RenderPoint[] (투영된 2D 픽셀)
    ↓ (assign_colors_to_points)
Screen (ncurses 캔버스)
```

### 색상 코딩 전략

Palette.hpp에 여러 색상 팔레트가 정의되어 있다:
- Rainbow, Chain colors, Protein colors (기본 9 + 흐린 6)
- Interface (마젠타/흐림), Aligned (밝은 초록/흐림)
- pLDDT (4단계 파란-주황 그라데이션)
- Conservation (10단계 파란-빨간 그라데이션)

---

## 5. 핵심 모듈 및 역할

| 모듈 | 역할 | 핵심 함수 |
|------|------|-----------|
| **Protein** | 구조 데이터 로드 및 변환 | `load_init_atoms()`, `compute_interface()`, `apply_ut_to_init_atoms()` |
| **Parameters** | CLI 인자 파싱 | `Parameters::Parameters(argc, argv)` |
| **StructureMaker** | 이차 구조 렌더링 | `calculate_ss_points()` (나선 지그재그, 베타 시트 리본, Catmull-Rom 스플라인) |
| **SSPredictor** | 거리/비틀림 각도 기반 구조 예측 | `run()`, `run_chain()` |
| **FoldseekParser** | Foldseek 결과 파일 파싱 (12/17/21/29열 형식) | `load()`, 형식 자동 감지 |
| **FoldMasonParser** | FoldMason JSON/FASTA MSA 파싱 | `load_json()`, `load_fasta()`, `compute_column_entropy()` |
| **MSAParser** | MSA 파일 파싱 (FASTA/A3M) | `load()`, `compute_conservation()` |
| **PDBDownloader** | 다양한 DB에서 Foldseek 히트 자동 다운로드 | DB 타입 감지 (PDB, AlphaFold, GMGCL, TED, BFMD 등) |
| **Screen** | 메인 렌더링 및 인터랙션 | `draw_screen()`, `handle_input()`, `compute_interface_all()` |
| **Panel** | 정보 패널 UI | `draw_panel()`, `set_foldseek_hit_info()`, `update_hover_info()` |
| **Camera** | 스크린샷 캡처 | `screenshot()`, PNG 이미지 생성 |

---

## 6. 빌드 시스템

```bash
git clone --recurse-submodules https://github.com/steineggerlab/StrucTTY.git
cd StrucTTY && mkdir build && cd build
cmake ..
make -j 10
# 출력: build/StrucTTY
```

- macOS: Homebrew ncurses 수동 설치 필요 (CMake 플래그로 include/lib 경로 지정)
- Windows: MSYS2 MinGW64 GitHub Actions 워크플로우 (`.github/workflows/windows-mingw64.yml`)

---

## 7. 테스트 현황

**공식 테스트 스위트 없음.** 현재 테스트 방식:
- 벤치마크 모드(`-b` 플래그): 200회 워밍업 + 2000회 측정 이벤트
- `example/` 디렉토리의 예제 파일 (10개 PDB/mmCIF)
- 벤치마크 스크립트를 통한 성능 측정

---

## 8. 진입점 및 메인 흐름

**`src/structty.cpp` (main 함수):**

1. ncurses 초기화 및 CLI 인자 파싱 (Parameters 클래스)
2. 입력 파일별 단백질 로드 → `screen.set_protein()`
3. 변환 적용:
   - UT 매트릭스 (`-ut` 플래그)
   - 뷰포트 정규화
   - 인터페이스/정렬 영역 계산
4. 고급 기능 로드 (해당 시):
   - Foldseek 히트 (`-fs`): 히트 목록 로드 → 첫 히트 자동 로드
   - MSA 보존성 (`--msa`): MSA 파싱 → 스코어 계산
   - FoldMason (`--foldmason`): JSON/FASTA MSA → 중첩 → 보존성
5. 벤치마크 모드 또는 인터랙티브 루프 실행

**키보드 컨트롤:**
- `0-6`: 제어할 단백질 선택
- `W/A/S/D`: 이동 (±Y, ±X)
- `X/Y/Z`: 축 회전
- `R/F`: 줌 인/아웃
- `N/P`: Foldseek 히트 네비게이션 (다음/이전)
- `Q`: 종료

---

## 9. 개발 현황

### 완료된 기능

- **Feature 1** — 인터페이스 영역 검출 (CA-CA 거리 < 8Å 임계값)
- **Feature 2** — pLDDT 신뢰도 색상화 (B-factor 시각화)
- **Feature 4a** — UTMatrix 정렬 색상화 (`-ut` 플래그)
- **Feature 4b** — Foldseek 기반 정렬 색상화 (`-fs` + 좌표 공간 수정)
- **Feature 5** — MSA 보존성 스코어링 (엔트로피 기반)
- **Feature 6** — 마우스 호버 잔기 정보 패널
- **Feature 7** — 이차 구조 렌더링 개선 (나선 지그재그, 베타 시트 리본, Catmull-Rom 스플라인)
- **Feature 3 (Mode A)** — Foldseek 단일 파일 히트 네비게이션

### 진행 중 / 부분 구현

- **Feature 3 (Mode B)** — 다중 파일 Foldseek 네비게이션
- **Feature 8 (Mode B)** — 다중 파일 FoldMason

### 수정된 버그

- **BUG-A** (커밋 fdb2024): `-m aligned` 모드 좌표 공간 불일치 수정
- **BUG-B**: BFMD DBType 추가 (미지원 DB 다운로드 URL 메시지)

### 미해결 이슈

- 다중 파일 모드(Foldseek/FoldMason) 미완성
- 회귀 테스트 스위트 부재

---

## 10. 코드 컨벤션

| 항목 | 규칙 |
|------|------|
| 클래스명 | PascalCase (Protein, StructureMaker, FoldseekParser) |
| 함수명 | snake_case (load_init_atoms, compute_interface) |
| 멤버 변수 | 파라미터에 trailing underscore, 필드에 접두사 없음 |
| 상수 | UPPER_CASE (RAINBOW, INTERFACE_COLOR) |
| 헤더 전용 유틸 | Benchmark.hpp, Palette.hpp, Common.hpp |
| 에러 처리 | CLI 파싱은 예외 기반, 파일 연산은 bool 반환 + out-parameter |

**성능 고려사항:**
- 고빈도 객체에서 `std::string` 회피 (RenderPoint는 `char[4]` 사용)
- 이중 버퍼 설계: init_atoms (Å 공간) vs screen_atoms (뷰포트 공간)
- Catmull-Rom 스플라인 보간으로 부드러운 백본 렌더링

---

## 11. 설정 파일 요약

| 파일 | 용도 |
|------|------|
| `CMakeLists.txt` | 빌드 설정 (C++17, Curses, 서브모듈) |
| `.gitignore` | build/, .vscode/, benchmark/ 제외 |
| `.gitmodules` | gemmi & lodepng 서브모듈 링크 |
| `plan.md` | 8-기능 구현 로드맵 + 3 버그 리포트 |
| `research_align.md` | `-m aligned` 좌표 공간 불일치 근본 원인 분석 |
| `research_foldmason.md` | FoldMason JSON/FASTA 형식 명세 |
| `THIRD_PARTY_NOTICES.md` | 라이선스 표기 (Gemmi: MPL-2.0, LodePNG: zlib) |

---

## 12. 요약

StrucTTY는 파싱·변환·렌더링의 깔끔한 관심사 분리와 확장 가능한 기능 설계를 결합한, 잘 설계된 단백질 구조 시각화 도구이다. 8-기능 로드맵의 약 90%가 완료되었으며, Foldseek 통합을 위한 견고한 에러 처리와 고급 분자생물학 시각화 모드를 갖추고 있다. 단일/이중 구조 시각화에 대해 프로덕션 수준이며, 다중 파일 정렬 시나리오를 적극적으로 개발 중이다.
