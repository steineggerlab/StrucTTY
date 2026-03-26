# StrucTTY 프로젝트 리서치

---

## 1. 프로젝트 개요

**StrucTTY**는 C++17로 작성된 터미널 네이티브 단백질 구조 시각화 도구. 터미널에서 직접 단백질 3D 구조를 인터랙티브하게 시각화하고 비교할 수 있다. 최대 6개의 단백질을 동시에 표시하며, Foldseek/FoldMason 등 생물정보학 도구와 통합된다.

**핵심 특징**:
- 유니코드 Braille 문자로 서브픽셀 렌더링 (셀당 2×4 논리 픽셀)
- ncurses 기반 터미널 제어
- 3단계 깊이 안개(depth fog) 시스템
- 7가지 색상 모드 (protein, chain, rainbow, pLDDT, interface, conservation, aligned)
- Foldseek/FoldMason 결과 직접 로드 및 탐색

---

## 2. 디렉토리 구조

```
/home/user/StrucTTY/
├── src/
│   ├── structty.cpp                    # 메인 진입점
│   ├── structure/                      # 단백질 데이터 구조 & 파서
│   │   ├── Atom.hpp                    # 단일 원자 표현
│   │   ├── Protein.hpp/.cpp            # 단백질 클래스 (체인/원자 관리, 변환)
│   │   ├── Parameters.hpp/.cpp         # CLI 인자 파싱
│   │   ├── StructureMaker.hpp/.cpp     # 2차 구조 시각화 (나선/시트 기하)
│   │   ├── SSPredictor.hpp/.cpp        # 2차 구조 예측
│   │   ├── FoldseekParser.hpp/.cpp     # Foldseek .m8 파서
│   │   ├── FoldMasonParser.hpp/.cpp    # FoldMason JSON/FASTA MSA 파서
│   │   ├── MSAParser.hpp/.cpp          # MSA 파일 파서 (FASTA/A3M)
│   │   └── PDBDownloader.hpp/.cpp      # 다중 DB PDB 파일 다운로드
│   ├── visualization/                  # 렌더링 & 디스플레이
│   │   ├── Screen.hpp/.cpp             # 메인 스크린/렌더러
│   │   ├── Panel.hpp/.cpp              # 정보 패널
│   │   ├── Camera.hpp/.cpp             # 카메라 & 스크린샷
│   │   ├── Palette.hpp                 # 색상 정의 (250+ ncurses 색상 쌍)
│   │   └── RenderPoint.hpp             # 렌더링 포인트
│   └── utils/
│       ├── Common.hpp                  # 타임스탬프 유틸리티
│       ├── Curses.hpp                  # ncurses 헤더 감지
│       └── Benchmark.hpp               # 성능 측정
├── lib/                                # 외부 의존성
│   ├── gemmi/                          # Gemmi (MPL-2.0) - mmCIF/PDB 파싱
│   └── lodepng/                        # LodePNG (zlib) - PNG 인코딩
├── example/                            # 예제 데이터
│   ├── *.cif                           # mmCIF 구조 파일 (1NPL, 3HGM, 9N47, 1CJK 등)
│   ├── chainfile_9N47.tsv              # 체인 선택 예제
│   ├── foldseek_result/                # Foldseek .m8 출력 (여러 DB)
│   ├── foldmason_result/               # FoldMason JSON + FASTA
│   └── msa_result/                     # MSA 파일 예제
├── CMakeLists.txt                      # CMake 빌드 설정
├── README.md                           # 사용자 문서
└── THIRD_PARTY_NOTICES.md              # 라이선스 정보
```

---

## 3. 빌드 시스템

### 요구 사항
- C++17 (GCC ≥7.1 또는 Clang ≥5.0)
- CMake ≥3.15
- Linux/macOS (Windows는 Git Bash/MinTTY)

### 의존성
| 라이브러리 | 라이선스 | 용도 | 통합 방식 |
|-----------|---------|------|----------|
| **ncurses** | Public Domain | 터미널 제어/렌더링 | 시스템 패키지 (REQUIRED) |
| **gemmi** | MPL-2.0 | mmCIF/PDB 파싱 | `gemmi::gemmi_cpp` CMake 타겟 |
| **lodepng** | zlib | PNG 인코딩 (스크린샷) | 정적 라이브러리 |

### 빌드 방법
```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## 4. 핵심 데이터 구조

### Atom (`src/structure/Atom.hpp`)
```cpp
struct Atom {
    float x, y, z;                      // 3D 좌표
    char structure = 'x';               // 'x' (coil) | 'H' (helix) | 'S' (sheet)
    float bfactor = 0.0f;               // pLDDT 신뢰도
    bool is_interface = false;           // 계면 잔기 여부
    bool is_aligned = false;            // 정렬 영역 여부
    float conservation_score = -1.0f;   // 0.0-1.0 또는 -1.0 (미설정)
    int residue_number = -1;
    std::string residue_name;           // "GLU", "ALA" 등
};
```

### Protein (`src/structure/Protein.hpp/.cpp`)
- `init_atoms`: 원래 Å 좌표 (map<chainID, vector<Atom>>)
- `screen_atoms`: 변환된 화면 좌표
- `chain_res_count`: 체인별 잔기 수
- `bounding_box`: AABB
- 주요 메서드:
  - `load_data()`: PDB/mmCIF 로드 (gemmi 사용)
  - `normalize()`: 중심화 + 스케일링
  - `compute_interface(threshold)`: CA-CA 거리 기반 계면 계산
  - `compute_aligned_regions_from_aln()`: 정렬 영역 표시
  - `apply_ut_to_init_atoms(U, T)`: UT 변환 적용
  - `apply_conservation_scores()`: 보존도 점수 적용
  - `set_rotate/shift/scale()`: 변환 설정

### Screen (`src/visualization/Screen.hpp/.cpp`)
- 메인 렌더러 (최대 6개 단백질)
- `logicalPixels`: Braille 2W×4H 버퍼
- 깊이 안개: `depth_base_min_z`, `depth_base_max_z`로 캘리브레이션
- 주요 메서드:
  - `project()`: 3D→2D 투영 + Braille 렌더
  - `assign_colors_to_points()`: 모드별 색상 할당
  - `print_screen_braille()`: 터미널 출력
  - `draw_line()`: Bresenham 래스터화
  - `handle_input()`: 키보드/마우스 입력 처리
  - `load_next_hit()`: Foldseek hit 로드
  - `apply_foldmason_superposition()`: FoldMason 중첩

### RenderPoint (`src/visualization/RenderPoint.hpp`)
```cpp
struct RenderPoint {
    int x, y;
    float depth;
    int color_id;
    std::string chainID;
    char structure;             // 'H', 'S', 'x'
    bool is_interface;
    bool is_aligned;
    float bfactor;
    float conservation_score;
    int residue_number;
    char residue_name[4];
    int depth_band;             // 0=near(밝음), 1=mid(보통), 2=far(어두움)
};
```

### Panel (`src/visualization/Panel.hpp/.cpp`)
- 화면 우측 정보 패널
- 체인 정보, Foldseek hit 메타데이터, FoldMason 정보, 호버 잔기 상세

### Camera (`src/visualization/Camera.hpp/.cpp`)
- 3D→2D 투영 (원근 투영)
- PNG 스크린샷 기능 (lodepng)

### Palette (`src/visualization/Palette.hpp`)
- 250+ ncurses 색상 쌍 정의
- 모드별(protein/chain/rainbow/pLDDT/conservation/interface/aligned) 색상
- 깊이 안개용 near/mid/far 색상 변형

---

## 5. 데이터 흐름 (입력 → 렌더링)

```
입력 파일 (PDB/mmCIF)
    ↓
Parameters (CLI 인자 파싱)
    ↓
Protein::load_data() → init_atoms (원래 Å 좌표)
    ↓
[선택] UTMatrix 적용 → apply_ut_to_init_atoms(U, T)
    ↓
normalize() → screen_atoms (중심화 + 스케일링)
    ↓
[선택] 계면 계산 → is_interface 플래그
[선택] 정렬 계산 → is_aligned 플래그
[선택] 보존도 계산 → conservation_score
[선택] 2차 구조 예측/로드 → structure 문자
    ↓
Camera 투영 → 3D→2D 화면 좌표 + 깊이
    ↓
draw_line() → Bresenham 래스터화 (CA-CA 백본)
    ↓
깊이 밴드 분류: near(0.0-0.25) / mid(0.25-0.60) / far(0.60-1.0)
    ↓
색상 할당 → 모드별 + 깊이 밴드별 색상 쌍
    ↓
Braille 매핑 → 2×4 논리 픽셀 → 유니코드 Braille 문자
    ↓
ncurses 출력 + Panel 정보 표시
```

---

## 6. 렌더링 파이프라인

### Braille 서브픽셀 렌더링
- 터미널 셀 1개 = 2×4 논리 픽셀 그리드
- `logicalPixels[2*width × 4*height]` 버퍼
- 유니코드 Braille 문자로 변환: ⠿ (all on) ~ ⠀ (all off)
- 일반 문자 대비 8배 해상도

### 투영 공식
```cpp
const float FOV = 90.0f;
float fovRads = 1.0f / tan((FOV / zoom_level) * 0.5 / 180 * PI);
int screenX = (int)((x/z * fovRads + 1.0) * 0.5 * logical_w);
int screenY = (int)((1.0 - y/z * fovRads) * 0.5 * logical_h);
```

### 깊이 안개 (Depth Fog)
3단계 깊이 기반 색상으로 깊이감 제공:
```
band 0 (near, 0.0-0.25): 밝고 선명한 색상
band 1 (mid, 0.25-0.60): 기본 색상
band 2 (far, 0.60-1.0): 어둡고 색조 유지된 색상
```

---

## 7. 색상 모드 (7가지)

### protein (기본)
- 9가지 기본색 (olive, turquoise, navy, purple, pink, coral, brown, orange, red)
- 단백질당 1색 (순환)
- `-s` 플래그: helix=yellow(41), sheet=cyan(42), coil=dim protein color

### chain
- 15가지 체인 색상 (pairs 21-35)
- 체인별 고유 색상

### rainbow
- 20단계 그라데이션 (pairs 51-70)
- CA 인덱스 → 색조 순차 매핑

### plddt (AlphaFold 신뢰도)
| pLDDT 범위 | 색상 | pair |
|-----------|------|------|
| ≥90 | 파랑 | 71 |
| 70-90 | 시안 | 72 |
| 50-70 | 노랑 | 73 |
| <50 | 주황 | 74 |

### interface (계면)
- 계면 잔기 (CA-CA <8.0Å, 다른 체인): 마젠타 (pair 43)
- 비계면: 올리브 dim (pair 44)

### conservation (보존도)
- 10단계 그라데이션 (pairs 75-84): 파랑(가변) → 빨강(보존)
- Shannon 엔트로피 또는 MSA 보존도 기반

### aligned (정렬)
- 정렬된 잔기: 밝은 단백질 색상 (pairs 101-109)
- 비정렬: 회색 dim (pair 110)

---

## 8. 입력 포맷

| 포맷 | 파서 | 파일 타입 | 내용 |
|------|------|----------|------|
| PDB/mmCIF | gemmi | .pdb, .cif | 원자 좌표, 체인, B-factor, HELIX/SHEET |
| Foldseek .m8 | FoldseekParser | 12/17/21/29 컬럼 | query/target, evalue, qaln/taln, U(3×3), T(3) |
| FoldMason JSON | FoldMasonParser | .json | MSA 엔트리, Cα 좌표, 3Di 구조 |
| FoldMason FASTA | FoldMasonParser | .fa/.fasta | gap 포함 AA 서열 |
| MSA | MSAParser | .fasta, .a3m, .m8 | 다중 서열 (보존도 계산용) |
| 체인 선택 | Parameters | .tsv | 인덱스 + 체인 ID |
| UT 행렬 | Parameters | .tsv | 회전(3×3) + 이동(3) |

### Foldseek .m8 컬럼 포맷
- **12 컬럼**: 기본 (정렬 없음)
- **17 컬럼**: qaln/taln 포함
- **21 컬럼**: alis 포맷 (prob, evalue, lengths, taxinfo) ← easy-search 기본 출력
- **29 컬럼**: 전체 (U, T, qaln, taln, lddt, qtmscore)

---

## 9. CLI 옵션

```bash
./StrucTTY <input_files...> [OPTIONS]

입력:
  최대 6개 PDB/mmCIF 파일

옵션:
  -m, --mode <MODE>        색상 모드: protein|chain|rainbow|plddt|interface|conservation|aligned
  -c, --chains <FILE>      TSV: 인덱스 + 체인 ID
  -s, --structure          2차 구조 표시 (나선/시트)
  -p, --predict            2차 구조 예측 (파일에 없을 때)
  -ut, --utmatrix <FILE>   TSV: 회전(3×3) + 이동(3) 정렬용
  --msa <FILE>             MSA 파일 (보존도 모드용)
  -fs, --foldseek <FILE>   Foldseek .m8 결과 (hit 탐색)
  --db-path <DIR>          PDB 디렉토리 (Foldseek hit 로드용)
  -fm, --foldmason <FILE>  FoldMason JSON/FASTA MSA
  -b, --benchmark          벤치마크 모드 (2000 프레임)
  -n, --nopanel            정보 패널 숨기기
  --help                   도움말
```

---

## 10. 키보드 조작

| 키 | 동작 |
|----|------|
| 0-9 | 제어할 단백질 선택 (0=전체, 1-9=개별) |
| W/w | 위로 이동 (+Y) |
| A/a | 왼쪽으로 이동 (-X) |
| S/s | 아래로 이동 (-Y) |
| D/d | 오른쪽으로 이동 (+X) |
| X/x | X축 회전 |
| Y/y | Y축 회전 |
| Z/z | Z축 회전 |
| R/r | 줌 인 |
| F/f | 줌 아웃 |
| Q/q | 종료 |
| 마우스 | 호버 시 잔기 정보 패널 표시 |

---

## 11. 고급 기능

### Foldseek 통합 (`FoldseekParser`)
- `.m8` 파일 파싱 (다중 컬럼 포맷 지원)
- qaln/taln 정렬 문자열로 정확한 영역 매칭
- U/T 변환 행렬로 구조 중첩
- hit 인터랙티브 탐색 (load_next_hit)

### FoldMason MSA 중첩 (`FoldMasonParser`)
- JSON: Cα 좌표 + 3Di 구조 포함 엔트리
- FASTA: 정렬된 서열만
- 중첩: 정렬된 잔기 쌍 → Kabsch SVD
- 보존도: 컬럼별 Shannon 엔트로피 (gap 인지)

### 2차 구조 시각화 (`StructureMaker`)
- 나선(Helix): 원통 기하 (16단계 원, ±2.5Å 반경)
- 시트(Sheet): 리본 기하 (±4단계 오프셋, 0.28Å/단계)
- 예측: `SSPredictor` (거리/비틀림각 투표)

### PDB 다운로드 (`PDBDownloader`)
- 10종 DB 자동 감지: PDB, AlphaFold, ESMAtlas, CATH, BFVD, GMGCL, TED, BFMD 등
- 캐싱: `~/.cache/structty/pdb/`
- curl/wget 기반 다운로드

### 스크린샷 (`Camera`)
- lodepng를 통한 PNG 내보내기
- Braille 렌더 → RGBA 이미지 변환

### 벤치마크 (`Benchmark`)
- CSV 로깅: FPS, 지연시간, 프레임당 시간
- 워밍업 + 측정 실행 (기본 2000 프레임)

---

## 12. 메인 실행 흐름 (`structty.cpp`)

```
main()
  1. Parameters 파싱 (CLI 인자)
  2. ncurses 초기화
  3. Screen 객체 생성
  4. 단백질 파일 로드 (set_protein)
  5. [선택] UTmatrix 로드 및 적용
  6. 단백질 정규화 (중심, 스케일)
  7. [선택] 계면 계산 (-m interface)
  8. [선택] 정렬 영역 계산 (-m aligned, -fs, -fm)
  9. [선택] Foldseek hit 로드 → 변환 적용
  10. [선택] MSA 보존도 로드 (-m conservation + --msa)
  11. [선택] FoldMason 로드 → 중첩 + 보존도
  12. 벤치마크 모드? → 2200 프레임 (200 워밍업 + 2000 측정)
  13. 인터랙티브 루프:
      - draw_screen()
      - handle_input() → 회전/이동/줌
      - 'q' 누를 때까지 반복
  14. 정리 (endwin())
```

---

## 13. 수정 대상 파일 요약 (전체)

| 파일 | 역할 | 크기 |
|------|------|------|
| `src/structty.cpp` | 메인 진입점, 초기화, 이벤트 루프 | ~중형 |
| `src/structure/Atom.hpp` | 원자 데이터 구조 | ~소형 |
| `src/structure/Protein.hpp/.cpp` | 단백질 로드/변환/계산 | ~대형 |
| `src/structure/Parameters.hpp/.cpp` | CLI 파싱 | ~중형 |
| `src/structure/StructureMaker.hpp/.cpp` | 2차 구조 기하 | ~중형 |
| `src/structure/SSPredictor.hpp/.cpp` | 2차 구조 예측 | ~중형 |
| `src/structure/FoldseekParser.hpp/.cpp` | Foldseek 파서 | ~중형 |
| `src/structure/FoldMasonParser.hpp/.cpp` | FoldMason 파서 | ~대형 |
| `src/structure/MSAParser.hpp/.cpp` | MSA 파서 | ~중형 |
| `src/structure/PDBDownloader.hpp/.cpp` | PDB 다운로드 | ~중형 |
| `src/visualization/Screen.hpp/.cpp` | 메인 렌더러 | ~대형 |
| `src/visualization/Panel.hpp/.cpp` | 정보 패널 | ~중형 |
| `src/visualization/Camera.hpp/.cpp` | 카메라/스크린샷 | ~중형 |
| `src/visualization/Palette.hpp` | 색상 팔레트 정의 | ~중형 |
| `src/visualization/RenderPoint.hpp` | 렌더 포인트 구조 | ~소형 |
| `src/utils/Common.hpp` | 타임스탬프 유틸 | ~소형 |
| `src/utils/Curses.hpp` | ncurses 감지 | ~소형 |
| `src/utils/Benchmark.hpp` | 성능 측정 | ~소형 |

---

## 14. 예제 데이터 (`example/`)

| 파일/디렉토리 | 타입 | 용도 |
|-------------|------|------|
| `*.cif` | mmCIF | 샘플 단백질 구조 (1NPL, 3HGM, 9N47, 9FL9, 3A0C, AF-*, 8KGM, 1CJK, 1M1C) |
| `chainfile_9N47.tsv` | TSV | 체인 선택: 0 A,B,C / 1 K,G / 2 M |
| `utfile_1npl_3a0c.tsv` | TSV | UT 행렬 (회전 3×3 + 이동 3) |
| `foldseek_result/` | Foldseek | 다중 DB hit 결과 (BFVD, afdb50, afdb-swissprot 등) |
| `foldmason_result/` | FoldMason | JSON MSA + FASTA 서열 |
| `msa_result/` | MSA | A3M/FASTA/M8 서열 정렬 |

---

## 15. 키 7/8/9 개별 단백질 선택 실패 분석

### 버그 위치
`src/visualization/Screen.cpp:982`

### 문제 코드
```cpp
// Screen.cpp:967-985
switch(key){
    case 48: case 49: case 50: case 51: case 52:
    case 53: case 54: case 55: case 56: case 57:  // '0'~'9'
        if (key == 48){
            structNum = -1;       // '0' → 전체 선택
        }
        else if (key - 48 <= data.size()) {  // ← 여기가 문제
            structNum = key - 49; // '1'~'9' → 개별 선택 (인덱스 0~8)
        }
        break;
}
```

### 근본 원인

조건식 `key - 48 <= data.size()` 자체는 사실 경계 검사가 맞다 (`key - 49 < data.size()`와 동치). 문제는 **조건이 실패했을 때 `else` 절이 없다**는 것.

data.size()가 6일 때:

| 키 | key-48 | 조건 (6개 로드) | 결과 |
|----|--------|----------------|------|
| 1~6 | 1~6 | 1~6 ≤ 6? YES | structNum = 0~5 ✓ |
| 7 | 7 | 7 ≤ 6? **NO** | structNum **변경 안 됨** ✗ |
| 8 | 8 | 8 ≤ 6? **NO** | structNum **변경 안 됨** ✗ |
| 9 | 9 | 9 ≤ 6? **NO** | structNum **변경 안 됨** ✗ |

조건이 실패하면 `structNum`이 이전 값 그대로 유지된다. 기본값이 `-1`(전체 선택)이므로, 키 7/8/9를 눌러도 '0'을 누른 것과 동일하게 전체가 움직인다.

### 결론

**7/8/9번째 단백질이 `data`에 실제로 로드되어 있어야 한다.** `data.size() >= 7/8/9`일 때만 해당 키가 작동한다. 현재 로드된 단백질 수가 7개 미만이면 조건이 실패하는 것은 정상 동작이다.

### 확인 필요

1. **현재 몇 개의 단백질이 로드되어 있는 상태인가?** 7개 이상 로드된 상태에서도 안 되는 것인지, 아니면 6개 이하에서 7/8/9를 시도하는 것인지
2. 만약 7개 이상이 실제로 로드되어 있다면, `data.size()` 값이 기대와 다를 수 있음 (이 경우 별도 조사 필요)

### 수정 방향 (시나리오별)

**시나리오 A: 실제로 7개 이상 로드된 상태에서 안 되는 경우**
→ `data.size()` 값 확인 필요. 디버그 출력으로 확인 가능

**시나리오 B: 범위 밖 키를 눌렀을 때 "무시" 대신 피드백이 필요한 경우**
→ `else` 절 추가하여 범위 밖이면 `structNum = -1`로 전체 선택 fallback:
```cpp
else if (key - 48 <= (int)data.size()) {
    structNum = key - 49;
} else {
    // 범위 밖이면 전체 선택으로 fallback (또는 무시)
    structNum = -1;
}
```
