# libstructty-render 추출 구현 계획

## 배경

Martin이 Foldseek 연동을 위한 두 가지 접근 방식을 거부했다:
1. ncurses를 Foldseek에 링크 — 너무 비대
2. StrucTTY를 subprocess로 호출 — 수용 불가

Foldseek 통합 계획은 확정되어 있으며, 렌더링 엔진을 ncurses 독립 standalone 라이브러리로 추출해야 한다.

---

## 목표 및 범위

### 목표
`libstructty-render`라는 ncurses 의존성 없는 정적 라이브러리를 만들어, Foldseek이 단일 단백질 구조를 ANSI escape code로 터미널에 표시할 수 있게 한다.

### 포함 범위
- 3D → 2D 원근 투영 (Camera 로직)
- Z-버퍼 / 깊이 정렬
- 색상 모드별 컬러 할당 (protein/chain/rainbow/plddt/interface/conservation/aligned)
- 3-밴드 깊이 안개 (near/mid/far)
- 점자(Braille) 문자 생성
- ANSI escape code stdout 출력 백엔드
- CMake `INTERFACE` 헤더 + 정적 라이브러리 타겟

### 제외 범위 (StrucTTY 내부에 유지)
- 패널 레이아웃 및 렌더링 (`Panel.cpp`)
- 키보드/마우스 입력 처리 (`getch`, `handle_input`)
- 멀티 구조 비교 UI
- ncurses TUI 라이프사이클 (`initscr`, `endwin`)
- Foldseek 히트 탐색, FoldMason 중첩 UI

### Foldseek 측 기대 동작
```
libstructty-render 링크
→ render_structure_ansi(atoms, width, height, mode)
→ ANSI 이스케이프 문자열 반환 또는 stdout 출력
→ 인터랙티브 없음, 패널 없음
```

---

## 영향받는 파일 목록

### 수정 대상
| 파일 | 변경 내용 |
|------|---------|
| `src/visualization/RenderPoint.hpp` | `color_id` 의미 재정의 (ncurses pair# → 논리 팔레트 인덱스) |
| `src/visualization/Palette.hpp` | ANSI 변환 테이블/함수 추가 (`palette_to_ansi_fg()`) |
| `src/visualization/Screen.cpp/hpp` | 렌더링 코어 추출 후 Renderer 위임으로 리팩터 |
| `src/visualization/Camera.cpp/hpp` | 투영 로직을 ncurses 무관하게 분리 |
| `CMakeLists.txt` (루트) | `libstructty-render` 타겟 추가, StrucTTY는 기존 유지 |
| `lib/CMakeLists.txt` | lodepng는 StrucTTY 전용으로 분리 (render lib에 불포함) |

### 신규 생성 파일
| 파일 | 역할 |
|------|------|
| `src/render/Renderer.hpp` | 렌더링 코어 클래스 헤더 (ncurses 미포함) |
| `src/render/Renderer.cpp` | 투영 + Z-버퍼 + 컬러 할당 + 점자 생성 구현 |
| `src/render/AnsiOutput.hpp` | ANSI 출력 백엔드 헤더 |
| `src/render/AnsiOutput.cpp` | escape code 생성 및 stdout 출력 구현 |
| `include/structty_render.h` | Foldseek이 include할 공개 C/C++ API 헤더 |
| `src/render/CMakeLists.txt` | render 라이브러리 빌드 설정 |

---

## 구현 단계

### Phase 0: 설계 확정 및 감사 (코드 수정 없음)

- [x] **0-1. ncurses 호출 전수 감사**
  - `Screen.cpp`에서 `attron`, `attroff`, `mvaddstr`, `waddstr`, `COLOR_PAIR`, `init_pair`, `mvprintw` 등 ncurses 심볼 전수 목록화
  - 렌더링 파이프라인 경로(투영→Z버퍼→컬러→점자 출력)에 속하는 것과 TUI 관리에 속하는 것을 분류
  - 결과를 이 파일 Appendix에 기록 → **Appendix A 참조**

- [x] **0-2. `color_id` 의존 경로 파악**
  - `RenderPoint.color_id`가 ncurses pair 번호로 사용되는 모든 지점 목록화
  - `Palette.hpp`의 컬러 쌍 번호와 xterm-256 색상 번호 매핑 관계 문서화
  - ANSI 변환에 사용할 수 있는 xterm-256 ID2RGBA 테이블이 이미 `Palette.hpp`에 있음을 확인
  - **핵심 발견: `Camera::renderPoint2image()`에 이미 완전한 color_id→xterm-256 변환 테이블이 구현되어 있음 → Appendix B 참조**

- [x] **0-3. 공개 API 표면 확정**
  - Foldseek 측과 합의된 호출 시그니처를 `include/structty_render.h`에 초안 작성
  - 단일 구조 / 다중 구조 여부 결정
  - 색상 모드 선택권 범위 결정 (전체 7종 vs Foldseek 필요 서브셋)
  - **결정: 단일 구조, 전체 7종 색상 모드 지원, 시그니처는 Appendix C 참조**

---

### Phase 1: 컬러 시스템 이식

이 Phase의 핵심은 기존 StrucTTY의 ncurses 컬러 쌍 시스템을 그대로 유지하면서, ANSI 출력에 필요한 변환 경로를 추가하는 것이다. `RenderPoint.color_id`의 값 자체는 바꾸지 않는다.

- [x] **1-1. Palette.hpp에 ANSI 변환 함수 추가**
  - `palette_to_xterm256_fg(int color_id) -> int` 구현
    - color_id(ncurses pair 번호) → xterm-256 전경색 번호 반환
    - 기존 `Palette.hpp`의 xterm-256 ID2RGBA 테이블 활용
  - `palette_to_ansi_fg_str(int color_id) -> std::string` 구현
    - `"\033[38;5;{n}m"` 포맷 반환 (xterm-256 전경색 escape)
  - `palette_to_ansi_reset() -> const char*` 구현
    - `"\033[0m"` 반환

- [x] **1-2. 깊이 밴드 → ANSI 밝기 매핑 확정**
  - `depth_band` (0=near, 1=mid, 2=far)에 대한 ANSI 밝기 조절 방식 결정
    - near: 기본 색상 그대로
    - mid: 기본 색상 그대로
    - far: xterm-256 회색 계열 혼합 또는 dim escape (`\033[2m`) 사용 결정
  - `RenderPoint.depth_band` 필드가 이미 있으므로 추가 구조 변경 불필요

- [x] **1-3. 변환 함수 단위 검증**
  - 9개 단백질 색상 (pair 1-9) ANSI 변환 출력 육안 확인
  - 7개 색상 모드 각각에 대해 xterm-256 번호 매핑이 올바른지 확인
  - 깊이 near/mid/far 세 밴드에서 색상 변화가 의도대로인지 확인

---

### Phase 2: 렌더링 코어 추출 (Screen.cpp 분리)

Screen.cpp 1942 lines에서 렌더링 파이프라인에 해당하는 로직을 `Renderer` 클래스로 이관한다. Screen은 Renderer를 내부적으로 사용하는 구조로 변경한다.

- [x] **2-1. Renderer 클래스 헤더 설계 (`src/render/Renderer.hpp`)**
  ```cpp
  class Renderer {
  public:
      Renderer(int width, int height, const std::string& mode);
      
      // 원자 목록을 받아 logicalPixels 버퍼 생성
      void render(const std::vector<RenderAtom>& atoms, bool show_structure);
      
      // 채워진 버퍼 접근
      const std::vector<RenderPoint>& get_pixels() const;
      int get_logical_width() const;   // 2 * width
      int get_logical_height() const;  // 4 * height
      
  private:
      int width_, height_;
      std::string mode_;
      std::vector<RenderPoint> logical_pixels_;
      
      void project_atoms(...);
      void apply_zbuffer(...);
      void assign_colors(...);
      void generate_braille(...);
  };
  ```

- [x] **2-2. RenderAtom 구조체 정의 (`src/render/Renderer.hpp`)**
  - Protein/Atom과 무관한 순수 데이터 구조체 (ncurses, gemmi 미포함)
  ```cpp
  struct RenderAtom {
      float x, y, z;
      char structure;           // 'x', 'H', 'S'
      float bfactor;
      bool is_interface;
      bool is_aligned;
      float conservation_score; // -1 = unset
      int residue_number;
      char residue_name[4];
      char chain_id;
      int protein_index;        // 0-8, 어느 단백질 소속인지
  };
  ```

- [x] **2-3. `project_atoms()` 추출**
  - Screen.cpp의 3D→2D 원근 투영 로직 이관
  - Camera.cpp의 `focal_offset`, FOV 상수 의존성 파악 후 파라미터화
  - ncurses 호출 없음 검증

- [x] **2-4. `apply_zbuffer()` 추출**
  - 논리 픽셀 격자 Z-버퍼 및 깊이 정렬 로직 이관
  - `logicalPixels` 버퍼 초기화 포함

- [x] **2-5. `assign_colors()` 추출**
  - 7개 색상 모드별 컬러 할당 로직 이관
  - `color_id`에 논리 팔레트 인덱스 저장 (기존 값 그대로)
  - `depth_band` 설정 로직 포함

- [x] **2-6. `generate_braille()` 추출**
  - 논리 픽셀 2×4 → 점자 유니코드 문자 생성 로직 이관
  - `RenderPoint.pixel` 채우기 포함

- [x] **2-7. Screen.cpp 리팩터**
  - `draw_screen()`이 내부적으로 `Renderer::render()` 호출 후 ncurses로 출력하는 구조로 변경
  - Screen의 기존 기능(패널, 입력, Foldseek 탐색 등) 변경 없음
  - `Protein::screen_atoms` → `RenderAtom` 변환 어댑터 함수 추가 (`Screen::to_render_atoms()`)

---

### Phase 3: ANSI 출력 백엔드 구현

- [ ] **3-1. AnsiOutput 클래스 구현 (`src/render/AnsiOutput.cpp/hpp`)**
  ```cpp
  class AnsiOutput {
  public:
      // Renderer의 픽셀 버퍼를 ANSI escape string으로 변환
      static std::string to_ansi_string(
          const std::vector<RenderPoint>& pixels,
          int logical_width,
          int logical_height
      );
      
      // stdout으로 직접 출력 (Foldseek 사용 케이스)
      static void print_to_stdout(
          const std::vector<RenderPoint>& pixels,
          int logical_width,
          int logical_height
      );
  };
  ```

- [ ] **3-2. 점자 문자 + ANSI 컬러 결합 로직 구현**
  - 각 터미널 문자 위치에서 `RenderPoint` 2×4 블록 수집
  - 전경색: `palette_to_ansi_fg_str(color_id)` 적용
  - 깊이 밴드 far의 경우 dim 처리 (`\033[2m`) 또는 색상 변형 적용
  - 점자 유니코드 출력 후 reset (`\033[0m`)

- [ ] **3-3. 줄 바꿈 및 커서 처리**
  - 각 논리 픽셀 행 → 터미널 행 매핑 (4 논리 픽셀 행 = 1 터미널 행)
  - 출력 후 커서 위치 초기화 여부 결정 (Foldseek UX 요구사항에 따라)

---

### Phase 4: 공개 API 헤더 및 CMake 타겟

- [ ] **4-1. 공개 API 헤더 작성 (`include/structty_render.h`)**
  ```cpp
  // C++ API (내부용 + Foldseek)
  #pragma once
  #include "Renderer.hpp"
  #include "AnsiOutput.hpp"
  
  namespace structty {
  
  // 단일 호출 편의 함수
  // atoms: RenderAtom 벡터
  // width/height: 터미널 문자 단위 크기
  // mode: "protein" | "chain" | "rainbow" | "plddt" | "conservation" | "aligned"
  std::string render_to_ansi(
      const std::vector<RenderAtom>& atoms,
      int width,
      int height,
      const std::string& mode = "protein",
      bool show_structure = false
  );
  
  void render_to_stdout(
      const std::vector<RenderAtom>& atoms,
      int width,
      int height,
      const std::string& mode = "protein",
      bool show_structure = false
  );
  
  } // namespace structty
  ```

- [ ] **4-2. `src/render/CMakeLists.txt` 작성**
  ```cmake
  add_library(structty_render STATIC
      Renderer.cpp
      AnsiOutput.cpp
  )
  
  target_include_directories(structty_render
      PUBLIC
          ${PROJECT_SOURCE_DIR}/include
          ${CMAKE_CURRENT_SOURCE_DIR}
          ${PROJECT_SOURCE_DIR}/src/visualization  # RenderPoint.hpp, Palette.hpp
  )
  
  # ncurses 링크 없음 — 명시적 검증
  target_link_libraries(structty_render
      PRIVATE
          gemmi::gemmi_cpp  # RenderAtom이 Protein과 완전 분리되면 제거 가능
  )
  ```

- [ ] **4-3. 루트 CMakeLists.txt 수정**
  - `add_subdirectory(src/render)` 추가
  - `StrucTTY` 타겟에 `structty_render` 링크 추가
  - ncurses는 `StrucTTY`만 링크하고 `structty_render`는 미링크임을 주석으로 명시

- [ ] **4-4. ncurses 미링크 빌드 검증**
  - `structty_render` 타겟을 `ldd` 또는 `nm`으로 검사하여 ncurses 심볼 없음 확인
  - Linux: `nm -u libstructty_render.a | grep -i ncurses` 결과 비어야 함

---

### Phase 5: Foldseek 연동 검증

- [ ] **5-1. 최소 연동 테스트 바이너리 작성**
  - `example/render_test.cpp`: `structty_render.h`만 include하여 `.cif` 파일 하나를 로드하고 ANSI 출력하는 독립 프로그램 작성
  - CMake 옵션(`BUILD_RENDER_TEST=ON`)으로 선택적 빌드

- [ ] **5-2. 예제 파일로 출력 검증**
  - `example/1CJK-assembly1.cif` 로드 → `render_to_stdout()` 호출 → 터미널 출력 육안 확인
  - 7가지 색상 모드 각각 출력 확인

- [ ] **5-3. Foldseek 빌드 시스템 연동 확인**
  - Foldseek이 요구하는 CMake `find_package` 또는 `add_subdirectory` 방식 확인 후 지원 방식 결정
  - 필요 시 `structty_renderConfig.cmake` 패키지 파일 작성

---

## 인터페이스 변경사항

### RenderPoint.hpp — 변경 없음
`color_id`의 값(ncurses pair 번호)은 그대로 유지한다. 의미만 "논리 팔레트 인덱스"로 재해석한다. ANSI 변환은 출력 레이어에서 담당.

### Palette.hpp — 추가
```cpp
// 기존 코드 유지, 아래 함수 추가:

// color_id(ncurses pair 번호) → xterm-256 전경색 번호
int palette_to_xterm256_fg(int color_id);

// color_id → ANSI 전경색 escape string ("\033[38;5;Nm")
std::string palette_to_ansi_fg_str(int color_id);
```

### Screen.hpp — 추가
```cpp
// 기존 멤버 변경 없음, 아래 추가:
private:
    std::vector<RenderAtom> to_render_atoms() const;
    // screen_atoms를 RenderAtom 벡터로 변환 (어댑터)
```

### Screen.cpp — draw_screen() 내부 변경
```cpp
// Before: Screen 내부에서 직접 투영 + Z버퍼 + 컬러 + 점자 + ncurses 출력

// After:
void Screen::draw_screen(bool no_panel) {
    auto atoms = to_render_atoms();           // 어댑터
    renderer_.render(atoms, show_structure_); // 코어 위임
    
    // ncurses 출력 (StrucTTY 전용, 라이브러리 미포함)
    const auto& pixels = renderer_.get_pixels();
    // ... mvaddwstr 등 ncurses 호출 유지 ...
    
    if (!no_panel) draw_panel();
}
```

### 신규: RenderAtom (src/render/Renderer.hpp)
위 Phase 2-2 참조. Protein/Atom과 완전 분리된 순수 데이터 구조체.

### 신규: Renderer (src/render/Renderer.hpp)
위 Phase 2-1 참조.

### 신규: AnsiOutput (src/render/AnsiOutput.hpp)
위 Phase 3-1 참조.

### 신규: 공개 API (include/structty_render.h)
위 Phase 4-1 참조.

---

## 테스트 계획

### T1. 컬러 변환 정확성
- 각 색상 모드별로 `palette_to_xterm256_fg()` 반환값이 원래 ncurses 쌍의 전경색 xterm-256 번호와 일치하는지 확인
- 검증 대상: protein 9색, chain 15색, rainbow 20색조, plddt 4밴드, conservation 10단계, interface, aligned

### T2. 렌더링 동등성 (StrucTTY ncurses 출력 vs ANSI 출력)
- 동일 구조 파일(예: `1CJK-assembly1.cif`)에 대해 StrucTTY의 ncurses 출력과 새 ANSI 출력이 동일한 점자 문자 패턴을 가지는지 확인
- `RenderPoint.pixel` 값의 분포가 동일해야 함

### T3. 깊이 밴드 시각 확인
- near/mid/far 밴드가 ANSI 출력에서 구분되는지 육안 확인
- 앞면 구조가 뒷면보다 밝게 보여야 함

### T4. 7가지 색상 모드 ANSI 출력
- protein, chain, rainbow, plddt, interface, conservation, aligned 각각에 대해 `render_to_stdout()` 호출 후 출력 확인
- 예제 파일: `example/1CJK-assembly1.cif` (단일), `example/1NPL-assembly1.cif` + `example/3A0C-assembly1.cif` (비교)

### T5. ncurses 미의존 빌드 검증
- `structty_render` 단독 빌드 (StrucTTY 없이)가 성공하는지 확인
- `ldd` / `nm`으로 ncurses 심볼 미포함 확인

### T6. StrucTTY 기존 기능 회귀 없음
- Phase 2 리팩터 후 기존 StrucTTY 빌드 및 실행이 정상인지 확인
- 특히 ncurses 컬러 출력, 마우스 호버, 패널 표시, Foldseek 히트 탐색이 이전과 동일하게 동작하는지 확인

### T7. 최소 연동 테스트 바이너리 (Phase 5-1)
- `example/render_test.cpp`가 ncurses 링크 없이 빌드되고, 터미널에 구조를 출력하는지 확인

---

## 예상 리스크 및 주의사항

### R1. color_id → ANSI 변환 테이블 불완전 (높음)
`Palette.hpp`의 ncurses pair 정의와 xterm-256 번호 사이의 매핑을 수동으로 추출해야 한다. pair 번호 250개 이상에 대한 전경색 xterm-256 번호 목록이 Palette.hpp에 명시적으로 정리되어 있지 않으면, `init_pair()` 호출을 역추적해야 한다.

**대응:** Phase 0-2에서 이 매핑을 먼저 완전히 문서화한 후 Phase 1로 진행.

### R2. Screen.cpp 분리 작업 범위 증가 가능 (높음)
1942 lines의 God class에서 렌더링 코어만 분리하는 과정에서, 예상보다 많은 상태 변수가 얽혀 있을 수 있다. 특히 `screen_width`, `screen_height`, `proteins[]`, `camera_` 등이 투영 로직에 직접 참조된다.

**대응:** Phase 2-3 착수 전, 투영 함수에서 사용되는 모든 멤버 변수 목록을 작성하고 파라미터로 전달 가능한지 판단.

### R3. StructureMaker가 렌더링 파이프라인에 포함될 경우 (중간)
2차 구조 시각화(나선 실린더, 시트 리본) 기하학 생성 로직(`StructureMaker.cpp`)이 렌더링 코어에 포함되어야 하는지, 아니면 호출부(Screen)에서 미리 원자로 변환해서 넘기는지 결정이 필요하다.

**대응:** `show_structure=false`가 Foldseek의 초기 사용 케이스라면 우선 제외. `RenderAtom`에 기하학 점들도 포함 가능하도록 확장성은 열어두기.

### R4. Foldseek 빌드 시스템 요구사항 불확실 (중간)
Foldseek이 `add_subdirectory` 방식을 요구하는지, 설치된 패키지(`find_package`)를 요구하는지에 따라 CMake 설정이 달라진다.

**대응:** Phase 5-3에서 Foldseek 측과 빌드 방식 확인 후 진행.

### R5. 점자 문자와 ANSI escape 결합 시 터미널 호환성 (낮음)
일부 터미널에서 xterm-256 전경색과 Unicode Braille 문자를 함께 쓸 때 렌더링 이슈가 발생할 수 있다. StrucTTY가 현재 ncurses로 해결하고 있는 wide character 처리를 ANSI 경로에서는 직접 다뤄야 한다.

**대응:** Phase 3-3에서 UTF-8 로케일 설정(`setlocale(LC_ALL, "")`) 호출 포함.

---

## 작업 순서 요약

```
Phase 0 (설계/감사) ✓ 완료
  ├─ 0-1: ncurses 호출 전수 감사 ✓
  ├─ 0-2: color_id 의존 경로 파악 + 변환 테이블 초안 ✓
  └─ 0-3: 공개 API 표면 확정 ✓

Phase 1 (컬러 이식)
  ├─ 1-1: Palette.hpp ANSI 변환 함수 추가
  ├─ 1-2: 깊이 밴드 ANSI 처리 방식 확정
  └─ 1-3: 변환 함수 단위 검증

Phase 2 (코어 추출 — 가장 큰 공수)
  ├─ 2-1: Renderer 헤더 설계
  ├─ 2-2: RenderAtom 구조체 정의
  ├─ 2-3: project_atoms() 추출
  ├─ 2-4: apply_zbuffer() 추출
  ├─ 2-5: assign_colors() 추출
  ├─ 2-6: generate_braille() 추출
  └─ 2-7: Screen.cpp 리팩터 + 회귀 확인 (T6)

Phase 3 (ANSI 출력)
  ├─ 3-1: AnsiOutput 클래스 구현
  ├─ 3-2: 점자 + 컬러 결합 로직
  └─ 3-3: 줄 바꿈 및 커서 처리

Phase 4 (CMake 타겟)
  ├─ 4-1: 공개 API 헤더
  ├─ 4-2: src/render/CMakeLists.txt
  ├─ 4-3: 루트 CMakeLists.txt 수정
  └─ 4-4: ncurses 미링크 검증 (T5)

Phase 5 (연동 검증)
  ├─ 5-1: 최소 연동 테스트 바이너리
  ├─ 5-2: 예제 파일 출력 검증 (T1~T4, T7)
  └─ 5-3: Foldseek 빌드 연동 확인
```

---

## Appendix A: ncurses 호출 전수 감사 결과

### 분류 기준
- **C: 렌더링 경계** — 라이브러리에서 ANSI 경로로 교체 대상
- **TUI**: StrucTTY에 유지, 라이브러리 불포함

### Screen.cpp ncurses 호출 목록

| 라인 | 심볼 | 위치(함수) | 분류 | 비고 |
|------|------|----------|------|------|
| 18 | `start_color()` | Screen constructor | TUI | ncurses 컬러 서브시스템 초기화 |
| 19 | `use_default_colors()` | Screen constructor | TUI | 터미널 기본 배경색(-1) 허용 |
| 20,148 | `init_pair()` ×250+ | `init_color_pairs()`, `set_protein()` | TUI | 250개 이상 컬러 쌍 사전 할당 |
| 28 | `keypad(stdscr, TRUE)` | Screen constructor | TUI | 함수키/화살표 키 활성화 |
| 29 | `mousemask()` | Screen constructor | TUI | 마우스 이벤트 활성화 |
| 30 | `mouseinterval()` | Screen constructor | TUI | 클릭 인터벌 0으로 설정 |
| 751 | `getmaxyx(stdscr, rows, cols)` | `draw_screen()` | TUI | 패널 높이 계산용 터미널 크기 조회 |
| 778 | `erase()` | `draw_screen()` | TUI | 화면 지우기 |
| 795 | `refresh()` | `draw_screen()` | TUI | 버퍼 화면 출력 |
| 826 | `getmaxyx(stdscr, rows, cols)` | `print_screen_braille()` | **C** | 렌더링 범위 제한 — 라이브러리에서 `width`/`height` 파라미터로 대체 |
| 867 | `attron(COLOR_PAIR(best_color_id))` | `print_screen_braille()` | **C** | 전경색 설정 — ANSI `\033[38;5;{n}m`으로 대체 |
| 868 | `mvaddstr(row, tx, utf8)` | `print_screen_braille()` | **C** | 커서이동+문자출력 — 위치 기반 string 빌딩으로 대체 |
| 869 | `attroff(COLOR_PAIR(best_color_id))` | `print_screen_braille()` | **C** | 색상 해제 — ANSI `\033[0m`으로 대체 |
| 935 | `refresh()` | `update_hover_info()` | TUI | 호버 패널 부분 갱신 후 출력 |
| 941 | `getch()` | `handle_input()` | TUI | 키 입력 대기 |
| 957 | `getmouse(&event)` | `handle_input_impl()` | TUI | 마우스 이벤트 읽기 |
| 976 | `flushinp()` | `handle_input_impl()` | TUI | 입력 버퍼 비우기 |

### Panel.cpp ncurses 호출 목록 (모두 TUI 분류)

| 심볼 | 용도 |
|------|------|
| `move(r, c)` | 커서 이동 |
| `clrtoeol()` | 현재 줄 끝까지 지우기 |
| `addnstr(s, k)` | 문자열 출력 |
| `addch('-')` | 구분선 문자 출력 |
| `attron(COLOR_PAIR(...))` | 파일명/체인 강조 색상 시작 |
| `attroff(COLOR_PAIR(...))` | 색상 종료 |

### 결론
**라이브러리에서 교체해야 하는 ncurses 호출은 `print_screen_braille()` 내 4개뿐이다.**  
나머지는 모두 TUI 레이어에 속하며 StrucTTY에 그대로 유지된다.

---

## Appendix B: color_id → xterm-256 변환 테이블

### 핵심 발견
`Camera::renderPoint2image()` (`src/visualization/Camera.cpp`:41-71)에 이미 color_id → xterm-256 변환의 완전한 구현이 존재한다. Phase 1의 `palette_to_xterm256_fg()` 함수는 이 로직을 Palette.hpp로 추출하기만 하면 된다.

### 전체 매핑 표

| color_id 범위 | xterm-256 소스 | 설명 |
|--------------|---------------|------|
| 1–9 | `PROTEIN_COLORS[cid-1]` | 단백질 vivid 9색 (mid) |
| 11–19 | `PROTEIN_DIM_COLORS[cid-11]` | 단백질 dim 9색 (coil, mid) |
| 21–35 | `CHAIN_COLORS[cid-21]` | 체인 15색 (mid) |
| 41 | `226` | 나선 노란색 |
| 42 | `51` | 시트 청록색 |
| 43 | `INTERFACE_COLOR` (201) | 인터페이스 마젠타 (mid) |
| 44 | `INTERFACE_DIM_COLOR` (58) | 비인터페이스 dim (mid) |
| 45 | `ALIGNED_COLOR` (46) | 정렬 밝은 초록 (mid) |
| 46 | `ALIGNED_DIM_COLOR` (58) | 비정렬 dim (mid) |
| 51–70 | `RAINBOW[cid-51]` | 무지개 20색조 (mid) |
| 71–74 | `PLDDT_COLORS[cid-71]` | pLDDT 4밴드 (mid) |
| 75–84 | `CONSERVATION_COLORS[cid-75]` | 보존도 10단계 (mid) |
| 101–109 | `PROTEIN_BRIGHT_COLORS[cid-101]` | 정렬 bright 9색 (mid) |
| 110 | `ALIGNED_NONALIGNED_DIM` (240) | 비정렬 회색 (mid) |
| 120–128 | `PROTEIN_NEAR_COLORS[cid-120]` | 단백질 near (depth_band=0) |
| 130–144 | `CHAIN_NEAR_COLORS[cid-130]` | 체인 near |
| 145–159 | `CHAIN_FAR_COLORS[cid-145]` | 체인 far (depth_band=2) |
| 160–179 | `RAINBOW_NEAR[cid-160]` | 무지개 near |
| 180–199 | `RAINBOW_FAR[cid-180]` | 무지개 far |
| 200–208 | `PROTEIN_FAR_COLORS[cid-200]` | 단백질 far |
| 209–212 | `PLDDT_NEAR[cid-209]` | pLDDT near |
| 213–216 | `PLDDT_FAR[cid-213]` | pLDDT far |
| 217–226 | `CONSERVATION_NEAR[cid-217]` | 보존도 near |
| 227–236 | `CONSERVATION_FAR[cid-227]` | 보존도 far |
| 237 | `INTERFACE_NEAR_COLOR` (213) | 인터페이스 near |
| 238 | `INTERFACE_DIM_NEAR_COLOR` (243) | 비인터페이스 near |
| 239 | `INTERFACE_FAR_COLOR` (90) | 인터페이스 far |
| 240 | `INTERFACE_DIM_FAR_COLOR` (236) | 비인터페이스 far |
| 241–249 | `ALIGNED_BRIGHT_NEAR[cid-241]` | 정렬 bright near |
| 250 | `ALIGNED_NONALIGNED_FAR` (235) | 비정렬 어두운 far |
| 그 외 | `231` (흰색) | fallback |

### 확인 사항
- xterm-256 ID2RGBA LUT은 `Palette.hpp:204`에 이미 완전히 정의되어 있음 (256 RGBA 엔트리)
- ANSI 출력: `\033[38;5;{xterm_id}m` — xterm-256 전경색 직접 지정 가능
- 배경은 기본값 유지 (`-1` = 터미널 기본 배경)

---

## Appendix C: 공개 API 확정 초안

### 결정 사항
- **단일 구조 우선** (Foldseek 초기 사용 케이스)
- **전체 7종 색상 모드 지원** (protein, chain, rainbow, plddt, interface, conservation, aligned)
- **`show_structure=false` 기본** (Foldseek 사용 케이스; StructureMaker 의존성 제외 가능)
- **출력 형식**: ANSI escape string 반환 또는 stdout 직접 출력

### 확정 시그니처 (`include/structty_render.h`)

```cpp
#pragma once
#include <string>
#include <vector>

namespace structty {

struct RenderAtom {
    float x, y, z;
    char  structure;           // 'x'=coil, 'H'=helix, 'S'=sheet
    float bfactor;             // pLDDT 신뢰도 (0-100)
    bool  is_interface;
    bool  is_aligned;
    float conservation_score;  // -1.0 = 미설정, 0.0–1.0 = 보존도
    int   residue_number;
    char  residue_name[4];     // 3-letter code + null
    char  chain_id;
    int   protein_index;       // 0-8 (멀티 구조 지원 예약)
};

// ANSI escape string 반환 (Foldseek 임베딩용)
// width/height: 터미널 문자 단위
// mode: "protein"|"chain"|"rainbow"|"plddt"|"interface"|"conservation"|"aligned"
std::string render_to_ansi(
    const std::vector<RenderAtom>& atoms,
    int width,
    int height,
    const std::string& mode = "protein",
    bool show_structure = false
);

// stdout 직접 출력 (단일 호출 편의)
void render_to_stdout(
    const std::vector<RenderAtom>& atoms,
    int width,
    int height,
    const std::string& mode = "protein",
    bool show_structure = false
);

} // namespace structty
```

### Phase 1 구현 지침 (Phase 0 발견 결과 반영)

1. `palette_to_xterm256_fg(int color_id) -> int`는 Camera.cpp의 if-else 블록을 Palette.hpp로 이동하기만 하면 됨
2. `palette_to_ansi_fg_str(int color_id) -> std::string`은 위 함수 결과를 `"\033[38;5;" + to_string(n) + "m"` 포맷으로 감싸는 one-liner
3. `depth_band` 처리: near/mid는 색상 그대로, far는 동일 색상 그대로 (이미 Palette에서 어두운 xterm 번호로 사전 선택됨) — `\033[2m` dim escape 불필요
4. 렌더링 경계 교체는 `print_screen_braille()` 내 3줄(attron+mvaddstr+attroff)만 ANSI 경로로 분기하면 됨
