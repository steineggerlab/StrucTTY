# Plan: Multimer 결과 시각화 선택 단위를 chain → complex로 변경

## 목표 및 범위

`foldseek easy-multimer-search` 실행 후 StrucTTY를 열었을 때 숫자키(1/2/0)로 선택하는 단위를 **chain 단위**(12AS_A, 12AS_B)에서 **complex 단위**(12AS)로 변경한다.

현재 문제는 두 층에 걸쳐 있다:

1. **CLI 진입 불가**: StrucTTY에는 이미 `multimer_mode_`와 `set_multimer_report()` 인프라가 구현되어 있으나, `Parameters`에 `--report-format` / `--query-db` 플래그가 없어 CLI에서 절대로 진입할 수 없다. `opts.report_format`이 항상 `false`이므로 `structty.cpp:57`의 분기가 실행되지 않는다.

2. **숫자키 complex 단위 미지원**: 현재 `Screen::handle_input_impl`의 숫자키 핸들러는 `structNum = key - 49`로 `data[]` 인덱스를 직접 사용한다. multimer mode에서 `data[]`는 query/target의 체인 단위로 채워지므로, "1" = data[0] = 첫 번째 query 체인이 된다. complex 단위로 동작하려면 "1" = query 체인 전체(data[0..mm_query_chain_count_-1])가 되어야 한다.

**범위 밖**: Foldseek 쪽 `easymultimersearch.sh` 변경(별도 plan.md — `../foldseek_work/plan.md` 참조), FoldMason, MSA, 일반 monomer 검색 워크플로.

---

## 영향받는 파일 목록

| 파일 | 변경 유형 | 이유 |
|------|-----------|------|
| `src/structure/Parameters.hpp` | 수정 | `--report-format`, `--query-db` 필드 및 getter 추가 |
| `src/structure/Parameters.cpp` | 수정 | 두 CLI 플래그 파싱 추가, input file 필수 조건 완화 |
| `src/main.cpp` | 수정 | `opts.report_format`, `opts.foldseek_query_db` 연결 |
| `src/visualization/Screen.cpp` | 수정 | multimer mode에서 숫자키 → complex 범위 선택으로 변경 |

---

## 구현 단계

- [x] 단계 1: `Parameters` — `--report-format` / `--query-db` CLI 플래그 추가
- [x] 단계 2: `main.cpp` — 새 Parameters 필드를 `RunOptions`에 연결
- [x] 단계 3: `Screen.cpp` — multimer mode에서 숫자키 complex 단위 선택 구현
- [ ] 단계 4: 수동 검증

---

## 단계별 세부 작업

### 단계 1: Parameters — CLI 플래그 추가

**`src/structure/Parameters.hpp`에 추가할 내용:**

```cpp
// private 멤버 추가
bool report_format = false;       // --report-format: foldseek_file이 _report TSV
std::string foldseek_query_db = ""; // --query-db: query complex DB 경로

// public getter 추가
bool get_report_format() { return report_format; }
std::string get_foldseek_query_db() { return foldseek_query_db; }
```

**`src/structure/Parameters.cpp`에 추가할 내용:**

`print_help()`에 두 옵션 문서 추가:
```cpp
std::cout << "  --report-format         Input is a Foldseek multimer _report file (14-col TSV)\n";
std::cout << "  --query-db <PATH>       Query complex Foldseek DB (required with --report-format)\n";
```

파싱 루프에 추가:
```cpp
} else if (!strcmp(argv[i], "--report-format")) {
    report_format = true;
} else if (!strcmp(argv[i], "--query-db")) {
    if (i + 1 < argc) {
        foldseek_query_db = argv[++i];
    } else {
        throw std::runtime_error("Error: Missing value for --query-db.");
    }
}
```

input file 필수 조건 완화 (`Parameters.cpp:115`):
```cpp
// Before
if (in_file.size() == 0) {
    std::cerr << "Error: Need input file dir" << std::endl;
    arg_okay = false;
    return;
}

// After
if (in_file.size() == 0 && !report_format) {
    std::cerr << "Error: Need input file dir" << std::endl;
    arg_okay = false;
    return;
}
```

`--report-format` 사용 시 `--foldseek`와 `--query-db` 모두 필요 → `Parameters.cpp` 말미에 검증 추가:
```cpp
if (report_format && (foldseek_file.empty() || foldseek_query_db.empty())) {
    std::cerr << "Error: --report-format requires both --foldseek <_report> and --query-db <DB>." << std::endl;
    arg_okay = false;
}
```

### 단계 2: main.cpp — RunOptions 연결

```cpp
// 기존 코드 뒤에 추가 (RunOptions 구조체 필드에 대응)
opts.report_format    = params.get_report_format();
opts.foldseek_query_db = params.get_foldseek_query_db();
```

`RunOptions.report_format`과 `RunOptions.foldseek_query_db`는 `include/structty.h`에 이미 선언되어 있으므로 헤더 변경 불필요.

### 단계 3: Screen.cpp — multimer mode 숫자키 complex 단위 선택

**private 헬퍼 선언** (`Screen.hpp` private 섹션):
```cpp
// multimer mode에서 complex index → data[] 범위 반환
// complex_idx=0: query, complex_idx=1: target. 반환값 {lo, hi} (inclusive)
std::pair<int,int> mm_complex_range(int complex_idx) const;
```

**헬퍼 구현** (`Screen.cpp`):
```cpp
std::pair<int,int> Screen::mm_complex_range(int complex_idx) const {
    if (complex_idx == 0) {
        return {0, mm_query_chain_count_ - 1};
    } else {  // complex_idx == 1: target
        return {mm_query_chain_count_, (int)data.size() - 1};
    }
}
```

**숫자키 핸들러 변경** (`Screen.cpp` — `handle_input_impl`, 현재 `case 48..57` 블록):

```cpp
// select protein / complex
case 48:  // '0' — deselect
    structNum = -1;
    break;
case 49:
case 50:
case 51:
case 52:
case 53:
case 54:
case 55:
case 56:
case 57:
    if (multimer_mode_) {
        // multimer mode: 1 = query complex (idx 0), 2 = target complex (idx 1)
        int requested = key - 49;  // 0-based
        if (requested < 2) structNum = requested;
        // 3 이상은 현재 복합체 수 초과 — 무시
    } else {
        if (key - 48 <= (int)data.size()) {
            structNum = key - 49;
        }
    }
    break;
```

**이동/회전 연산에 complex 범위 적용** — `A/D/W/S` 키와 rotation 키의 `structNum` 사용 패턴을 multimer mode에서 루프로 대체:

현재 패턴 (예: A 키):
```cpp
if (structNum != -1) apply_pan(structNum, -pan_step_x, 0.0f);
else for (int i = 0; i < (int)data.size(); i++) apply_pan(i, -pan_step_x, 0.0f);
```

변경 후 (해당 케이스 전부 동일 패턴으로 적용):
```cpp
if (multimer_mode_ && structNum >= 0) {
    auto [lo, hi] = mm_complex_range(structNum);
    for (int i = lo; i <= hi; i++) apply_pan(i, -pan_step_x, 0.0f);
} else if (structNum != -1) {
    apply_pan(structNum, -pan_step_x, 0.0f);
} else {
    for (int i = 0; i < (int)data.size(); i++) apply_pan(i, -pan_step_x, 0.0f);
}
```

rotation 키 (X/Y/Z 등):
```cpp
// 현재
if (structNum != -1) data[structNum]->set_rotate(1, 0, 0);

// 변경 후
if (multimer_mode_ && structNum >= 0) {
    auto [lo, hi] = mm_complex_range(structNum);
    for (int i = lo; i <= hi; i++) data[i]->set_rotate(1, 0, 0);
} else if (structNum != -1) {
    data[structNum]->set_rotate(1, 0, 0);
}
```

적용 위치 전체 목록 (Screen.cpp):
- `apply_pan` 사용처: A(1153), D(1159), W(1165), S(1171)
- `set_rotate` 사용처: X/x 회전(1178-1190), Y/y 회전(1191-1203), Z/z 회전(1204-1216) 등 — 실제 라인은 빌드 후 확인

### 단계 4: 수동 검증

아래 테스트 계획 참조.

---

## 인터페이스 변경사항

### `Parameters.hpp` — 필드 및 getter 추가

```cpp
// Before (private)
string foldseek_db = "";

// After (private, 아래 추가)
bool report_format = false;
std::string foldseek_query_db = "";

// After (public getter, 아래 추가)
bool get_report_format() { return report_format; }
std::string get_foldseek_query_db() { return foldseek_query_db; }
```

### `Screen.hpp` — private 헬퍼 선언 추가

```cpp
// Step 7 multimer 헬퍼에 추가
std::pair<int,int> mm_complex_range(int complex_idx) const;
```

### CLI 인터페이스 — 새 플래그

| 플래그 | 타입 | 설명 |
|--------|------|------|
| `--report-format` | flag (bool) | `--foldseek`로 전달된 파일이 `_report` 14-컬럼 TSV임을 명시 |
| `--query-db <PATH>` | string | query complex Foldseek DB 경로 (`--report-format`과 함께 필수) |

### Foldseek 측 연계 변경 (별도 plan)

StrucTTY 측 변경이 완료되면 `../foldseek_work/data/easymultimersearch.sh`의 StrucTTY 호출을 아래로 변경해야 한다 (해당 내용은 `../foldseek_work/plan.md`를 업데이트할 것):

```sh
# Before
StrucTTY "${QUERY_INPUT}" "${TARGET_INPUT}" -ut "${UT_FILE}" --foldseek "${OUTPUT}" --db "${TARGET}"

# After  (-ut TSV 불필요, _report를 직접 전달)
StrucTTY --foldseek "${OUTPUT}_report" --report-format --query-db "${QUERY}" --db "${TARGET}"
```

- `-ut` TSV 생성 awk 블록 전체 제거 (U/T는 `_report` 파일에 내장됨)
- `QUERY_INPUT`/`TARGET_INPUT` positional 인자 제거 (DB에서 직접 로드)
- `QUERY` = foldseek query DB 경로 (chain createdb 후의 `${TMP_PATH}/query` 또는 원본 DB)

---

## 테스트 계획

| 테스트 케이스 | 검증 항목 | 유형 |
|---|---|---|
| `--report-format` 없이 기존 호출 | 기존 동작 유지 (no regression) | 수동 |
| `--report-format`만 있고 `--query-db` 없을 때 | 오류 메시지 출력 후 종료 | 수동 |
| `--report-format --foldseek <_report> --query-db <DB> --db <DB>` | multimer_report 분기 진입, query/target 체인 모두 로드 | 수동 |
| 숫자키 "1" 누름 (multimer mode) | query 복합체 전체 체인(12AS_A, 12AS_B 등) 함께 이동/회전 | 수동 |
| 숫자키 "2" 누름 (multimer mode) | target 복합체 전체 체인 함께 이동/회전 | 수동 |
| 숫자키 "0" 누름 (multimer mode) | 선택 해제, 전체 이동/회전 | 수동 |
| 숫자키 "1" 누름 (monomer mode) | 기존 동작 유지 (data[0] 단일 선택) | 수동 |
| n/p 키로 hit 탐색 (multimer mode) | 복합체 단위로 target hit 전환 (target 체인 모두 교체) | 수동 |
| `]/[` 키로 query 전환 (multimer mode) | 다른 query complex로 전환 | 수동 |
| `mm_complex_range(0)` 단위 테스트 | `{0, mm_query_chain_count_-1}` 반환 | unit |
| `mm_complex_range(1)` 단위 테스트 | `{mm_query_chain_count_, data.size()-1}` 반환 | unit |

---

## 예상 리스크 및 주의사항

- **`in_file.size() == 0` 완화**: `--report-format` 시 input file 없이 실행 가능하게 되므로, `Screen::set_protein`이 호출되지 않는다. `data[]`가 비어 있는 상태로 이벤트 루프에 진입하는 경로가 없어야 한다. `set_multimer_report`가 실패할 경우 `multimer_mode_ = false`로 복귀하므로, 이 경우 빈 화면이 표시된다. → 추가 방어 코드 필요: multimer mode 진입 실패 시 사용자에게 명확한 오류 메시지 출력 후 종료.

- **`data[]` 경계 검사**: `mm_complex_range`가 반환하는 `hi` 인덱스가 `data.size()-1`을 초과하지 않도록 보장해야 한다. target hit가 로드되지 않은 상태(empty target DB)에서 "2"를 누르면 `mm_query_chain_count_ > data.size()` 상황이 발생할 수 있다 → 헬퍼에서 `hi >= (int)data.size()` 조건 방어.

- **Foldseek `QUERY` DB 경로 확인**: `easymultimersearch.sh`에서 `--query-db "${QUERY}"`로 전달되는 `QUERY`는 createdb 실행 결과(`${TMP_PATH}/query`)인데, `REMOVE_TMP`가 활성화된 경우 StrucTTY 실행 전에 DB가 삭제될 수 있다. Foldseek 쪽에서 StrucTTY 호출을 cleanup 이전에 배치하는지 확인 필요 (현재 코드 상 StrucTTY 호출 블록이 `REMOVE_TMP` 블록보다 앞에 있어 문제 없음 — `easymultimersearch.sh:59`).

- **`structty.h`의 `RunOptions` 변경 불필요**: `report_format`과 `foldseek_query_db` 필드는 이미 `include/structty.h:22-23`에 선언되어 있다. 공개 헤더 변경 없음.

- **rotation 연산의 축 일관성**: 체인 간에 `set_rotate` 호출이 누적되지 않도록 주의 필요. Foldseek complex는 이미 U/T로 superpose된 상태이므로, 사용자가 complex 전체를 회전할 때 각 체인의 회전 행렬이 독립적으로 쌓이면 complex가 분해된다. → complex 선택 시 모든 체인에 동일한 회전 delta를 적용하는 현재 패턴이 올바르다 (각 chain이 동일한 set_rotate 호출을 받음).
