# FoldMason 결과 표시 연구

> 조사 목적: StrucTTY에서 FoldMason 다중 구조 정렬(MSTA) 결과를 시각화하기 위한 파일 포맷, 데이터 구조, 구현 방향 파악

---

## 1. FoldMason 개요

FoldMason은 Steinegger Lab이 개발한 **대규모 다중 단백질 구조 정렬(Multiple Protein Structure Alignment, MSTA)** 도구다.
- Foldseek의 3Di 구조 알파벳 기반으로 수십만 개 구조를 정렬
- Foldseek(pairwise) vs FoldMason(multiple): FoldMason은 여러 구조를 동시에 정렬하여 진화적 보존성을 분석
- 논문: Science (2026), bioRxiv 2024.08.01.606130
- 웹 서버: https://search.foldseek.com/foldmason

---

## 2. CLI 주요 모듈

| 모듈 | 설명 |
|------|------|
| `foldmason easy-msa <구조파일들> result tmpFolder` | 전체 워크플로 (DB 생성 → 정렬 → 보고서 생성) |
| `foldmason createdb <구조파일들> myDb` | 구조 파일을 내부 DB로 변환 |
| `foldmason structuremsa myDb result` | DB에서 MSTA 수행 |
| `foldmason msa2lddt myDb result_aa.fa` | 정렬 품질(LDDT 점수) 계산 |
| `foldmason msa2lddtreport myDb result_aa.fa result.html` | HTML 시각화 보고서 생성 |
| `foldmason msa2lddtjson myDb result_aa.fa result.json` | JSON 데이터 출력 |
| `foldmason refinemsa myDb result.fasta refined.fasta` | 반복적 MSA 개선 |

### easy-msa 주요 파라미터

| 파라미터 | 기본값 | 설명 |
|----------|--------|------|
| `--gap-open` | 10 | 갭 개시 패널티 |
| `--gap-extend` | 1 | 갭 연장 패널티 |
| `--refine-iters` | 0 | 초기 정렬 후 refinement 반복 횟수 |
| `--output-mode` | 0 | 0=아미노산, 1=3Di 알파벳 |
| `--report-mode` | 0 | 0=없음, 1=HTML, 2=JSON |
| `--precluster` | - | 대용량 데이터셋 사전 클러스터링 |

---

## 3. 출력 파일 포맷

### 3.1 기본 출력 파일

```
easy-msa 실행 후 생성되는 파일:
  result_aa.fa      — 아미노산 서열 MSA (FASTA 형식)
  result_3di.fa     — 3Di 구조 알파벳 MSA (FASTA 형식, 내부명 "ss")
  result.nw         — 가이드 트리 (Newick 형식)

--report-mode 1 추가 시:
  result.html       — 인터랙티브 HTML 시각화

--report-mode 2 추가 시:
  result.json       — 웹 서버용 JSON 데이터
```

### 3.2 FASTA 포맷 (\_aa.fa / \_3di.fa)

```fasta
>3A0C-assembly1cif_A
------DFKDK--PFGKAWAPGWHADPDAWIFHCHQQRKTFIDRNN---RTPDIQPLGNV-AGGWMWG...
>3A0C-assembly1cif_B
------DFKDK--AFGKAWAPGWHDDPDAWIFHCHQQRKTFIDDNN---RTPDIPPLGNV-AGRWMWT...
>L7RCY6pdb
DDDDDDDDDDDDDDDDDDDDDDDDDPPPDPPPPPPPPVPPDDDPPPPPDPPPPDQADPPPPPFAQQEEP...
```

- 각 시퀀스는 정렬된 위치를 공유하므로 **동일한 열 길이**
- 갭 문자: `-`
- `_aa.fa`: 표준 20가지 아미노산 1문자 코드 (A, C, D, ..., Y)
- `_3di.fa`: 3Di 구조 알파벳 문자 (아미노산 코드와 같은 문자 집합을 재사용하나 의미 다름)

### 3.3 3Di 구조 알파벳 (\_3di.fa)

FoldMason 내부에서 3Di 파일을 **"ss" (secondary structure)** 라고 부른다.
내부 소스: `seqDbr3Di`, `cigars_ss`, JSON의 `ss` 필드 모두 3Di 데이터를 가리킴.

- **3Di 알파벳 = 20개 상태**, Foldseek의 VQ-VAE 인코더로 학습된 구조적 상태
- 각 잔기에 대해 **공간적으로 가장 가까운 잔기 j와의 backbone 기하학**을 인코딩
- 아미노산 문자(A-Y 중 20개)를 재사용하여 표현 — 문자 자체가 아미노산을 뜻하지 않음
- 보존도가 높은 단백질 코어에 정보가 집중됨 (코일/루프 영역은 낮음)

**사용자 예시 (`foldmason_ss result`)의 해석:**

```
>L7RCY6pdb
DDDDDDDDDDDDDDDDDDDDDDDDDPPPDPPPPPPPPVPPDDDPPPPPDPPPPDQADPPPPP...
```

→ 이것은 `_3di.fa` 파일의 내용 (3Di 구조 알파벳)
→ `D`, `P`, `V`, `Q`, `A`, ... 등은 3Di 구조 상태를 나타내는 문자
→ 단백질이 D로 시작하는 긴 연속 구간 = 흔한 구조 상태(3Di 상태 D)가 연속적으로 배치된 것
→ DSSP의 H/E/C와 달리, 3Di는 tertiary interaction geometry를 인코딩

### 3.4 JSON 포맷 (`--report-mode 2`)

```json
{
  "entries": [
    {
      "name": "3A0C-assembly1cif_A",
      "aa":   "------DFKDK--PFGKAWAPG...",
      "ss":   "------DDPDP--DPPDDPDDD...",
      "ca":   [[x,y,z], ...]
    },
    {
      "name": "L7RCY6pdb",
      "aa":   "MKTLLLT...",
      "ss":   "DDDDDDD...",
      "ca":   [[x,y,z], ...]
    }
  ],
  "scores": [0.82, 0.79, 0.91, ...],
  "tree": "(3A0C-assembly1cif_A:0.1,L7RCY6pdb:0.2);",
  "statistics": {
    "db": "/path/to/myDb",
    "msaFile": "result_aa.fa",
    "msaLDDT": 0.85,
    "cmdString": "foldmason easy-msa ..."
  }
}
```

**필드 설명:**

| 필드 | 타입 | 설명 |
|------|------|------|
| `entries[].name` | string | 구조 헤더 (파일명 또는 PDB ID) |
| `entries[].aa` | string | 정렬된 아미노산 서열 (갭 포함) |
| `entries[].ss` | string | 정렬된 3Di 구조 알파벳 (갭 포함) |
| `entries[].ca` | array | C-alpha 좌표 배열 (옵션) |
| `scores` | float[] | 각 정렬 열의 LDDT 점수 |
| `tree` | string | Newick 형식 가이드 트리 (옵션) |
| `statistics.msaLDDT` | float | 전체 MSA LDDT 점수 |

---

## 4. LDDT 품질 점수

FoldMason은 MSA 품질 측정에 **LDDT (Local Distance Difference Test)** 를 사용:

- 모든 구조 쌍에 대해 정렬된 잔기의 pairwise LDDT 계산
- 갭이 있는 열은 계산에서 제외
- **per-column LDDT**: 각 정렬 열의 평균 품질 → JSON `scores` 배열
- **msaLDDT**: 전체 정렬의 평균 LDDT 점수 (0~1)
- StrucTTY 시각화 시 per-column LDDT로 열별 보존도 색상 표시 가능

---

## 5. FoldMason vs Foldseek 차이점

| 항목 | Foldseek | FoldMason |
|------|----------|-----------|
| 정렬 유형 | Pairwise (1:1) | Multiple (N:N) |
| 출력 형식 | `.m8` 탭 구분 테이블 | FASTA MSA + 트리 + JSON |
| 쿼리 | 단일 구조 | 여러 구조 동시 입력 |
| 결과 | 유사 구조 목록 + 정렬 정보 | 공통 정렬 프레임의 MSA |
| 시각화 | 없음 (텍스트 테이블) | HTML 인터랙티브 리포트 |
| 용도 | 구조 검색 | 진화 보존성 분석, 계통수 |

---

## 6. StrucTTY 구현 방향 (사전 조사)

### 6.1 입력 파일 타입

```
--foldmason result.json    → JSON 파일 로드 (entries + scores + tree)
--foldmason result_aa.fa   → FASTA MSA 파일 로드 (aa 서열만)
```

JSON이 aa/ss/ca/scores 모두 포함하므로 JSON 파싱을 우선 지원하는 것이 유리.
FASTA 파일만 있을 경우 `_aa.fa` + `_3di.fa` 쌍을 함께 받거나 aa만 처리.

### 6.2 파싱 대상 데이터

- `entries[i].name`: 각 구조 식별자 → Foldseek hit 탐색의 target ID와 동일 형식
- `entries[i].aa`: 갭 포함 아미노산 정렬 → 잔기-열 매핑용
- `entries[i].ss`: 갭 포함 3Di 구조 알파벳 → 열별 구조 보존도 계산 가능
- `entries[i].ca`: C-alpha 좌표 → 구조 겹침(superposition) 가능
- `scores[]`: per-column LDDT → 색상 표시 기준

### 6.3 시각화 아이디어

1. **MSA 컬럼 뷰어**: 수평 스크롤로 정렬 열 탐색, 잔기별 conservation 색상
2. **per-column LDDT 막대 그래프**: JSON `scores` 배열로 각 열의 구조 보존도 표시
3. **3Di 보존도 표시**: ss 컬럼별 다양도(entropy)로 구조적 보존 정도 시각화
4. **N/P 키로 시퀀스 전환**: Foldseek hit 탐색처럼 entries 간 전환

### 6.4 구현 시 주의사항

- `_aa.fa`와 `_3di.fa`는 **동일한 열 길이**를 가져야 함 (갭 포함)
- 갭(`-`) 위치는 C-alpha 좌표가 없으므로 `ca` 배열 인덱스 주의
- `ca` 필드는 옵션 — 없을 경우 구조 겹침 불가, 서열 기반 시각화만
- entries 수가 많을 경우(수백~수천) 렌더링 성능 고려 필요
- JSON 파일은 수십 MB가 될 수 있음 → 스트리밍 파서 검토

### 6.5 신규 파서 클래스 설계안

```cpp
// src/structure/FoldMasonParser.hpp

struct FoldMasonEntry {
    std::string name;
    std::string aa;       // gap-included aa alignment
    std::string ss;       // gap-included 3Di alignment
    std::vector<std::array<float,3>> ca;  // C-alpha coords (optional)
};

class FoldMasonParser {
public:
    // JSON 파일 파싱
    static bool parse_json(const std::string& path,
                           std::vector<FoldMasonEntry>& entries,
                           std::vector<float>& scores,
                           std::string& tree,
                           std::string& error_msg);

    // FASTA MSA 파일 파싱 (_aa.fa or _3di.fa)
    static bool parse_fasta(const std::string& path,
                            std::vector<FoldMasonEntry>& entries,
                            std::string& error_msg);

    // 열별 Shannon entropy 계산 (aa 또는 ss 컬럼 기반)
    static std::vector<float> compute_column_entropy(
        const std::vector<FoldMasonEntry>& entries,
        bool use_ss = false);
};
```

---

## 7. 웹 서버 API (참고)

- URL: https://search.foldseek.com/foldmason
- POST로 PDB 파일 업로드 또는 기존 JSON 결과 업로드 가능
- 결과를 JSON으로 다운로드하여 CLI와 동일한 포맷으로 로컬 처리 가능

---

## 8. 참고 자료

- [FoldMason GitHub (steineggerlab/foldmason)](https://github.com/steineggerlab/foldmason)
- [FoldMason README](https://github.com/steineggerlab/foldmason/blob/main/README.md)
- [Science 논문 (2026)](https://www.science.org/doi/10.1126/science.ads6733)
- [bioRxiv 프리프린트](https://www.biorxiv.org/content/10.1101/2024.08.01.606130v2.full)
- [FoldMason 웹 서버](https://search.foldseek.com/foldmason)
- [Foldseek GitHub (구조 알파벳 원본)](https://github.com/steineggerlab/foldseek)
