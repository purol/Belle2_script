# FEI 선택 단계 및 BCS 추적

## Category 2·3 원인 추적용 확장 실행 (2026-11-20)

새 패키지 `fei_category23_diagnosis`에서는 다음 명령을 사용합니다.
이전에 완료한 BCS 실행을 덮어쓰지 않고 별도 출력 폴더를 만듭니다.

```bash
cd fei_category23_diagnosis
python3 -u trace_fei_stages.py --scope full --details --workers 4 --output fei_category23_trace
```

아래의 일반 `full` 실행에 `--details`를 추가한 것입니다. 같은 7개 ROOT snapshot을
한 번씩 읽으며, 상세 집계를 위해 원본 ROOT를 다시 열지 않습니다.
모든 24,860개 파일을 확인하고, 최종 생존 사건이 없는 파일도 포함합니다.
입력은 `rb`로 열고 새 출력 폴더에만 씁니다. 기존 출력 폴더를 지정하면 중단합니다.
실행 중에는 `RUNNING`, 정상 완료 시에는 `COMPLETE`가 생깁니다.
검사 문제가 있으면 `REVIEW_REQUIRED`가 생기며 물리적 원인으로 해석하면 안 됩니다.

추가 결과:

- `category23_cutflow.csv`: Category 2·3의 단계별 Data/MC 수율. MC는 원래 표본 정규화만 적용.
- `category23_cut_retention.csv`: 같은 사건을 연결한 각 cut의 Data/MC 잔류율과 상대 잔류율.
- `category23_diagnostic_summary.txt`: 관측된 비율 감소가 가장 큰 구간. 유의도나 물리적 원인 판정은 아님.
- `selection_event_trace.csv.gz`: 각 단계의 후보 수, tag 모드, 최고 FEI 점수,
  다른 category 경쟁 여부, correct/incorrect/unknown truth 구성.
- `selection_transitions.csv.gz`: 어떤 cut에서 대상 후보를 모두 잃었는지,
  사건 전체가 사라졌는지 또는 다른 후보만 남았는지 구분.
- `selection_candidate_trace.csv.gz`: 첫 저장 단계에서 Category 2 또는 3 control 후보가
  있었던 사건의 모든 후보를 이후 단계까지 저장. 대상 후보가 사라진 뒤의 경쟁 후보도 보존.
  Mbc, deltaE, Eecl, missing theta, signal momentum, tag mode, FEI 점수와
  저장돼 있는 경우의 ROE track/pi0/KS 개수, matched/unmatched Eecl 등을 포함.
- `selection_components.csv`: 단계별 경쟁 여부·최고 점수 tag 모드·MC truth의 교차 집계.
  최고 점수가 여러 tag 모드에서 같으면 모드 집합을 저장하며 임의로 하나를 고르지 않음.
- `extra_branch_coverage.csv`: 추가 관측량의 파일별 존재 여부를 표본/단계별로 집계.
  없는 값은 NaN이며 측정값 0으로 취급하지 않음.
- `representative_branch_inventory.csv`: 각 표본의 첫 파일에 대한 tree/branch 목록.
  전체 생산 설정이 같다는 증거로 취급하지 않음.

상세 집계는 두 category 각각, `control_K_Kpi`와 `mode0` 각각에 대해 기존 전체
cutflow의 사건 수·정규화 수율·후보 수와 일치하는지 검사합니다. 인접 단계의 사건
연결도 앞뒤 snapshot 모두와 일치해야 합니다. Category/view 사이에는 사건이 겹칩니다.
각 파일의 모든 후보를 함께 읽으므로 대상 category 밖의 경쟁 후보를 누락하지 않습니다.
BCS는 원래 최대 점수 규칙과 저장된 실제 winner를 검사하며 다시 무작위 선택하지 않습니다.

이 실행으로 선택 단계와 후보/MC 성분을 좁힐 수 있지만, 첫 저장 단계 이전의
재구성·ROE/vertex 선택, FEI 학습 및 생산 설정, 저장되지 않은 daughter track/cluster
식별자는 복원할 수 없습니다. 검출효율/붕괴모델의 특정 물리적 원인을 자동 확정하지 않습니다.
Tag truth는 대상 후보 중 correct가 하나라도 있는지를 구분하며, 실제 winner의
truth 또는 signal 쪽 J/psi의 진위를 뜻하지 않습니다.

Python 3.8 이상에서 실행하도록 작성했으며 numpy/pandas/uproot 외 추가 분석 패키지는
필요하지 않습니다. 테스트에만 pytest가 필요합니다. 상세 출력과 읽는 단계가 늘어나므로
이전 BCS 전용 실행보다 오래 걸립니다. 출력은 gzip으로 압축하며 중간 재시작은 지원하지 않습니다.

완료 후 `fei_category23_trace` 폴더 전체를 가져오면 원본 ROOT를 다시 읽지 않고
cut·tag 모드·경쟁 여부·truth에 따른 원인을 교차 확인할 수 있습니다.

2026-11-20에 받은 실제 `Analysis_main.cc`, `Loader.h`와 중간 출력 폴더 구조를 기준으로 작성했습니다.
기존 ROOT 파일을 읽어 조사하며, 원본 분석 코드나 ROOT 파일을 수정하지 않습니다.
ROOT/basf2를 다시 실행하거나 새 ntuple을 만들 필요가 없습니다.

## 서버에서 실행

빠른 버전의 압축을 풀고 `fei_stage_trace_fast` 폴더에서 실행합니다. Python 3.8 이상을 사용합니다.
먼저 BCS 전후에 집중하는 실행을 권장합니다.

```bash
python3 -m pip install -r requirements.txt
python3 -u trace_fei_stages.py --scope bcs --workers 4 --output fei_bcs_trace_fast
```

현재 서버 경로는 코드에 이미 지정되어 있습니다.

- MC: `/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_Jpsi`
- Data: `/home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_data_Jpsi`
- 버전: `validation_v009`

`reference/`에는 앞서 제공한 실제 FEI export의 `files.csv`, `selected_events.csv`, `reference.csv`가 포함되어 있습니다.
이 manifest의 파일별 정규화를 사용하고, 같은 원본 파일에 해당하는 중간 출력끼리 연결합니다.
해당 export의 24,860개 파일을 조사합니다. 원래 필요한 파일이 누락됐는지 확인하려고
최종 선택 사건이 0개였던 파일도 읽습니다. 최종 생존 사건이 있는 파일만 골라 읽지 않습니다.
진행 상황은 100개 원본 파일마다 출력됩니다.

`--scope bcs`는 원본 파일당 `before_Dveto_cut`과 `final_output` 2개를 읽습니다.
`--workers 4`는 최대 4개의 읽기 프로세스를 사용합니다. 부모 프로세스만 결과를 기록하며,
입력 ROOT 파일은 모든 프로세스에서 읽기 전용입니다.
동시에 준비하는 원본 파일 수는 worker 수의 2배로 제한하고 결과를 원래 순서대로 모읍니다.

전체 cut 단계까지 조사하려면 다음 명령을 사용합니다.

```bash
python3 -u trace_fei_stages.py --scope full --workers 4 --output fei_full_trace_fast
```

`--scope full`은 저장된 7개 단계 전체를 읽고 같은 BCS 진단도 포함합니다.
BCS 결과가 그 이전 문제를 가리킬 때 전체 단계를 조사하면 됩니다.
`--workers`를 생략하면 순차 실행(1개), `--scope`를 생략하면 전체 단계 실행입니다.
출력 폴더가 이미 있으면 덮어쓰지 않고 중단하므로, 재실행 시 다른 출력 이름을 사용합니다.

### 이전 버전이 실행 중인 경우

새 압축은 다른 폴더(`fei_stage_trace_fast`)에 풀어 기존 실행 파일을 보존합니다.
전환하기로 했다면 기존 터미널에서 Ctrl+C로 중단한 후 위의 새 명령을 실행합니다.
중단은 읽던 원본 ROOT를 변경하지 않지만, 기존 결과 폴더는 미완료 상태로 남습니다.
현재 버전은 중간 결과에서 이어서 실행하는 기능이 없으며, 새 실행은 처음부터 읽습니다.
기존 출력 폴더를 지우거나 같은 출력 이름을 재사용할 필요는 없습니다.

프로세스 수와 읽는 snapshot 수를 줄이는 효과는 서버의 CPU/파일시스템 상태에 따라 달라집니다.
정확한 속도 향상 배율은 보장하지 않습니다. 로그의 rough_remaining은 현재 평균 속도를
단순 환산한 값으로, 다른 표본으로 넘어가면 달라질 수 있습니다.

파일 목록만 빠르게 확인하려면 다음 명령을 사용합니다. ROOT 내용을 읽기 전에 경로/파일 대응을 확인합니다.

```bash
python3 trace_fei_stages.py --preflight-only
```

경로를 옮겼다면 `--mc-base`, `--data-base`를 지정할 수 있습니다.
다른 표본/버전에 적용할 경우에는 그 실행에서 생성한 FEI export를 `--reference`로 주어야 합니다.
폴더 이름만 바꿔 서로 다른 생산 버전의 파일을 연결해서는 안 됩니다.

## 무엇을 비교하는가

| 저장된 단계 | 다음 단계까지 적용되는 조건 |
|---|---|
| before_Mbc_cut | tag Mbc > 5.27 |
| before_delE_cut | -0.2 < tag deltaE < 0.2 |
| before_Eecl_cut | Eecl_v200 < 1.3 |
| before_missing_momentum_theta_cut | 0.297 < theta_missing < 2.618 |
| before_psig_cut | 0.5 < p_sig(CMS) < 2.96 |
| before_Dveto_cut | 1.84 < Bsig_M < 1.89 후보 제거, 이후 BCS |
| final_output | 실제 최종 BCS 및 IsBCSValid 처리 결과 |

`before_Dveto_cut`에 원본과 동일한 엄격 부등식을 적용해 `after_Dveto_reconstructed`를 만듭니다.
다음으로 모든 후보를 경쟁시켜 `after_max_probability_reconstructed`를 만듭니다.
Category 3이나 Mode 0를 먼저 골라서 BCS를 수행하지 않습니다.
최대 확률 비교는 Loader가 사용하는 연속된 `(experiment, run, event, ncandidates)` 블록 및 정확한 동률 조건을 재현합니다.
물리적 사건 집계에는 파일 내 `(experiment, run, event)`를 사용하며,
같은 사건의 ncandidates 값이 여러 종류이거나 블록이 흩어져 있으면 따로 보고합니다.

무작위 BCS의 seed/hash를 Python에서 새로 만들지 않습니다. `final_output`에 저장된 실제 후보를 읽습니다.
따라서 다음을 구분할 수 있습니다.

- `kept_target`: 해당 category/mode에 남음.
- `lost_at_max_probability`: 더 높은 점수의 다른 후보가 있어서 탈락.
- `lost_among_equal_maxima`: 동률 최대 후보에는 포함됐지만 실제 선택된 후보는 다른 category/mode.
- `audit_issue_no_unique_valid_winner`: 최종 후보가 유일하지 않거나 최대 후보와 대응하지 않는 등 추가 검토 필요.

`final_output`에는 원본 `IsBCSValid`의 data 중복 처리도 포함됩니다. 중복/불연속 블록/최종 사건 손실 등이 있으면
이를 정상적인 무작위 탈락으로 단정하지 않고 audit issue로 표시합니다.

## 결과를 읽는 순서

1. `summary.txt`와 `COMPLETE` 또는 `REVIEW_REQUIRED`를 확인합니다.
2. `data_mc_stage_ratios.csv`에서 `category=3`, `view=mode0`를 봅니다.
   어떤 단계에서 비율이 낮아지는지 Category 1·2와 비교합니다.
3. `bcs_migration_summary.csv`에서 `origin_category=3`, `view=mode0`를 봅니다.
   최대 점수 경쟁에서 밀리는지, 동률 선택에서 밀리는지, 선택된 category/mode가 무엇인지 봅니다.
4. `bcs_event_trace.csv.gz`와 `bcs_candidate_trace.csv.gz`에서 해당 사건의 경쟁 후보들을 확인합니다.
   파일 ID는 `source_files.csv`를 통해 실제 ROOT 파일 경로와 연결됩니다.

`stage_summary.csv`는 후보 수와 중복을 제거한 사건 수를 모두 제공합니다.
BCS 전에는 한 사건에 여러 category/mode 후보가 있을 수 있어 각 category의 사건 수를 더하면 중복됩니다.
BCS migration의 origin cohort도 서로 겹칠 수 있습니다. 각 origin 내의 행들을 비교해야 합니다.

### 가중치와 오차

단계 비교에는 **기존 manifest의 표본 정규화만** 적용합니다. Detector/FEI 보정은 추가하지 않습니다.
같은 후보의 가중치를 단계마다 바꿔서 생기는 차이를 피하기 위한 정의입니다.
따라서 최종 비율은 기존 보정계수 0.58/0.23 자체와 다릅니다.
실제 제공된 export가 재현된다면 normalization-only 비율은 다음 값과 대응합니다.

- Neutral: 0.9127
- Category 1: 0.8299
- Category 2: 0.5498
- Category 3: 0.2167

단계별 통계 오차는 사건 수와 정규화 가중치의 sumw2로 계산합니다.
0개 data의 대칭 오차는 0으로 쓰지 않고 빈 값으로 둡니다.
같은 사건을 재사용하는 단계/영역은 서로 상관되어 있으며, systematic이나 유의도는 계산하지 않습니다.

## 일관성 검사와 한계

- Tree 길이와 experiment/run/event/candidate/ncandidates를 비교합니다. 다른 후보를 행 번호만으로 붙이지 않습니다.
- 필수 branch와 유효한 값, probability 범위를 확인합니다. 선택된 파일만 검사하는 방식이 아닙니다.
- 각 단계에서 원래 cut을 적용한 후보 목록과 실제 다음 snapshot의 후보 목록/변수 값을 비교합니다.
- 최종 후보가 D-veto 후의 실제 후보 및 최대 점수 후보에 속하는지 확인합니다.
- final_output의 control 선택 후보가 이전 FEI export의 선택 후보 식별자 및 개수와 일치하는지 확인합니다.
- 파일 목록이 reference와 다르면 중단합니다. Merge 하위 폴더는 읽지 않습니다.
- Data의 같은 사건 ID가 서로 다른 입력 파일에 겹치면 표시합니다. MC는 원본 파일 namespace를 유지하여 다른 생산의 ID를 합치지 않습니다.
- 선택 전 사건이 0개인 파일도 처리합니다. 일치하지 않는 파일을 조용히 제외하지 않습니다.

`COMPLETE`는 모든 파일 처리와 검사가 통과한 경우입니다.
`REVIEW_REQUIRED`가 나오면 `audit_issues.csv`를 먼저 검토해야 하며 수율/이동 결과는 잠정적입니다.
오류로 중단되면 `RUNNING`만 남고 완료 표시가 생기지 않습니다. 이 경우 터미널 오류와 생성된 결과를 보존해 주세요.

가장 이른 `before_Mbc_cut`도 최초 ROE/vertex 조건과 그 앞의 Jpsi 선택/공식 skim 이후입니다.
이 단계부터 이미 차이가 있으면 이번 파일만으로 그 앞의 손실까지 복원할 수는 없습니다.
관측된 손실 단계와 실제 물리적 원인은 구분해서 해석해야 합니다.

## 검증

```bash
python3 -m pip install pytest
python3 -m pytest test_trace_fei_stages.py -q
```

실제 형식의 모의 ROOT 파일로 전체 입출력과 선택 단계, candidate 이동을 확인하는 테스트를 포함합니다.
실제 서버 파일을 이 Windows 환경에서 읽은 것은 아닙니다. 제공된 파일명 규칙이 서버 파일명과 다르면
첫 preflight에서 예시 파일명을 포함한 오류가 나오며 임의 추측으로 연결하지 않습니다.

ROOT 읽기는 [Uproot 공식 문서](https://uproot.readthedocs.io/en/stable/basic.html)의 file-handle 및 array API를 사용합니다.
