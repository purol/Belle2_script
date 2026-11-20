# FEI 계수 원인 진단

사용자가 제공한 `20261113/Jpsi/THStack_Jpsi_FEI_efficiency.c`를 기준으로 작성했습니다.
원본의 선택 조건, 입력 경로, 보정 가중치 계산식, Category 구분을 유지하고 진단 출력을 추가했습니다.
원본 스냅샷은 `source/`에 있으며 SHA-256은
`b0a7b899aa0bb63289157b840500f01b8fe9affa9b3276857d639606cf40768b`입니다.

## 실행

1. 아래 파일을 **원본 ROOT 매크로를 실행하던 작업 디렉터리**에 복사합니다.
   - `THStack_Jpsi_FEI_diagnostics.C`
   - `FEIDiagnosticExport.h`
   - `analyze_fei_diagnostics.py`
   - `requirements.txt`
2. 원본과 같은 ROOT 환경 및 `constants.h`, `ObtainWeight.h`, `correctors.h`, `base.h`,
   보정 데이터 파일을 사용합니다. 다른 버전의 보정 파일로 교체하지 마세요.
3. 다음 명령을 실행합니다.

```bash
root -l -b -q 'THStack_Jpsi_FEI_diagnostics.C("fei_diagnostics_output")'
python3 analyze_fei_diagnostics.py fei_diagnostics_output
```

Python 의존성이 없다면 해당 분석 환경에서 다음 명령으로 설치합니다.

```bash
python3 -m pip install -r requirements.txt
```

ROOT 매크로의 헤더를 별도 디렉터리에 두는 환경이라면 기존의 ROOT include 설정을 그대로 사용합니다.
`ROOT_INCLUDE_PATH`에 그 디렉터리를 추가해도 됩니다. 보정 파일을 상대 경로로 읽는 기존 코드 때문에
작업 디렉터리는 원본 실행 때와 같아야 합니다. 원본의 `NEWFEIcal.txt`는 덮어쓰지 않습니다.

완료된 결과 디렉터리는 재사용하지 않습니다. 재실행 시 `fei_diagnostics_output_2`처럼 새 이름을 지정합니다.
불완전한 export에는 `COMPLETE`가 생성되지 않으며 Python 분석기는 그 결과를 거부합니다.

### 입력 경로 변경이 필요한 경우

기본값은 첨부 원본의 `storage_b1/bsub/Analysis/MizuhashiRD_LS_*_Jpsi` 경로입니다.
다른 경로를 사용하려면 기본 디렉터리만 지정합니다.

```bash
export FEI_DIAG_MC_BASE=/path/to/MizuhashiRD_LS_MC_Jpsi
export FEI_DIAG_DATA_BASE=/path/to/MizuhashiRD_LS_data_Jpsi
root -l -b -q 'THStack_Jpsi_FEI_diagnostics.C("fei_diagnostics_output_2")'
```

각 경로 뒤에 원본과 같은 `SIGNAL_analysis`/`CHG_analysis`/... 및
`validation_v009/final_output_root_after_MVA_Application`이 붙습니다.
전체 경로 구조가 다르면 매크로의 `Jpsi_*_before_dirname_value`를 수정하면 됩니다.

## 먼저 볼 결과

`fei_diagnostics_output/analysis/summary.txt`와 `category_ratios.png`부터 확인합니다.
나중에 결과를 함께 검토하려면 **결과 디렉터리 전체**를 보존하면 됩니다.

| 출력 | 확인하는 질문 |
|---|---|
| `NEWFEIcal_reproduced.txt` | 같은 원본 계산으로 계수를 재현하는가? |
| `analysis/reproduction_check.csv` | 원본 histogram의 수율 및 sumw2가 사건별 export와 일치하는가? |
| `analysis/category_summary.csv` | 보정 전 수율비와 공식 FEI 계수 적용 후 수율비는 각각 얼마인가? |
| `analysis/tag_mode_ratios.csv` | Category 2·3의 어떤 tag mode가 MC-data 차이에 크게 기여하는가? |
| `analysis/signal_mode_ratios.csv` | Xs mode ID 0/1/2에 따라 수율비가 다른가? |
| `analysis/mc_composition.csv` | SIGNAL/CHG/MIX/continuum 중 어느 성분이 큰가? |
| `analysis/truth_by_sample_and_mode.csv` | 모드별로 올바른 태그·오재구성 태그·미확인의 기여는 얼마인가? |
| `analysis/weight_sensitivity.csv` | 특정 보정을 제외하면 수율비가 크게 바뀌는가? |
| `analysis/weight_ranges.csv` | 0, 음수 또는 극단적인 보정값이 있는가? |
| `analysis/cutflow_ratios.csv` | 입력 ntuple에서 시작하는 국소 선택 전후의 차이는 어떤가? |
| `analysis/threshold_scans.csv` | 남아 있는 사건에 BDT/FEI 확률 조건을 더 강하게 적용하면 비율이 바뀌는가? |
| `analysis/event_identifiers.csv` | 같은 event ID가 여러 selected row에 나타나는가? |
| `files.csv` | 입력 파일, sample, 정규화 가중치, tree 정렬 확인 정보 |
| `branch_availability.csv` | 각 파일에 어떤 진단 변수가 존재하는가? |
| `selected_events.csv` | 사건별 원본 가중치, 개별 보정값, 진단 변수 |

Category 0은 모든 neutral tag, Category 1/2/3은 원본의 charged tag 분류입니다.
원본과 마찬가지로 selected row를 세며 사건 중복을 임의로 제거하지 않습니다.

## 해석 범위

- 이 코드는 **원인 후보를 구분할 관측값**을 만듭니다. 큰 차이를 보이는 항목만으로 원인을 확정하지 않습니다.
- 원본 계수는 **관측된 수율 / 전체 MC 수율**입니다. 배경을 제거한 순수 tagging efficiency 측정으로 바꾸지 않습니다.
- 공식 FEI 비교에서는 SIGNAL/CHG/MIX에 원본 환경의 `Corrector_FEI` 계수를 적용하고 continuum은 그대로 둡니다.
  Category 평균은 실제 사건 가중치와 모드 구성을 사용합니다.
- truth를 모르는 데이터에 MC truth 선택을 적용하지 않습니다. 데이터 전체를 correct-tag MC만으로 나눈 가짜 보정계수를 만들지 않습니다.
- `Btag_isSignal`이 없거나 값이 0/1이 아니면 `unknown`입니다. `Bsig_isSignal`도 저장된 변수 자체로 해석해야 합니다.
- 보정 하나를 제외한 결과는 민감도 비교입니다. 이를 새로운 systematic uncertainty로 자동 채택하지 않습니다.
- cutflow는 `input ntuple → signal mode → Bsig_M < 1.0 → tag option`입니다.
  **입력 파일 생성 전에 적용된 skim/선택을 복원하지 못합니다.**
  이 표는 정규화 가중치만 사용합니다. 최종 보정까지 포함한 원본 비율은 `category_summary.csv`입니다.
- threshold scan도 이미 저장된 사건 안에서 조건을 강화하는 비교입니다. 존재하지 않는 이전 단계의 효율을 주장하지 않습니다.
- 일부 파일에서 해당 변수가 없거나 유효하지 않으면 그 Category의 shape/scan 비교를 생략합니다.
  부분적으로 존재하는 MC만 모아 data/MC를 계산하지 않습니다.
- 독립적인 Jpsi 후보 목록인 `Jpsi` tree를 Upsilon/Btag와 행 번호로 연결하지 않습니다.
- 필수 tree/branch 누락, 형식 불일치, tree 길이 불일치, 공통 event 식별자 불일치가 있으면 중단합니다.
  공통 식별자가 없으면 길이만 확인할 수 있으며 `files.csv`에 남습니다.
- 오차는 데이터 계수 통계와 MC `sumw2`만 반영합니다. 원본 재현용 대칭 오차는 delta method입니다.
  그림의 비대칭 오차는 Garwood 데이터 구간과 MC 통계 오차를 근사적으로 합친 것이며 정확한 ratio 신뢰구간이나 유의도가 아닙니다.
  같은 사건을 사용하는 비교 시나리오/scan 점들은 서로 상관되어 있습니다.
  여러 후보가 같은 사건에 속하면 row 기반 오차의 독립성 가정을 별도로 검토해야 합니다.

## 검증과 재생성

Python 수치·integrity 테스트:

```bash
python3 -m pytest test_fei_diagnostics.py -q
```

ROOT가 설치된 환경에서 exporter 자체를 모의 tree로 점검할 수 있습니다.
이 테스트는 분석 전용 calibration header 없이 실행됩니다.

```bash
root -l -b -q 'test_FEIDiagnosticExport.C("fei_diagnostics_smoke")'
python3 analyze_fei_diagnostics.py fei_diagnostics_smoke
```

매크로 재생성:

```bash
python3 build_diagnostic_macro.py
```

빌더는 `source/`의 보존 원본에 명시적인 진단 hook을 넣고 예상 anchor가 다르면 중단합니다.
추출한 입력과 가중치가 원본과 같은지 검토할 수 있도록 자동 생성 파일도 함께 제공합니다.
실제 입력을 실행하기 전까지 계수 재현 및 물리적 원인 분석은 완료되지 않은 상태입니다.

ROOT 입력 확인은 [TTree API](https://root.cern.ch/doc/master/classTTree.html)의
`SetBranchAddress` 반환 상태와 `GetEntry` 동작을 기준으로 구현했습니다.

## 사건별 결과를 이용한 교차 진단

ROOT 실행에서 나온 CSV가 있으면 다시 ROOT를 실행하지 않고 아래처럼 분석할 수 있습니다.
EXPORT_DIRECTORY에는 selected_events.csv, files.csv, reference.csv,
input_cutflow.csv, branch_availability.csv가 있어야 합니다.

```bash
python3 analyze_fei_joint.py EXPORT_DIRECTORY JOINT_OUTPUT_DIRECTORY
python3 -m pytest test_fei_joint.py -q
```

`joint_signal_mode_probability.csv`는 signal-side mode와 FEI probability를 함께 나눕니다.
0.01 이하와 초과를 서로 겹치지 않게 구분하며, 누락된 값은 unknown으로 남깁니다.
`category3_joint.png`/`.pdf`에는 Category 3의 data와 MC tag truth 성분을 표시합니다.
`application_diagnostic.csv`는 전체 data/MC 비율을 correct tag에만 적용하고
나머지에는 official FEI factor를 유지할 때 같은 제어 표본의 수율을 계산합니다.
이는 동일 데이터를 재사용하는 적용 방식 점검이며 독립적인 closure 검증이나
새로운 correct-tag 보정계수의 측정이 아닙니다. 오차나 유의도를 배정하지 않습니다.
다른 카테고리나 잘못 재구성된 tag에 대한 원인도 이 표만으로 확정하지 않습니다.
원본 수율 및 sumw2 재현, 교차 분할의 합 보존, 사용 가능한 사건 식별자의 반복 여부를
확인하고 입력 CSV의 SHA256을 validation.json에 기록합니다.
