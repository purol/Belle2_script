초기 FEI 선택 추적 (2026-11-20)

목적
Category 2/3의 수율 차이가 Analysis_main 입력에 이미 있는지,
초기 ROE/vertex 선택 중 어느 단계에서 생기는지 구분한다.
특히 Category 3의 tag_mode 30/31/32가 원본의 다른 signal mode에는 있는지도 확인한다.

실행
새 압축을 서버에 풀고 생성된 fei_initial_diagnosis 폴더에서:

python3 -u trace_fei_initial.py --workers 16 --output fei_initial_trace_w16

기존 numpy/pandas/uproot 환경을 그대로 사용한다. 추가 패키지 설치는 필요 없다.
workers를 생략하면 4개다. 원본 ntuple은 이전 중간 파일보다 클 수 있으므로
메모리 사용량과 서버 상태에 따라 줄일 수 있다. 대기 중인 작업은 worker 수로 제한한다.
Ctrl+C 중단은 원본을 수정하지 않지만 이어서 실행하는 기능은 없다.
재실행할 때는 새 출력 폴더 이름을 사용한다. 기존 폴더를 덮어쓰지 않는다.

설정한 원본 입력 위치
MC: /home/belle2/junewoo/storage_ghi/Ntuple/MizuhashiRD_LS_MC_Jpsi/{SAMPLE}/validation
Data: /home/belle2/junewoo/storage_ghi/Ntuple/MizuhashiRD_LS_data_Jpsi/SIGNAL/validation

대조할 저장 결과
MC: /home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_MC_Jpsi/{SAMPLE}_analysis/validation_v009/before_Mbc_cut
Data: /home/belle2/junewoo/storage_b1/bsub/Analysis/MizuhashiRD_LS_data_Jpsi/SIGNAL_analysis/validation_v009/before_Mbc_cut

Jpsi/{SAMPLE} 보조 폴더와 log 폴더를 원본 입력으로 혼합하지 않는다.
validation 폴더의 비재귀 *.root 파일명을 기존 manifest의 24,860개 source stem과
정확히 대조한다. 누락/추가 파일이 있으면 예시 파일명을 출력하고 중단한다.
경로를 추정해 다른 폴더의 파일을 대신 사용하지 않는다.

원본과 before_Mbc의 두 ROOT만 읽는다. 입력은 모두 rb로 열며 원본을 수정하지 않는다.
중간 분석의 7개 snapshot 또는 BCS를 다시 실행할 필요는 없다.
경로가 달라졌다면 --raw-mc-base, --raw-data-base, --mc-base, --data-base로 지정할 수 있다.
파일 대응만 확인할 때는 --preflight-only를 추가한다. 이 경우 결과 폴더도 만들지 않는다.

순서대로 재현하는 조건
1. nROE_Tracks(looseMask) < 0.5
2. ROE pi0 수 < 0.5
3. ROE KS 수 < 0.5
4. Bsig_daughter_0_extraInfo_mychiProb > -0.5
5. Btag_chiProb > -0.5

제공된 Analysis_main.cc 및 Loader.h의 순서, branch 이름, 엄격 부등식을 사용한다.
NaN 비교는 원래 C++처럼 통과하지 못한다. NaN을 0으로 채우지 않는다.
필수 branch가 없으면 중단하며, optional 관측량 부재는 NaN/coverage로 기록한다.

검증
각 파일에서 위 조건을 적용한 후보들의 ID와 핵심 관측량, 초기 cut 관측량을
실제 before_Mbc_cut 후보와 대조한다. 전체 수율과 후보 수도 이전에 완료한
fei_category23_trace_w16의 첫 단계 결과와 일치해야 한다.
Data 사건 ID의 파일 간 중복, tree 간 정렬, 후보 ID 중복도 검사한다.
원래 분석과 연결되지 않는 결과를 조용히 사용하지 않는다.

결과
- summary.txt, run_info.json: 완료 상태, 검증 결과, 코드/입력 reference 해시.
- initial_data_mc_ratios.csv: 초기 선택 단계별 Data/MC 비율.
- initial_stage_summary.csv: 단계별 category/view 사건 수·후보 수·정규화 수율.
- focus_mode_groups.csv: mode 30/31/32, mode 22, mode 18/19의 사건 수.
  all_modes와 control_K_Kpi를 함께 보면 원본 채널 부재와 control 표본 부재를 구분할 수 있다.
- tag_mode_presence.csv: 모든 tag mode의 존재 여부와 단계별 감소.
- tag_signal_mass_composition.csv: tag mode x signal mode x mass<1/그 외/비유한값 집계.
- initial_candidate_trace.csv.gz: 원본 Category 2/3 control 사건 및 원본 Category 3
  mode 30/31/32 사건의 모든 후보. 후보별 최초 탈락 cut과 before_Mbc 존재 여부 포함.
  in_saved_before_mbc는 BCS winner 표시가 아니다.
- initial_event_trace.csv.gz: 각 단계에서 대상 control 후보가 남은 사건의 구성.
- initial_event_losses.csv.gz: 원본 control 사건이 대상 후보를 모두 잃은 첫 단계.
- file_checks.csv, audit_issues.csv, source_files.csv, branch_coverage.csv,
  representative_branch_inventory.csv, raw_nonfinite_counts.csv: 재현 및 입력 상태 확인.

완료되면 fei_initial_trace_w16 폴더 전체를 이 컴퓨터로 가져오면 된다.
COMPLETE가 있어야 검사를 통과한 결과다. REVIEW_REQUIRED이면 먼저 오류를 검토한다.
RUNNING만 있으면 중단되었거나 아직 실행 중인 불완전 결과다.

해석 범위
수율에는 기존 manifest의 표본 정규화만 적용하며 detector/FEI 보정을 새로 넣지 않는다.
같은 사건이 category/view/tag-mode presence 집계에 여러 번 등장할 수 있으므로
이 집계들을 단순 합하지 않는다. MC truth unknown은 incorrect로 바꾸지 않는다.
입력 ntuple은 이미 FEI 및 앞선 skim/reconstruction을 거친 자료다.
입력부터 차이가 있으면 그보다 앞선 생산 설정/재구성도 확인해야 하며,
이 코드가 특정 검출기 효과나 붕괴모델을 물리적 원인으로 자동 확정하지는 않는다.
