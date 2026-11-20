Category 2 / 3 원인 추적용 패키지 (2026-11-20)

새 압축을 서버에 풀고, 생성된 fei_category23_diagnosis 폴더에서 실행:

python3 -u trace_fei_stages.py --scope full --details --workers 4 --output fei_category23_trace

numpy/pandas/uproot는 이전 실행에서 사용한 환경을 그대로 사용하면 됩니다.
입력 경로와 24,860개 파일 manifest는 이미 설정돼 있습니다.
기존 BCS 결과나 입력 ROOT를 수정하지 않습니다. 출력 폴더가 이미 있으면 중단합니다.

이번에 구분하려는 내용

1) Category 2
BCS 직전부터 Data/MC가 0.602로 낮음. 어떤 앞선 cut에서 비율이 변하는지 확인.
BCS에서는 경쟁 후보가 있는 사건의 잔류율 차이를 tag 모드, 점수, truth와 연결.
경쟁 없는 사건의 비중이 MC에서 더 높다는 설명은 Category 2에는 맞지 않음.

2) Category 3
BCS 이전에 다른 category 경쟁 후보가 없는 사건이 Data에서 부족함.
이 성분이 어떤 cut, tag 모드 및 MC truth 집단에서 줄어드는지 확인.
후보가 사라지는 과정과 경쟁 후보만 사라지는 과정을 따로 기록.

3) 검증
각 cut을 실제 코드와 대조하고, 저장된 최종 BCS winner를 그대로 추적.
사건 중복 및 입력 tree 정렬, 최종 control 후보의 기존 export 재현 여부 검사.
상세 분해와 사건 연결이 기존 전체 cutflow에 합계 일치하는지도 검사.
COMPLETE이면 검사 통과. REVIEW_REQUIRED 또는 RUNNING만 있으면 결과를 확정하지 말 것.

완료 후 확인 및 전달

fei_category23_trace/summary.txt
fei_category23_trace/category23_diagnostic_summary.txt

위 파일은 진단 안내와 관측된 비율 변화 요약입니다. 물리적 원인을 자동 판정하지 않습니다.
추가 교차 분석을 위해 fei_category23_trace 폴더 전체를 이 컴퓨터로 가져오면 됩니다.
상세 CSV는 압축돼 있으며, 새 ROOT 파일을 만들거나 가져올 필요는 없습니다.

확실히 구분해야 할 한계

이 검사는 저장된 단계 안에서 차이가 생기는 구간과 관련 후보 성분을 찾습니다.
첫 저장 단계부터 차이가 있으면 더 앞선 FEI/재구성/초기 선택 및 생산 설정을 확인해야 합니다.
검출효율이나 붕괴모델 중 어느 것이 물리적 원인인지 이 한 번의 실행으로
반드시 확정된다고 약속할 수는 없습니다. 없는 관측량이나 Data truth를 추정해 채우지 않습니다.
