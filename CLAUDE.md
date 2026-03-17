# CUBE_HR5 Project Instructions

## 물리 공식 단위 변환 검증 (필수)

물리 공식을 코드 단위로 변환할 때, 구현 후 반드시 다음 3단계를 수행:

1. **차원 분석**: 모든 변수의 지수와 계수가 물리 차원과 일치하는지 확인
2. **기존 코드 cross-check**: 이미 물리 단위로 구현된 코드(예: coolfine1)가 있으면, 새 코드 단위 변환이 동일한 결과를 주는지 역산 검증
3. **수치 sanity check**: 알려진 기준값(nISM, T2_star 등)을 대입하여 예상 물리값과 일치하는지 출력/계산으로 확인

## VPATH 우선순위

소스 파일 수정 시 반드시 VPATH 우선순위를 확인:
`patch/cuda` > `patch/oct_tree` > `patch/Horizon5-master-2` > base dirs

동일 파일이 여러 디렉토리에 있으면 **cuda 버전을 먼저 수정**하고, 안전을 위해 다른 버전도 동기화.
