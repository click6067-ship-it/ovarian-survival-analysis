# ovarian-survival-analysis# Ovarian Cancer Survival Analysis

생존분석 과제용 레포지토리입니다.  
R과 Python(lifelines)을 사용해 `survival` 패키지의 `ovarian` 데이터를 분석했습니다.

## Files

- `notebook/ovarian_survival.ipynb`  
  - Python(lifelines)으로 Kaplan–Meier, log-rank, Cox PH 모델을 수행한 노트북
- `R/ovarian_survival.R`  
  - 동일 분석을 R(survival, survminer)로 수행한 스크립트
- `data/ovarian.csv`  
  - R의 `survival::ovarian` 데이터를 export한 파일 (재현용, 선택)

## Environment

- Python: `pandas`, `matplotlib`, `lifelines`
- R: `survival`, `survminer`

## How to run

1. `data/ovarian.csv` 기준으로:
   - R: `source("R/ovarian_survival.R")`
   - Python: 노트북 `notebook/ovarian_survival.ipynb` 실행
