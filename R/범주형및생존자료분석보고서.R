########################################################
# Survival analysis with ovarian cancer data
# - 데이터: survival 패키지 내장 ovarian
# - 주요 목표:
#   1) 전체 생존곡선 확인 (K-M)
#   2) 치료군(rx), 잔여종양(resid.ds)에 따른 생존곡선 비교
#   3) Cox 모형으로 치료군/잔여종양/나이의 hazard ratio 추정
#   4) 비례위험가정 검토 (Schoenfeld 잔차)
########################################################

#-------------------------------------------------------
# 0. 패키지 로드 및 데이터 확인
#-------------------------------------------------------

# install.packages("survival")
# install.packages("survminer")

library(survival)
library(survminer)

# survival 패키지 내장 데이터 불러오기
data(ovarian)

# 데이터 구조 확인
str(ovarian)
summary(ovarian)

# 변수 설명 (보고서용 요약 예시)
# - futime : 생존 또는 검열 시점까지의 시간(일)
# - fustat : 사건 발생 여부 (1 = 사망, 0 = 검열)
# - rx     : 치료군 (1 = control, 2 = treatment)
# - resid.ds : 수술 후 잔여 종양 여부 (1 = no residual disease, 2 = yes)
# - age    : 진단 시 나이
# - ecog.ps: ECOG 수행능력 점수 (값이 클수록 상태가 나쁨)


#-------------------------------------------------------
# 1. 전처리: Surv 객체 생성 + 범주형 변수 정리
#-------------------------------------------------------

# Surv 객체 생성

s_obj <- Surv(time = ovarian$futime, event = ovarian$fustat == 1)
head(s_obj)

# 범주형 변수 factor로 변환 (해석과 그래프 라벨을 위해)
ovarian$rx_factor <- factor(
  ovarian$rx,
  levels = c(1, 2),
  labels = c("control", "treatment")
)

ovarian$resid_factor <- factor(
  ovarian$resid.ds,
  levels = c(1, 2),
  labels = c("no_residual", "residual")
)


# ovarian$ecog_factor <- factor(ovarian$ecog.ps)


#-------------------------------------------------------
# 2. 전체 Kaplan–Meier 생존곡선
#-------------------------------------------------------


fit_all <- survfit(Surv(futime, fustat == 1) ~ 1, data = ovarian)

# 요약(각 시간점에서 생존확률, 표본 수 등)
summary(fit_all)

# 전체 생존곡선 시각화
ggsurvplot(
  fit_all,
  data = ovarian,
  risk.table = TRUE,         # 하단에 위험집단 수 표
  conf.int = TRUE,           # 95% 신뢰구간
  surv.median.line = "hv",   # 중앙 생존시간 수평/수직선 표시
  title = "Overall survival curve (ovarian cancer)",
  xlab = "Time (days)",
  ylab = "Survival probability"
)


#-------------------------------------------------------
# 3. 그룹별 K-M 곡선 + log-rank test
#   (1) 치료군(rx), (2) 잔여 종양(resid.ds)
#-------------------------------------------------------

# 3-1) 치료군에 따른 K-M 생존곡선


fit_rx <- survfit(Surv(futime, fustat == 1) ~ rx_factor, data = ovarian)

ggsurvplot(
  fit_rx,
  data = ovarian,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = TRUE,                # log-rank test p-value 출력
  legend.title = "Treatment",
  legend.labs = c("control", "treatment"),
  title = "Survival curves by treatment group",
  xlab = "Time (days)",
  ylab = "Survival probability"
)

# log-rank test 수치로 확인 (survdiff 사용)
survdiff_rx <- survdiff(Surv(futime, fustat == 1) ~ rx_factor, data = ovarian)
survdiff_rx


# 3-2) 잔여 종양 여부에 따른 K-M 생존곡선
fit_resid <- survfit(Surv(futime, fustat == 1) ~ resid_factor, data = ovarian)

ggsurvplot(
  fit_resid,
  data = ovarian,
  risk.table = TRUE,
  conf.int = TRUE,
  pval = TRUE,                # log-rank p-value
  legend.title = "Residual disease",
  legend.labs = c("no residual", "residual"),
  title = "Survival curves by residual disease status",
  xlab = "Time (days)",
  ylab = "Survival probability"
)

# 잔여 종양 여부에 대한 log-rank test
survdiff_resid <- survdiff(Surv(futime, fustat == 1) ~ resid_factor, data = ovarian)
survdiff_resid


#-------------------------------------------------------
# 4. Cox proportional hazards model
#    (치료군 + 잔여 종양 + 나이)
#-------------------------------------------------------


cox_model <- coxph(
  Surv(futime, fustat == 1) ~ rx_factor + resid_factor + age,
  data = ovarian
)

summary(cox_model)

# summary(cox_model)에서 볼 것:
# - coef : 회귀계수 (log hazard ratio)
# - exp(coef) : hazard ratio (HR)
# - p값 : 각 변수의 효과가 유의한지
# - 예: HR > 1 → 해당 변수 수준에서 사망위험 증가
#      HR < 1 → 사망위험 감소


#-------------------------------------------------------
# 5. 비례위험가정 검토: Schoenfeld 잔차
#-------------------------------------------------------


zph_test <- cox.zph(cox_model)
print(zph_test)

# p값 해석:
# - 귀무가설: 해당 공변량의 효과(HR)는 시간에 따라 변하지 않는다.
# - p > 0.05: 비례위험가정이 크게 위반되지 않았다고 볼 수 있음.
# - p < 0.05: 시간이 지남에 따라 HR이 변할 수 있어, 비례위험가정이 의심됨.

#시각화: 각 공변량에 대해 Schoenfeld 잔차를 시간에 따라 그림
# plot(zph_test)
# abline(h = 0, col = "red")

library(survival)
data(ovarian)

write.csv(ovarian, "ovarian.csv", row.names = FALSE)

