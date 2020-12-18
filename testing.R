library(TMB)
TMBtools::export_models()
devtools::load_all()

set.seed(1)

N <- 1000
df <- data.frame(
  age = runif(N, 20, 40),
  female = rbinom(N, size = 1, prob = 0.5),
  biomarker = runif(N, 0, 120)
)

survdf <- simulate_bisurv(
  dist1 = Exponential$new(lambda = 0.1),
  dist2 = Gompertz$new(lambda = 0.5, gamma = 1.5),
  formula1 = ~ age + female,
  formula2 = ~ age + biomarker,
  data = df,
  beta1 = c(age = 0.01, female = 0.5),
  beta2 = c(age = 0.01, biomarker = 0.005),
  maxt1 = 2,
  maxt2 = 4,
  copula = copula::frankCopula(param = 3)
)

testdf <- cbind(df, survdf)

true_beta1 <- c(Y1.ln_lambda = log(0.1), Y1.age = 0.01, Y1.female = 0.5)
true_beta2 <- c(Y2.ln_lambda = log(0.5), Y2.age = 0.01, Y2.biomarker = 0.005)
true_ancillary <- c(theta = 3, Y2.ln_gamma = log(1.5))

# bisurv(
#   formula1 = Surv(eventtime1, status1) ~ age + female,
#   formula2 = Surv(eventtime2, status2) ~ age + female + biomarker,
#   data = testdf
# )

.X1 <- model.matrix(~ age + female, data = testdf)
.beta1 <- rep(0, ncol(.X1))
names(.beta1) <- colnames(.X1)

.X2 <- model.matrix(~ age + biomarker, data = testdf)
.beta2 <- rep(0, ncol(.X2))
names(.beta2) <- colnames(.X2)

f <- TMB::MakeADFun(
  data = list(
    model = "bisurvreg",
    t1 = testdf$eventtime1,
    status1 = testdf$status1,
    t2 = testdf$eventtime2,
    status2 = testdf$status2,
    X1 = .X1,
    X2 = .X2,
    distribution1 = 1,
    distribution2 = 3,
    copula = 2
  ),
  parameters = list(
    beta1 = .beta1,
    beta2 = .beta2,
    ancillary = c(theta = cor(testdf$eventtime1, testdf$eventtime2), log_gamma = 1)
  ),
  silent = TRUE,
  DLL = "bicopula_TMBExports"
)

# f$fn(f$par)
# f$report()
# View(cbind(f$report()))

fit <- nlminb(start = f$par, objective = f$fn, gradient = f$gr, hessian = f$he, control = list(trace = 1))
