set.seed(1)

df <- data.frame(
  age = runif(10, 20, 40),
  female = rbinom(10, size = 1, prob = 0.5),
  egfr = runif(10, 0, 120)
)

simulate_bisurv(
  dist1 = Exponential$new(lambda = 0.1),
  dist2 = Weibull$new(lambda = 0.1, gamma = 1.5),
  formula1 = ~ age + female,
  formula2 = ~ age + egfr,
  data = df,
  beta1 = c(age = 0.01, female = 0.5),
  beta2 = c(age = 0.01, egfr = 0.005),
  maxt1 = 1,
  maxt2 = 3
)
