/// @file bisurvreg.hpp

// Hazard
template<class Type>
Type bisurvreg_hazard(
    Type Xb, // linear predictor X * beta
    Type t, // observed time
    Type anc, // ancillary parameter for Weibull and Gompertz distributions
    int distribution // 1 = Exponential, 2 = Weibull, 3 = Gompertz
  ) {
  Type h;
  if (distribution == 1) {
    // Exponential
    h = exp(Xb);
  } else if (distribution == 2) {
    // Weibull
    h = (exp(Xb) * pow(t, anc - 1.0)) * anc;
  } else if (distribution == 3) {
    // Gompertz
    h = exp(Xb) * exp(anc * t);
  }
  return h;
}
VECTORIZE4_ttti(bisurvreg_hazard)

// Survival
template<class Type>
Type bisurvreg_surv(
    Type Xb, // linear predictor X * beta
    Type t, // observed time
    Type anc, // ancillary parameter for Weibull and Gompertz distributions
    int distribution // 1 = Exponential, 2 = Weibull, 3 = Gompertz
  ) {
  Type S;
  if (distribution == 1) {
    // Exponential
    S = exp(-exp(Xb) * t);
  } else if (distribution == 2) {
    // Weibull
    S = exp(-exp(Xb) * pow(t, anc));
  } else if (distribution == 3) {
    // Gompertz
    S = exp((-exp(Xb) / anc) * (exp(anc * t) - 1.0));
  }
  return S;
}
VECTORIZE4_ttti(bisurvreg_surv)

// Copula CDF C
template<class Type>
Type bisurvreg_copula_C(
    Type u, // CDF of margin 1
    Type v, // CDF of margin 2
    Type theta, // Parameter 'theta'
    int family // 1 = Independence, 2 = Frank
) {
  Type C;
  if (family == 1) {
    // Independence copula
    C = u * v;
  } else if (family == 2) {
    // Frank
    if (theta != 0) {
      C = (-1.0 / theta) * log(1.0 + (((exp(-theta * u) - 1.0) * (exp(-theta * v) - 1.0)) / (exp(-theta) - 1.0)));
    } else {
      // Special case when theta = 0: independence copula
      C = u * v;
    }
  }
  return C;
}
VECTORIZE4_ttti(bisurvreg_copula_C)

// Copula partial derivative w.r.t. kth outcome: Ck
template<class Type>
Type bisurvreg_copula_Ck(
    Type k, // CDF of margin k
    Type other, // CDF of other margin
    Type theta, // Parameter 'theta'
    int family // 1 = Independence, 2 = Frank
) {
  Type Ck;
  if (family == 1) {
    // Independence copula
    // This *should* be f(k) * other
    // We need to be careful to pass the right f and not F from the main function
    Ck = k * other;
  } else if (family == 2) {
    // Frank
    if (theta != 0) {
      Ck = exp(-theta * k) * (exp(-theta * other) - 1) / ((exp(-theta) - 1) + (exp(-theta * k) - 1) * (exp(-theta * other) - 1));
    } else {
      Ck = k * other;
    }
  }
  return Ck;
}
VECTORIZE4_ttti(bisurvreg_copula_Ck)

// Copula density c
template<class Type>
Type bisurvreg_copula_c(
    Type u, // CDF of margin 1
    Type v, // CDF of margin 2
    Type theta, // Parameter 'theta'
    int family // 1 = Independence, 2 = Frank
) {
  Type c;
  if (family == 1) {
    // Independence copula
    // This *should* be f(u) * f(v)
    c = u * v;
  } else if (family == 2) {
    // Frank
    if (theta != 0) { // cannot have theta here, need to re-code by implementing this check as I()
      Type cden = (1 - exp(-theta)) - (1 - exp(-theta * u)) * (1 - exp(-theta * v));
      c = (theta * (1 - exp(-theta)) * exp(-theta * (u + v))) / (cden * cden);
    } else {
      // Special case when theta = 0: independence copula
      c = u * v;
    }
  }
  return c;
}
VECTORIZE4_ttti(bisurvreg_copula_c)

#ifndef bisurvreg_hpp
#define bisurvreg_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Likelihood for a bivariate copula survival model
template<class Type>
Type bisurvreg(objective_function<Type>* obj) {
  // Survival outcome
  DATA_VECTOR(t1);
  DATA_VECTOR(status1);
  DATA_VECTOR(t2);
  DATA_VECTOR(status2);

  // Covariates
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);

  // Distributions & copula family
  DATA_INTEGER(distribution1);
  DATA_INTEGER(distribution2);
  DATA_INTEGER(family);

  // Parameters
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(ancillary);
  // Process ancillary parameters:
  // 1st is theta
  // 2nd is ancillary of distribution1
  // 3rd is ancillary of distribution3
  Type theta = ancillary(0), anc1 = 0, anc2 = 0;
  if (distribution1 != 1 && distribution2 != 1) {
    anc1 = ancillary(1);
    anc2 = ancillary(2);
  } else if (distribution1 != 1 && distribution2 == 1) {
    anc1 = ancillary(1);
  } else if (distribution1 == 1 && distribution2 != 1) {
    anc2 = ancillary(1);
  }
  // We model the ancillary parameters on the log-scale
  anc1 = exp(anc1);
  anc2 = exp(anc2);

  // Linear predictors
  vector<Type> Xb1 = X1 * beta1;
  vector<Type> Xb2 = X2 * beta2;

  // Number of obs
  int n = t1.size();

  // Marginal survival functions
  vector<Type> u1 = bisurvreg_surv(Xb1, t1, anc1, distribution1);
  vector<Type> u2 = bisurvreg_surv(Xb2, t2, anc2, distribution2);

  // Marginal hazard functions
  vector<Type> h1 = bisurvreg_hazard(Xb1, t1, anc1, distribution1);
  vector<Type> h2 = bisurvreg_hazard(Xb2, t2, anc2, distribution2);

  // Marginal density functions
  vector<Type> f1 = h1 * u1;
  vector<Type> f2 = h2 * u2;

  // Copula components
  vector<Type> C = bisurvreg_copula_C(u1, u2, theta, family);
  vector<Type> C1(n), C2(n), c(n);
  if (theta != 0) {
    C1 = bisurvreg_copula_Ck(u1, u2, theta, family);
    C2 = bisurvreg_copula_Ck(u2, u1, theta, family);
    c = bisurvreg_copula_c(u1, u2, theta, family);
  } else {
    C1 = bisurvreg_copula_Ck(u1, f2, theta, family);
    C2 = bisurvreg_copula_Ck(u2, f1, theta, family);
    c = bisurvreg_copula_c(f1, f2, theta, family);
  }

  // Likelihood bits
  vector<Type> log_U1 = log(h1) + log(u1);
  vector<Type> log_U2 = log(h2) + log(u2);
  vector<Type> log_L = status1 * status2 * (log(c) + log_U1 + log_U2);
  log_L += status1 * (1.0 - status2) * (log_U1 + log(C1));
  log_L += (1.0 - status1) * status2 * (log_U2 + log(C2));
  log_L += (1.0 - status1) * (1.0 - status2) * log(C);
  // The following removes the time unit from the likelihood:
  log_L += status1 * log(t1) + status2 * log(t2);

  REPORT(log_U1);
  REPORT(log_U2);
  REPORT(log_L);

  // Process infinite values
  for (int i = 0; i < n; i++) {
    if (R_finite(asDouble(log_L(i)))) {
      // do nothing
    } else {
      log_L(i) = 0.0;
    }
  }

  // Likelihood
  Type out = -sum(log_L);

  // Return
  return out;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
