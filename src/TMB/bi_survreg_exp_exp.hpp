/// @file bi_survreg_exp_exp.hpp

#ifndef bi_survreg_exp_exp_hpp
#define bi_survreg_exp_exp_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// Likelihood for a bivariate copula survival model
template<class Type>
Type bi_survreg_exp_exp(objective_function<Type>* obj) {
  // Survival outcome
  DATA_VECTOR(t1);
  DATA_VECTOR(status1);
  DATA_VECTOR(t2);
  DATA_VECTOR(status2);
  // Data, covariates
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  // Parameters
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER(theta);

  // Number of obs
  int n = t1.size();

  // Marginal survival functions
  vector<Type> u1 = exp(-exp(X1 * beta1) * t1);
  vector<Type> u2 = exp(-exp(X2 * beta2) * t2);
  REPORT(u1);
  REPORT(u2);

  // Copula density
  vector<Type> C = (-1 / theta) * log(1 + (((exp(-theta * u1) - 1) * (exp(-theta * u2) - 1)) / (exp(-theta) - 1)));
  vector<Type> C1 = exp(-theta * u1) * (exp(-theta * u2) - 1) / ((exp(-theta) - 1) + (exp(-theta * u1) - 1) * (exp(-theta * u2) - 1));
  vector<Type> C2 = exp(-theta * u2) * (exp(-theta * u1) - 1) / ((exp(-theta) - 1) + (exp(-theta * u1) - 1) * (exp(-theta * u2) - 1));
  vector<Type> cden = (1 - exp(-theta)) - (1 - exp(-theta * u1)) * (1 - exp(-theta * u2));
  vector<Type> c = (theta * (1 - exp(-theta)) * exp(-theta * (u1 + u2))) / (cden * cden);
  REPORT(C);
  REPORT(C1);
  REPORT(C2);
  REPORT(cden);
  REPORT(c);

  // Likelihood bits
  vector<Type> log_U1 = X1 * beta1 - exp(X1 * beta1) * t1;
  vector<Type> log_U2 = X2 * beta2 - exp(X2 * beta2) * t2;
  vector<Type> log_L = status1 * status2 * (log(c) + log_U1 + log_U2);
  log_L += status1 * (1.0 - status2) * (log_U1 + log(C1));
  log_L += (1.0 - status1) * status2 * (log_U2 + log(C2));
  log_L += (1.0 - status1) * (1.0 - status2) * log(C);
  // The following removes the time unit from the likelihood:
  log_L += status1 * log(t1) + status2 * log(t2);

  REPORT(log_U1);
  REPORT(log_U2);
  REPORT(log_L);

  // // Process infinite values
  // for (int i = 0; i < n; i++) {
  //   if (isFinite(log_L(i))) {
  //     // do nothing
  //   } else {
  //     log_L(i) = 0.0;
  //   }
  // }

  // Likelihood
  Type out = -sum(log_L);

  // Return
  return out;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
