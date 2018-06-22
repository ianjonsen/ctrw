#include <TMB.hpp>

using namespace density;

template<class Type>
  Type objective_function<Type>::operator() ()
  {
    // DATA
    DATA_MATRIX(Y);	            //  (x, y) observations
    DATA_VECTOR(dt);            //  time diff in some appropriate unit. this should contain dt for both interp and obs positions.
    DATA_IVECTOR(isd);          //  indexes observations vs. interpolation points

    DATA_INTEGER(obs_mod);       //  indicates which obs error model to be used

    // for KF observation model
    DATA_VECTOR(m);             //  m is the semi-minor axis length
    DATA_VECTOR(M);             //  M is the semi-major axis length
    DATA_VECTOR(c);             //  c is the orientation of the error ellipse

    // for LS observation model
    DATA_MATRIX(K);                 // error weighting factors for LS obs model

    // PARAMETERS
    PARAMETER_VECTOR(l_sigma);  //  Innovation variance (link scale)
    PARAMETER(l_rho_p);         //  Innovation correlation (link scale)
    PARAMETER_MATRIX(X);        //  Predicted locations TP - length(X) should be same as length(dt) - i.e. both interp & obs pos.

    // for LS observation model
    PARAMETER_VECTOR(l_tau);        // error dispersion for LS obs model (log scale)
    PARAMETER(l_rho_o);             // error correlation


    // Tronsform parameters
    vector<Type> sigma = exp(l_sigma);
    Type rho_p = Type(2.0) / (Type(1.0) + exp(-l_rho_p)) - Type(1.0);
    vector<Type> tau = exp(l_tau);
    Type rho_o = Type(2.0) / (Type(1.0) + exp(-l_rho_o)) - Type(1.0);

    // 2 x 2 covariance matrix for innovations
    matrix<Type> cov(2, 2);
    matrix<Type> cov_dt(2, 2);            // tmp matrix for dt * cov calcs withn process loop
    // 2 x 2 covariance matrix for observations
    matrix<Type> cov_obs(2, 2);

    cov(0, 0) = sigma(0) * sigma(0);
    cov(0, 1) = rho_p * sigma(0) * sigma(1);
    cov(1, 0) = cov(0, 1);
    cov(1, 1) = sigma(1) * sigma(1);

    parallel_accumulator<Type> jnll(this);       // (Complete data) negative log likelihood
    MVNORM_t<Type> nll_proc(cov);	              // Multivariate Normal for process
    MVNORM_t<Type> nll_obs;                     // Multivariate Normal for observations

    // process model
    for(int i = 1; i < X.rows(); ++i) {
      cov_dt = pow(dt(i), 2) * cov;
      nll_proc.setSigma(cov_dt);
      jnll += nll_proc(X.row(i) - X.row(i - 1));
    }

    // observation model
    for(int i=0; i < Y.rows(); ++i) {
      if(isd(i) == 1) {
        if(obs_mod == 0) {
          // Argos Least Squares observations
          Type s = tau(0) * K(i,0);
          Type q = tau(1) * K(i,1);
          cov_obs(0,0) = pow(s, 2);
          cov_obs(1,1) = pow(q, 2);
          cov_obs(0,1) = s * q * rho_o;
          cov_obs(1,0) = cov_obs(0,1);
        } else {
          // Argos Kalman Filter observations
          Type s2c = sin(c(i)) * sin(c(i));  // sin2(c)
          Type c2c = cos(c(i)) * cos(c(i));  // cos2(c)
          Type M2  = (M(i)/sqrt(2)) * (M(i)/sqrt(2));
          Type m2 = (m(i)/sqrt(2)) * (m(i)/sqrt(2));

          cov_obs(0,0) = M2 * s2c + m2 * c2c;
          cov_obs(1,1) = M2 * c2c + m2 * s2c;
          cov_obs(0,1) = (0.5 * (pow(M(i),2) - pow(m(i),2))) * cos(c(i)) * sin(c(i));
          cov_obs(1,0) = cov_obs(0,1);
        }
        nll_obs.setSigma(cov_obs);   // set up i-th obs cov matrix

        jnll += nll_obs(Y.row(i) - X.row(i));   // innovations
      }
    }

    ADREPORT(rho_p);
    ADREPORT(sigma);
    ADREPORT(rho_o);
    ADREPORT(tau);

    return jnll;
  }
