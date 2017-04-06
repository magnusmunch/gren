#include <R.h>
#include <Rmath.h>
#include <vector>
#include <math.h>
#include <exception>
//#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <random>
#include <chrono>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

const double __PI = 3.141592653589793238462643383279502884197;
const double PISQ = __PI * __PI;
const double HALFPI = 0.5 * __PI;
const double FOURPISQ = 4 * __PI * __PI;
const double __TRUNC = 0.64;
const double __TRUNC_RECIP = 1.0 / __TRUNC;
const double trunc_schedule[] = { // seq(1,4,by=0.01) -> 301 entries.
0.64, 0.68, 0.72, 0.75, 0.78, 0.8, 0.83, 0.85, 0.87, 0.89, 
0.91, 0.93, 0.95, 0.96, 0.98,   1, 1.01, 1.03, 1.04, 1.06, 
1.07, 1.09,  1.1, 1.12, 1.13, 1.15, 1.16, 1.17, 1.19,  1.2,
1.21, 1.23, 1.24, 1.25, 1.26, 1.28, 1.29,  1.3, 1.32, 1.33,
1.34, 1.35, 1.36, 1.38, 1.39,  1.4, 1.41, 1.42, 1.44, 1.45, 
1.46, 1.47, 1.48,  1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.57, 
1.58, 1.59,  1.6, 1.61, 1.62, 1.63, 1.65, 1.66, 1.67, 1.68, 
1.69,  1.7, 1.71, 1.72, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 
1.8 , 1.81, 1.82, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.9, 
1.91, 1.92, 1.93, 1.95, 1.96, 1.97, 1.98, 1.99,    2, 2.01, 
2.02, 2.03, 2.04, 2.05, 2.07, 2.08, 2.09,  2.1, 2.11, 2.12, 
2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.21, 2.22, 2.23, 
2.24, 2.25, 2.26, 2.27, 2.28, 2.29,  2.3, 2.31, 2.32, 2.33, 
2.35, 2.36, 2.37, 2.38, 2.39,  2.4, 2.41, 2.42, 2.43, 2.44, 
2.45, 2.46, 2.47, 2.48, 2.49, 2.51, 2.52, 2.53, 2.54, 2.55, 
2.56, 2.57, 2.58, 2.59,  2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 
2.66, 2.68, 2.69,  2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 
2.77, 2.78, 2.79,  2.8, 2.81, 2.82, 2.83, 2.84, 2.85, 2.87, 
2.88, 2.89,  2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 
2.98, 2.99,    3, 3.01, 3.02, 3.03, 3.04, 3.06, 3.07, 3.08, 
3.09,  3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 
3.19,  3.2, 3.21, 3.22, 3.23, 3.24, 3.25, 3.27, 3.28, 3.29, 
 3.3, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 
 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.49,  3.5, 
3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59,  3.6, 
3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.71, 
3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79,  3.8, 3.81, 
3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89,  3.9, 3.91, 
3.92, 3.93, 3.95, 3.96, 3.97, 3.98, 3.99,    4, 4.01, 4.02, 
4.03, 4.04, 4.05, 4.06, 4.07, 4.08, 4.09,  4.1, 4.11, 4.12, 4.13};
const int grid_size = 81;
const double ygrid[] = {
  0.0625,0.06698584,0.07179365,0.07694653,0.08246924,
  0.08838835,0.09473229,0.1015315,0.1088188,0.1166291,
  0.125,0.1339717,0.1435873,0.1538931,0.1649385,
  0.1767767,0.1894646,0.2030631,0.2176376,0.2332582,
  0.25,0.2679434,0.2871746,0.3077861,0.329877,
  0.3535534,0.3789291,0.4061262,0.4352753,0.4665165,
  0.5,0.5358867,0.5743492,0.6155722,0.659754,
  0.7071068,0.7578583,0.8122524,0.8705506,0.933033,
  1,1.071773,1.148698,1.231144,1.319508,
  1.414214,1.515717,1.624505,1.741101,1.866066,
  2,2.143547,2.297397,2.462289,2.639016,
  2.828427,3.031433,3.24901,3.482202,3.732132,
  4,4.287094,4.594793,4.924578,5.278032,
  5.656854,6.062866,6.498019,6.964405,7.464264,
  8,8.574188,9.189587,9.849155,10.55606,
  11.31371,12.12573,12.99604,13.92881,14.92853,
  16};
const double vgrid[] = {
  -256,-222.8609,-194.0117,-168.897,-147.0334,
  -128,-111.4305,-97.00586,-84.4485,-73.51668,
  -63.99997,-55.71516,-48.50276,-42.22387,-36.75755,
  -31.99844,-27.85472,-24.24634,-21.10349,-18.36524,
  -15.97843,-13.89663,-12.07937,-10.49137,-9.101928,
  -7.884369,-6.815582,-5.875571,-5.047078,-4.315237,
  -3.667256,-3.092143,-2.580459,-2.124095,-1.716085,
  -1.350442,-1.022007,-0.7263359,-0.4595871,-0.2184366,
  0,0.1982309,0.3784427,0.5425468,0.6922181,
  0.828928,0.953973,1.068498,1.173516,1.269928,
  1.358533,1.440046,1.515105,1.584282,1.64809,
  1.706991,1.761401,1.811697,1.858218,1.901274,
  1.941143,1.978081,2.012318,2.044068,2.073521,
  2.100856,2.126234,2.149802,2.171696,2.192042,
  2.210954,2.228537,2.244889,2.260099,2.274249,
  2.287418,2.299673,2.311082,2.321703,2.331593,
  2.340804};
const double IYPI = 3.141592653589793238462643383279502884197;
const double tol    = 1e-8;






//**************************** RNG class ************************************//
class RNG {

 public:

  double norm  (double mean , double sd) {
    return R::rnorm(mean, sd);
  }
  
  double unif  (void) {
    return R::runif(0, 1);
  }
  
  double ltgamma(double shape, double rate, double trunc);
  
  double expon_rate(double rate) {
    return R::rexp(rate);
  }
  
  static double p_gamma_rate(double x, double shape, double rate, int use_log=0) {
    double scale = 1.0 / rate;
    return R::pgamma(x, shape, scale, 1, use_log);
  }
  
  double igauss(double mu, double lambda);
  
  static double Gamma (double x, int use_log=0) {
    double y = R::lgammafn(x);
    if (!use_log) y = exp(y);
    return y;
  }
  
  static double p_norm (double x, int use_log=0) {
    return R::pnorm(x, 0.0, 1.0, 1, use_log);
  } 
  
  static double p_igauss(double x, double mu, double lambda);
  double rtinvchi2(double scale, double trunc);
  double tnorm(double left);
  
  double alphastar(double left) {
    return 0.5 * (left + sqrt(left*left + 4));
  }
  
  double texpon_rate(double left, double rate) {
    return expon_rate(rate) + left;
  }
  
  double gamma_scale (double shape, double scale); // Gamma_Scale
  double gamma_rate(double shape, double rate);
  
}; // RNG

double RNG::ltgamma(double shape, double rate, double trunc) {
  double a = shape;
  double b = rate * trunc;
  if (shape ==1) return expon_rate(1) / rate + trunc;
  double d1 = b-a;
  double d3 = a-1;
  double c0 = 0.5 * (d1 + sqrt(d1*d1 + 4 * b)) / b;
  double x = 0.0;
  bool accept = false;
  while (!accept) {
    x = b + expon_rate(1) / c0;
    double u = unif();
    double l_rho = d3 * log(x) - x * (1-c0);
    double l_M   = d3 * log(d3 / (1-c0)) - d3;
      accept = log(u) <= (l_rho - l_M);
  }
  return trunc * (x/b);
}

double RNG::igauss(double mu, double lambda) {
    // See R code for specifics.
    double mu2 = mu * mu;
    double Y = norm(0.0, 1.0);
    Y *= Y;
    double W = mu + 0.5 * mu2 * Y / lambda;
    double X = W - sqrt(W*W - mu2);
    if (unif() > mu / (mu + X))
      X = mu2 / X;
    return X;
}

double RNG::p_igauss(double x, double mu, double lambda) {
    // z = 1 / mean
    double z = 1 / mu;
    double b = sqrt(lambda / x) * (x * z - 1);
    double a = sqrt(lambda / x) * (x * z + 1) * -1.0;
    double y = RNG::p_norm(b) + exp(2 * lambda * z) * RNG::p_norm(a);
    return y;
}

double RNG::rtinvchi2(double scale, double trunc) {
    double R = trunc / scale;
    // double X = 0.0;
    // // I need to consider using a different truncated normal sampler.
    // double E1 = r.expon_rate(1.0); double E2 = r.expon_rate(1.0);
    // while ( (E1*E1) > (2 * E2 / R)) {
    //   // Rprintf("E %g %g %g %g\n", E1, E2, E1*E1, 2*E2/R);
    //   E1 = r.expon_rate(1.0); E2 = r.expon_rate(1.0);
    // }
    // // Rprintf("E %g %g \n", E1, E2);
    // X = 1 + E1 * R;
    // X = R / (X * X);
    // X = scale * X;
    double E = tnorm(1/sqrt(R));
    double X = scale / (E*E);
    return X;
}

double RNG::tnorm(double left) {
  double rho, ppsl;
  //int count = 1;
  if (left < 0) { // Accept/Reject Normal
    while (true) {
      ppsl = norm(0.0, 1.0);
      if (ppsl > left) return ppsl;
    }
  }
  else { // Accept/Reject Exponential
    // return tnorm_tail(left); // Use Devroye.
    double astar = alphastar(left);
    while (true) {
      ppsl = texpon_rate(left, astar);
      rho  = exp( -0.5 * (ppsl - astar) * (ppsl - astar) );
      if (unif() < rho) return ppsl;
    }
  }
} // tnorm

double RNG::gamma_scale(double shape, double rate) {
  return R::rgamma(shape, 1.0 / rate);
}

double RNG::gamma_rate(double shape, double rate) {
  return R::rgamma(shape, rate);
}


//**************************** Poly-Gamma class******************************//
class PolyaGamma {

  // For sum of Gammas.
  int T;
  std::vector<double> bvec;

  public:

    // Constructors.
    PolyaGamma(int trunc = 200): T(trunc), bvec(T) {
      bvec.resize(T);
      for(int k=0; k < T; ++k){
        // + since we start indexing at 0.
        double d = ((double) k + 0.5);
        bvec[k] = FOURPISQ * d * d;
      }
    } // PolyaGamma

    static double pg_m1(double b, double z);
    static double pg_m2(double b, double z);
    double draw(int n, double z, RNG& r);
    double draw_sum_of_gammas(double n, double z, RNG& r);
    double draw_like_devroye(double z, RNG& r);
    double mass_texpon(double Z);
    double rtigauss(double Z, RNG& r);
    double a(int n, double x);

};

double PolyaGamma::pg_m1(double b, double z) {
  z = fabs(0.5*z);
  double m1 = 0.0;
  if (z > 1e-12)
    m1 = b * tanh(z) / z;
  else
    m1 = b * (1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6));
  return m1* 0.25;
}

double PolyaGamma::pg_m2(double b, double z) {
  z = fabs(0.5*z);
  double m2 = 0.0;
  if (z > 1e-12)
    m2 = (b+1) * b * pow(tanh(z)/z,2) + b * ((tanh(z)-z)/pow(z,3));
  else
    m2 = (b+1) * b * pow(1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6), 2) +
	   b * ((-1.0/3) + (2.0/15) * pow(z,2) - (17.0/315) * pow(z,4));
  return m2* 0.0625;
}

double PolyaGamma::draw(int n, double z, RNG& r) {
  double sum = 0.0;
  for (int i = 0; i < n; ++i)
    sum += draw_like_devroye(z, r);
  return sum;
} // draw

double PolyaGamma::draw_sum_of_gammas(double n, double z, RNG& r) {
  double x = 0;
  double kappa = z * z;
  for(int k=0; k < T; ++k)
    x += r.gamma_scale(n, 1.0) / (bvec[k] + kappa);
  return 2.0 * x;
} // draw_sum_of_gammas

double PolyaGamma::draw_like_devroye(double Z, RNG& r) {
  // Change the parameter.
  Z = fabs(Z) * 0.5;
  // Now sample 0.25 * J^*(1, Z := Z/2).
  double fz = 0.125 * __PI*__PI + 0.5 * Z*Z;
  // ... Problems with large Z?  Try using q_over_p.
  // double p  = 0.5 * __PI * exp(-1.0 * fz * __TRUNC) / fz;
  // double q  = 2 * exp(-1.0 * Z) * pigauss(__TRUNC, Z);
  double X = 0.0;
  double S = 1.0;
  double Y = 0.0;
  // int iter = 0; If you want to keep track of iterations.
  while (true) {
    // if (r.unif() < p/(p+q))
    if ( r.unif() < mass_texpon(Z) )
      X = __TRUNC + r.expon_rate(1) / fz;
    else
      X = rtigauss(Z, r);\
    S = a(0, X);
    Y = r.unif() * S;
    int n = 0;
    bool go = true;
    // Cap the number of iterations?
    while (go) {
      // Break infinite loop.  Put first so it always checks n==0.
      ++n;
      if (n%2==1) {
        S = S - a(n, X);
	      if ( Y<=S ) return 0.25 * X;
      }
      else {
	      S = S + a(n, X);
	      if ( Y>S ) go = false;
      }
    }
    // Need Y <= S in event that Y = S, e.g. when X = 0.
  }
} // draw_like_devroye

double PolyaGamma::mass_texpon(double Z) {
  double t = __TRUNC;
  double fz = 0.125 * __PI*__PI + 0.5 * Z*Z;
  double b = sqrt(1.0 / t) * (t * Z - 1);
  double a = sqrt(1.0 / t) * (t * Z + 1) * -1.0;
  double x0 = log(fz) + fz * t;
  double xb = x0 - Z + RNG::p_norm(b, 1);
  double xa = x0 + Z + RNG::p_norm(a, 1);
  double qdivp = 4 / __PI * ( exp(xb) + exp(xa) );
  return 1.0 / (1.0 + qdivp);
}

double PolyaGamma::rtigauss(double Z, RNG& r) {
  Z = fabs(Z);
  double t = __TRUNC;
  double X = t + 1.0;
  if (__TRUNC_RECIP > Z) { // mu > t
    double alpha = 0.0;
    while (r.unif() > alpha) {
      // X = t + 1.0;
      // while (X > t)
      //   X = 1.0 / r.gamma_rate(0.5, 0.5);
      // Slightly faster to use truncated normal.
      double E1 = r.expon_rate(1.0); 
      double E2 = r.expon_rate(1.0);
      while ( E1*E1 > 2 * E2 / t) {
	      E1 = r.expon_rate(1.0); E2 = r.expon_rate(1.0);
      }
      X = 1 + E1 * t;
      X = t / (X * X);
      alpha = exp(-0.5 * Z*Z * X);
    }
  }
  else {
    double mu = 1.0 / Z;
    while (X > t) {
      double Y = r.norm(0, 1.0); Y *= Y;
      double half_mu = 0.5 * mu;
      double mu_Y    = mu  * Y;
      X = mu + half_mu * mu_Y - half_mu * sqrt(4 * mu_Y + mu_Y * mu_Y);
      if (r.unif() > mu / (mu + X))
	      X = mu*mu / X;
    }
  }
  return X;
}

double PolyaGamma::a(int n, double x) {
  double K = (n + 0.5) * __PI;
  double y = 0;
  if (x > __TRUNC) {
    y = K * exp( -0.5 * K*K * x );
  }
  else if (x > 0) {
    double expnt = -1.5 * (log(0.5 * __PI)  + log(x)) + log(K) - 2.0 * (n+0.5)*(n+0.5) / x;
    y = exp(expnt);
    // y = pow(0.5 * __PI * x, -1.5) * K * exp( -2.0 * (n+0.5)*(n+0.5) / x);
    // ^- unstable for small x?
  }
  return y;
}





//**************************** Poly-GammaALT class***************************//
class PolyaGammaAlt {

  public:

    // Draw.
    double draw(double h, double z, RNG& r, int max_inner=200);
    double draw_abridged(double h, double z, RNG& r, int max_inner=200); 
    double w_left (double trunc, double h, double z);
    double w_right(double trunc, double h, double z);
    double rtigauss(double h, double z, double trunc, RNG& r);
    double a_coef_recursive(double n, double x, double h, double coef_h, double& gamma_nh_over_n);
    double g_tilde(double x, double h, double trunc);
    double pigauss(double x, double Z, double lambda);
    double rtinvchi2(double h, double trunc, RNG& r);

};

double PolyaGammaAlt::draw_abridged(double h, double z, RNG& r, int max_inner) {
  // Change the parameter.
  z = fabs(z) * 0.5;
  int    idx   = (int) floor((h-1.0)*100.0);
  double trunc = trunc_schedule[idx]; // trunc_schedule is a constant
  // Now sample 0.25 * J^*(1, z := z/2).
  double rate_z       = 0.125 * __PI*__PI + 0.5 * z*z;
  double weight_left  = w_left (trunc, h, z);
  double weight_right = w_right(trunc, h, z);
  double prob_right   = weight_right / (weight_right + weight_left);
  // Rprintf("prob_right: %g\n", prob_right);
  double coef1_h = exp(h * log(2.0) - 0.5 * log(2.0 * __PI));
  // double gamma_nh_over_n = RNG::Gamma(h);
  double gnh_over_gn1_gh = 1.0; // Will fill in value on first call to a_coef_recursive.
  int num_trials = 0;
  int total_iter = 0;
  while (num_trials < 10000) {
    num_trials++;
    double X = 0.0;
    double Y = 0.0;
    // if (r.unif() < p/(p+q))
    double uu = r.unif();
    if ( uu < prob_right )
      X = r.ltgamma(h, rate_z, trunc);
    else
      X = rtigauss(h, z, trunc, r);
    // double S  = a_coef(0, X, h);
    double S = a_coef_recursive(0.0, X, h, coef1_h, gnh_over_gn1_gh);
    double a_n = S;
    // double a_n2 = S2;
    // Rprintf("a_n=%g, a_n2=%g\n", a_n, a_n2);
    double gt =  g_tilde(X, h, trunc);
    Y = r.unif() * gt;
    // Rprintf("test gt: %g\n", g_tilde(trunc * 0.1, h, trunc));
    // Rprintf("X, Y, S, gt: %g, %g, %g, %g\n", X, Y, S, gt);
    bool decreasing = false;
    int  n  = 0;
    bool go = true;
    // Cap the number of iterations?
    while (go && n < max_inner) {
      total_iter++;
      ++n;
      double prev = a_n;
      // a_n  = a_coef(n, X, h);
      a_n = a_coef_recursive((double)n, X, h, coef1_h, gnh_over_gn1_gh);
      // Rprintf("a_n=%g, a_n2=%g\n", a_n, a_n2);
      decreasing = a_n <= prev;
      if (n%2==1) {
        S = S - a_n;
        if ( Y<=S && decreasing) return 0.25 * X;
      }
      else {
	      S = S + a_n;
	      if ( Y>S && decreasing) go = false;
      }
    }
    // Need Y <= S in event that Y = S, e.g. when X = 0.
  }
  // We should never get here.
  return -1.0;
} // draw

double PolyaGammaAlt::draw(double h, double z, RNG& r, int max_inner) {
  double n = floor( (h-1.0) / 4.0 );
  double remain = h - 4.0 * n;
  double x = 0.0;
  for (int i = 0; i < (int)n; i++) 
    x += draw_abridged(4.0, z, r);
  if (remain > 4.0)
    x += draw_abridged(0.5 * remain, z, r) + draw_abridged(0.5 * remain, z, r);
  else
    x += draw_abridged(remain, z, r);
  return x;
}

double PolyaGammaAlt::pigauss(double x, double z, double lambda) {
  // z = 1 / mean
  double b = sqrt(lambda / x) * (x * z - 1);
  double a = sqrt(lambda / x) * (x * z + 1) * -1.0;
  double y = RNG::p_norm(b) + exp(2 * lambda * z) * RNG::p_norm(a);
  return y;
}

double PolyaGammaAlt::w_left(double trunc, double h, double z) {
  double out = 0;
  if (z != 0) 
    out = exp(h * (log(2.0) - z)) * pigauss(trunc, z/h, h*h);
  else
    out = exp(h * log(2.0)) * (1.0 - RNG::p_gamma_rate(1/trunc, 0.5, 0.5*h*h));
  return out;
}

double PolyaGammaAlt::w_right(double trunc, double h, double z) {
  double lambda_z = PISQ * 0.125 + 0.5 * z * z;
  double p = exp(h * log(HALFPI / lambda_z)) * (1.0-RNG::p_gamma_rate(trunc, h, lambda_z));
  return p;
}

double PolyaGammaAlt::rtigauss(double h, double z, double trunc, RNG& r) {
  z = fabs(z);
  double mu = h/z;
  double X = trunc + 1.0;
  if (mu > trunc) { // mu > t
    double alpha = 0.0;
    while (r.unif() > alpha) {
      X = rtinvchi2(h, trunc, r);
      alpha = exp(-0.5 * z*z * X);
    }
    // Rprintf("rtigauss, part i: %g\n", X);
  }
  else {
    while (X > trunc) {
      X = r.igauss(mu, h*h);
    }
    // Rprintf("rtigauss, part ii: %g\n", X);
  }
  return X;
}

double PolyaGammaAlt::g_tilde(double x, double h, double trunc) {
  double out = 0;
  if (x > trunc) 
    out = exp(h * log(0.5 * __PI) + (h-1) * log(x) - PISQ * 0.125 * x - RNG::Gamma(h, true));
  else 
    out = h * exp( h * log(2.0) - 0.5 * log(2.0 * __PI * x * x * x) - 0.5 * h * h / x);
    // out = h * pow(2, h) * pow(2 * __PI * pow(x,3), -0.5) * exp(-0.5 * pow(h,2) / x);
  return out;
}

double PolyaGammaAlt::a_coef_recursive(double n, double x, double h, double coef_h, double& gnh_over_gn1_gh) {
  double d_n = 2.0 * (double) n + h;
  // gamma_nh_over_n *= (n + h - 1) / n;  // Can speed up further by separate function for a0 and an, n > 0.
  if (n != 0)
    gnh_over_gn1_gh *= (n + h - 1) / n;
  else
    gnh_over_gn1_gh  = 1.0;
  double coef       = coef_h * gnh_over_gn1_gh;
  double log_kernel = - 0.5 * (log(x * x * x) + d_n * d_n / x) + log(d_n);
  return coef * exp(log_kernel);
  // double out = exp(out) is a legal command.  Weird.
}

double PolyaGammaAlt::rtinvchi2(double h, double trunc, RNG& r) {
  double h2 = h * h;
  double R = trunc / h2;
  double X = 0.0;
  // I need to consider using a different truncated normal sampler.
  double E1 = r.expon_rate(1.0); double E2 = r.expon_rate(1.0);
  while ( (E1*E1) > (2 * E2 / R)) {
    // Rprintf("E %g %g %g %g\n", E1, E2, E1*E1, 2*E2/R);
    E1 = r.expon_rate(1.0); E2 = r.expon_rate(1.0);
  }
  // Rprintf("E %g %g \n", E1, E2);
  X = 1 + E1 * R;
  X = R / (X * X);
  X = h2 * X;
  return X;
}







//**************************** PolyaGammaSP ********************************//
struct FD {
  double val;
  double der;
};

struct Line {
  double slope;
  double icept;
};

class PolyaGammaSP {

  public:

    int draw(double& d, double h, double z, RNG& r, int max_iter=200);

  protected:
  
    void   delta_func(double x, double z  , FD& delta);
    double y_func(double v); // y = tan(sqrt(v)) / sqrt(v);
    double tangent_to_eta(double x, double z, double mid, Line& tl);
    double rtigauss(double mu, double lambda, double trunc, RNG& r);
    double sp_approx(double x, double n, double z);
    double phi_func  (double x, double mid, FD& phi);
    double cos_rt(double v);
    double v_eval(double y, double tol=1e-9, int max_iter=1000);
    void   fdf_eval(double v, void* params, double* fp, double* dfp);
    void   ydy_eval(double v, double* yp, double* dyp);
    double   y_eval(double v);

};

int PolyaGammaSP::draw(double& d, double n, double z, RNG& r, int maxiter) {
  z = 0.5 * fabs(z);
  double xl = y_func(-1*z*z);    // Mode of phi - Left point.
  double md = xl * 1.1;          // Mid point.
  double xr = xl * 1.2;          // Right point.
  // Rprintf("xl, md, xr: %g, %g, %g\n", xl, md, xr);
  // Inflation constants
  // double vmd  = yv.v_func(md);
  double vmd  = v_eval(md);
  double K2md = 0.0;
  if (fabs(vmd) >= 1e-6) 
    K2md = md*md + (1-md) / vmd;
  else
    K2md = md*md - 1/3 - (2/15) * vmd;
  double m2 = md * md;
  double al = m2*md / K2md;
  double ar = m2    / K2md;
  // Rprintf("vmd, K2md, al, ar: %g, %g %g %g\n", vmd, K2md, al, ar);
  // Tangent lines info.
  Line ll, lr;
  tangent_to_eta(xl, z, md, ll);
  tangent_to_eta(xr, z, md, lr);
  double rl = -1. * ll.slope;
  double rr = -1. * lr.slope;
  double il = ll.icept;
  double ir = lr.icept;
  // Rprintf("rl, rr, il, ir: %g, %g, %g, %g\n", rl, rr, il, ir);
  // Constants
  double lcn = 0.5 * log(0.5 * n / __PI);
  double rt2rl = sqrt(2 * rl);
  // Rprintf("sqrt(rl): %g\n", rt2rl);
  // Weights
  double wl, wr, wt, pl;
  wl = exp(0.5 * log(al) - n * rt2rl + n * il + 0.5 * n * 1./md) * 
    RNG::p_igauss(md, 1./rt2rl, n);
  wr = exp(0.5 * log(ar) + lcn - n * log(n * rr) + n * ir - n * log(md)) *
    // yv.upperIncompleteGamma(md, n, n*rr);
    RNG::Gamma(n) * (1.0 - RNG::p_gamma_rate(md, n, n*rr));
  // Rprintf("wl, wr: %g, %g\n", wl, wr);
  wt = wl + wr;
  pl = wl / wt;
  // Sample
  bool go  = true;
  int iter = 0;
  double X = 2.0;
  double F = 0.0;
  while(go && iter < maxiter) {
    // Put first so check on first pass.
    iter++;
    double phi_ev;
    if (r.unif() < pl) {
      X = rtigauss(1./rt2rl, n, md, r);
      phi_ev = n * (il - rl * X) + 0.5 * n * ((1.-1./X) - (1.-1./md));
      F = exp(0.5 * log(al) + lcn - 1.5 * log(X) + phi_ev);
    }
    else {
      X = r.ltgamma(n, n*rr, md);
      phi_ev = n * (ir - rr * X) + n * (log(X) - log(md));
      F = exp(0.5 * log(ar) + lcn + phi_ev) / X;
    }
    double spa = sp_approx(X, n, z);
    if (F * r.unif() < spa) go = false;
  }
  // return n * 0.25 * X;
  d = n * 0.25 * X;
  return iter;
}

double PolyaGammaSP::y_func(double v) {
  double tol = 1e-6;
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v > tol)
    y = tan(r) / r;
  else if (v < -1*tol)
    y = tanh(r) / r;
  else
    y = 1 + (1/3) * v + (2/15) * v * v + (17/315) * v * v * v;
  return y;
}

double PolyaGammaSP::tangent_to_eta(double x, double z, double mid, Line& tl) {
  FD phi, delta, eta;
  double v;
  v = phi_func(x, z, phi);
  delta_func(x, mid, delta);
  eta.val = phi.val - delta.val;
  eta.der = phi.der - delta.der;
  // Rprintf("v=%g\nphi=%g, phi.d=%g\ndelta=%g, delta.d=%g\neta=%g, eta.d=%g\n",
  //    v, phi.val, phi.der, delta.val, delta.der, eta.val, eta.der);
  tl.slope = eta.der;
  tl.icept = eta.val - eta.der * x;
  return v;
}

double PolyaGammaSP::rtigauss(double mu, double lambda, double trunc, RNG& r) {
  // mu = fabs(mu);
  double X = trunc + 1.0;
  if (trunc < mu) { // mu > t
    double alpha = 0.0;
    while (r.unif() > alpha) {
      X = r.rtinvchi2(lambda, trunc);
      alpha = exp(-0.5 * lambda / (mu*mu) * X);
    }
    // Rprintf("rtigauss, part i: %g\n", X);
  }
  else {
    while (X > trunc) {
      X = r.igauss(mu, lambda);
    }
    // Rprintf("rtigauss, part ii: %g\n", X);
  }
  return X;
}

double PolyaGammaSP::sp_approx(double x, double n, double z) {
  // double v  = yv.v_func(x);
  double v = v_eval(x);
  double u  = 0.5 * v;
  double z2 = z * z;
  double t  = u + 0.5 * z2;
  // double m  = y_func(-1 * z2);
  double phi = log(cosh(z)) - log(cos_rt(v)) - t * x;
  double K2  = 0.0;
  if (fabs(v) >= 1e-6) 
    K2 = x*x + (1-x) / v;
  else
    K2 = x*x - 1/3 - (2/15) * v;
  double log_spa = 0.5 * log(0.5 * n / __PI) - 0.5 * log(K2) + n * phi;
  return exp(log_spa);
}

double PolyaGammaSP::phi_func(double x, double z, FD& phi) {
  // double v = yv.v_func(x);
  double v = v_eval(x);
  double u = 0.5 * v;
  double t = u + 0.5 * z*z;
  phi.val = log(cosh(fabs(z))) - log(cos_rt(v)) - t * x;
  phi.der = -1.0 * t;
  return v;
}

void PolyaGammaSP::delta_func(double x, double mid, FD& delta) {
  if (x >= mid) {
    delta.val = log(x) - log(mid);
    delta.der = 1.0 / x;
  }
  else {
    delta.val = 0.5 * (1 - 1.0 / x) - 0.5 * (1 - 1.0 / mid);
    delta.der = 0.5 / (x*x);
  }
}

double PolyaGammaSP::cos_rt(double v) {
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v >= 0)
    y = cos(r);
  else
    y = cosh(r);
  return y;
}

double PolyaGammaSP::v_eval(double y, double tol, int max_iter) {
  double ylower = ygrid[0];
  double yupper = ygrid[grid_size-1];
  if (y < ylower) {
    return -1. / (y*y);
  } else if (y > yupper) {
    double v = atan(0.5 * y * IYPI);
    return v*v;
  }
  else if (y==1) return 0.0;
  double id = (log(y) / log(2.0) + 4.0) / 0.1;
  // Rprintf("y, id, y[id], v[id]: %g, %g, %g, %g\n", y, id, ygrid[(int)id], vgrid[(int)id]);
  // C++ default is truncate decimal portion.
  int idlow  = (int)id;
  int idhigh = (int)id + 1;
  double vl  = vgrid[idlow];  // lower bound
  double vh  = vgrid[idhigh]; // upper bound
  int    iter = 0;
  double diff = tol + 1.0;
  double vnew = vl;
  double vold = vl;
  double f0, f1;
  while (diff > tol && iter < max_iter) {
    iter++;
    vold = vnew;
    fdf_eval(vold, &y, &f0, &f1);
    vnew = vold - f0 / f1;
    vnew = vnew > vh ? vh : vnew;
    vnew = vnew < vl ? vl : vnew;
    diff = fabs(vnew - vold);
    // Rprintf("iter: %i, v: %g, diff: %g\n", iter, vnew, diff);
  }
  return vnew;
}

void PolyaGammaSP::fdf_eval(double v, void* params, double* fp, double* dfp) {
  double y = *((double*) params);
  ydy_eval(v, fp, dfp);
  *fp  -= y;
}

void PolyaGammaSP::ydy_eval(double v, double* yp, double* dyp) {
  // double r   = sqrt(fabs(v));
  double y = y_eval(v);
  *yp = y;
  if (fabs(v) >= tol) 
    *dyp = 0.5 * (y*y + (1-y) / v);
  else
    *dyp = 0.5 * (y*y - 1/3 - (2/15) * v);
}

double PolyaGammaSP::y_eval(double v) {
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v > tol)
    y = tan(r) / r;
  else if (v < -1*tol)
    y = tanh(r) / r;
  else
    y = 1 + (1/3) * v + (2/15) * v * v + (17/315) * v * v * v;
  return y;
}


vec rpg_hybrid(std::vector<double> & h, std::vector<double> & z, int num) {
  
  std::vector<double> x(num,0.0);
  RNG r;
  PolyaGamma dv;
  PolyaGammaAlt alt;
  PolyaGammaSP sp;
  
  for(int i=0; i < num; ++i){
    double b = h[i];
    if (b > 170) {
      double m = dv.pg_m1(b,z[i]);
      double v = dv.pg_m2(b,z[i]) - m*m;
      x[i] = r.norm(m, sqrt(v));
    }
    else if (b > 13) {
      sp.draw(x[i], b, z[i], r);
    }
    else if (b==1 || b==2) {
      x[i] = dv.draw((int)b, z[i], r);
    }
    else if (b > 1) {
      x[i] = alt.draw(b, z[i], r);
    }
    else if (b > 0) {
      x[i] = dv.draw_sum_of_gammas(b, z[i], r);
    }
    else {
      x[i] = 0.0;
    }
  }
  
  vec out = conv_to< colvec >::from(x);
  return out;
  
}


//**************************** parameter sample functions *****************************//
// class rpar {
//   
//   public:
//   
//     vec romegaC(vec beta, mat x, vec m, int n, int p, bool intercept);
//     vec rtauC(vec beta0, int p, double lambda1, vec lambda2, bool intercept);
//     vec rbetaC(mat x, vec kappa, int n, int p, vec tau, vec om, vec lambda2, bool intercept);
//   
// }; // rpar

// [[Rcpp::export]]
vec romegaC(vec beta, mat x, vec m, int n, int p, bool intercept) {
  
  vec evec(n,fill::ones);
  if(intercept) {
    x = join_rows(evec,x);
  }
  vec z = x*beta;
  
  std::vector<double> m2 = conv_to< std::vector<double> >::from(m);
  std::vector<double> z2 = conv_to< std::vector<double> >::from(z);
  vec out = rpg_hybrid(m2, z2, n);
  return out;
}

/*vec romegaC(vec beta, mat x, int n, int p, bool intercept) {
  
  std::vector<double> z(n, 0.0);
  std::vector<double> h(n, 1.0);
  for(int i=0; i<n; i++) {
    for(int j=0; j<p; j++) {
      z[i] += x(i,j)*beta(j);
    }
  }  
  vec out = rpg_hybrid(h, z, n);
  return out;
}*/

// [[Rcpp::export]]
//double rtauC(vec beta0, int p, double lambda1, vec lambda2, bool intercept) {
vec rtauC(vec beta0, int p, double lambda1, vec lambda2, bool intercept) {
// mat rtauC(vec beta0, int p, double lambda1, vec lambda2, bool intercept) {
  
  //std::random_device rd;
  std::mt19937 gen{static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count())};
  std::normal_distribution<double> dn(0.0,1.0);
  std::uniform_real_distribution<double> du(0.0,1.0);
  
  vec beta(p), tau(p);
  if(intercept) {
    beta = beta0.rows(1,p);
  } else {
    beta = beta0;
  }
	vec psi = pow(lambda1,2.0)/(4.0*lambda2);
	vec chi = lambda2 % square(beta);
  	
	double Y, Ysq, U, z, prob;
	
	for(int i=0; i<p; i++) {
	  bool isfin = false;
		while(isfin==false) {
		  Y = dn(gen);
	    Ysq = pow(Y,2.0);
	    U = du(gen);
		  z = sqrt(psi(i)/chi(i)) + Ysq/(2.0*chi(i)) - sqrt(sqrt(psi(i))*Ysq/pow(chi(i),1.5) + pow(Ysq,2.0)/(4.0*pow(chi(i),2.0)));
		  prob = 1.0/(1.0 + z*sqrt(chi(i)/psi(i)));
		  if(U <= prob)
			  tau(i) = 1.0/z + 1.0;
		  else
			  tau(i) = chi(i)*z/psi(i) + 1.0;

		  if(tau(i) <= 1.0)
		    tau(i) = 1.001;

		  isfin = (z!=0.0);
		}
	}
	
	//Y = dn(gen);

	// mat test(p,7);
	// test.submat(0,0,p-1,0) = psi;
	// test.submat(0,1,p-1,1) = chi;
	// test.submat(0,2,p-1,2) = Y;
	// test.submat(0,3,p-1,3) = U;
	// test.submat(0,4,p-1,4) = z;
	// test.submat(0,5,p-1,5) = prob;
	// test.submat(0,6,p-1,6) = tau;
	// return test;
	//return tau;
	return tau;

}

// // [[Rcpp::export]]
// vec rtauC(vec beta0, int p, double lambda1, vec lambda2, bool intercept) {
//   // mat rtauC(vec beta0, int p, double lambda1, vec lambda2, bool intercept) {
//   
//   vec beta(p);
//   if(intercept) {
//     beta = beta0.rows(1,p);
//   } else {
//     beta = beta0;
//   }
//   vec psi = pow(lambda1,2.0)/(4.0*lambda2);
//   vec chi = lambda2 % square(beta);
//   vec tau(p);
//   
//   vec Y(p,fill::randn);
//   vec Ysq = square(Y);
//   vec U(p,fill::randu);
//   
//   vec z = sqrt(psi/chi) + Ysq/(2.0*chi) - sqrt(sqrt(psi)%Ysq/(chi%sqrt(chi)) + square(Ysq)/(4.0*square(chi)));	
//   vec prob = 1.0/(1.0 + z%sqrt(chi/psi));
//   
//   for(int i=0; i<p; i++) {
//     if(U(i) <= prob(i))
//       tau(i) = 1.0/z(i) + 1.0;
//     else
//       tau(i) = chi(i)*z(i)/psi(i) + 1.0;
//     
//     if(tau(i) <= 1.0)
//       tau(i) = 1.001;
//     
//   }
//   
//   // mat test(p,7);
//   // test.submat(0,0,p-1,0) = psi;
//   // test.submat(0,1,p-1,1) = chi;
//   // test.submat(0,2,p-1,2) = Y;
//   // test.submat(0,3,p-1,3) = U;
//   // test.submat(0,4,p-1,4) = z;
//   // test.submat(0,5,p-1,5) = prob;
//   // test.submat(0,6,p-1,6) = tau;
//   // return test;
//   return tau;
//   
// }

// [[Rcpp::export]]
vec rbetaC(mat x, vec kappa, int n, int p, vec tau, vec om, vec lambda2, bool intercept) {
   
  vec evec(n,fill::ones);
  mat xaug = join_rows(evec,x);
  mat trx = x.t();
  vec hinv = 1/lambda2 - 1/(lambda2%tau);
  vec ominv = 1/om;
  int p0 = p;
  if(intercept) {
    p0 += 1;
  }
  
  double invtrom;
  mat Omadj(n,n);
  mat txhinv = trx.each_col() % hinv;
  mat xtrxom(n,n);
  mat Ainv(p,p);
  vec ainvxom(p), sigmapart(p);
  mat sigma(p0,p0);
  vec mu(p0);
	
	if (intercept) {
	  invtrom = 1/sum(om);
	  Omadj = diagmat(om) - invtrom*(om*om.t());
	  xtrxom = x*txhinv*Omadj;
	  xtrxom.diag() += 1;
	  Ainv= diagmat(hinv) - txhinv*Omadj*xtrxom.i()*txhinv.t();
	  ainvxom = Ainv*trx*om;
	  sigmapart = -invtrom*ainvxom;
	  sigma(0,0) = as_scalar(invtrom + pow(invtrom,2.0)*(om.t()*x*ainvxom));
	  sigma.submat(1,0,p,0) = sigmapart;
	  sigma.submat(0,1,0,p) = sigmapart.t();
	  sigma.submat(1,1,p,p) = Ainv;
	  mu = sigma*xaug.t()*kappa;
	} else {
	  xtrxom = x*txhinv + diagmat(ominv);
	  sigma = diagmat(hinv) - txhinv*xtrxom.i()*txhinv.t();
	  mu = sigma*trx*kappa;
	}
	
	vec beta0 = randn(p0);
	vec beta = chol(sigma)*beta0 + mu;
	return beta;

}



// [[Rcpp::export]]
List gibbsC(mat x, vec y, vec m, int n, int p, double lambda1, vec lambda2, vec b0, bool intercept, int K) {
  
  int p0 = p;
  if(intercept) {
    p0 += 1;
  }
  mat seq_beta(p0, K);
  mat seq_omega(n, K);
  mat seq_tau(p, K);
  
  vec kappa = y - 0.5*m;
  
  mat beta = b0;
  vec omega(n);
  vec tau(p);  
  // vec tau = b0;
  // mat test(p,7);
  // 
  for(int k=0; k<K; k++) {
    omega = romegaC(beta, x, m, n, p, intercept);
    tau = rtauC(beta, p, lambda1, lambda2, intercept);
    beta = rbetaC(x, kappa, n, p, tau, omega, lambda2, intercept);

    seq_tau.submat(0,k,p-1,k) = tau;
    seq_beta.submat(0,k,p0-1,k) = beta;
    seq_omega.submat(0,k,n-1,k) = omega;
    
    
  }
  
  return List::create(Named("beta") = seq_beta,
                      Named("tau") = seq_tau,
  // return List::create(Named("tau") = test,
                      Named("omega") = seq_omega);
}