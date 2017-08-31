/******************************* Preamble **************************
* code belonging to ENVB2.pdf					   *						
* version: 01							   *
* author: Magnus M?nch						   *
* created: 25-11-2016						   *
* last edited: 29-11-2016					   *
*******************************************************************/

/******************************* Notes *****************************
* 17-03-2017: works fine                                           *
* 25-11-2016: gives negative variances                             *
*******************************************************************/

#include <math.h>
#include <RcppArmadillo.h>
// #include <ctime>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List est_param(mat x, vec kappa, vec m, int n, int p, vec ciold, double phi, vec chiold, vec lambda2, bool intercept) {
    
    vec evec(n,fill::ones);
    mat xaug = join_rows(evec,x);
    mat trx = x.t();
    vec h = lambda2 + lambda2%sqrt(phi/chiold);  
    vec om = 0.5*m/ciold%tanh(ciold/2);
    vec ominv = 1/om;
    vec hinv = 1/h;
    int p0 = p;
    if(intercept) {
      p0 += p0;
    }

    double invtrom;
    mat Omadj(n,n);
    mat txhinv = trx.each_col() % hinv;
    mat xtrxom(n,n);
    mat Ainv(p,p);
    vec ainvxom(p);
    mat sigma(p + 1,p + 1);
    vec mu(p + 1);
    vec varmean(p0);
    
    // mat m1(p,n),m2(p,n), m3(n,p+1);
    // clock_t begin1, begin2;
    // clock_t end1, end2;
    
    if (intercept) {
      invtrom = 1/sum(om);
      Omadj = diagmat(om) - invtrom*(om*om.t());
      xtrxom = x*txhinv*Omadj;
      xtrxom.diag() += 1;
      Ainv= diagmat(hinv) - txhinv*Omadj*xtrxom.i()*txhinv.t(); // the slow part is the last multiplication
      
      // m1 = txhinv*Omadj*xtrxom.i();
      // 
      // // this is a slow part
      // begin1 = clock();
      // m2 = m1*txhinv.t();
      // end1 = clock();
      // 
      // Ainv = diagmat(hinv) - m2;

      ainvxom = Ainv*trx*om;
      sigma(0,0) = as_scalar(invtrom + pow(invtrom,2.0)*(om.t()*x*ainvxom));
      sigma.submat(1,0,p,0) = -invtrom*ainvxom;
      sigma.submat(0,1,0,p) = -invtrom*ainvxom.t();
      sigma.submat(1,1,p,p) = Ainv;
      mu = sigma*xaug.t()*kappa;
      
      // // this is a slow part
      // begin2 = clock();
      // m3 = xaug*sigma;
      // end2 = clock();
      // 
      // ciold = sqrt(diagvec(m3*xaug.t()) + square(xaug*mu));
      
      ciold = sqrt(diagvec(xaug*sigma*xaug.t()) + square(xaug*mu)); // the slow part is the xaug*sigma part
      
      varmean = diagvec(sigma) + square(mu);
      chiold = lambda2%varmean.rows(1,p);
    } else {
      xtrxom = x * txhinv + diagmat(ominv);
      sigma = diagmat(hinv) - txhinv*xtrxom.i()*txhinv.t();
      mu = sigma*trx*kappa;
      ciold = sqrt(diagvec(x*sigma*trx) + square(x*mu));
      varmean = diagvec(sigma) + square(mu);
      chiold = lambda2%varmean;
    }
    
    // double elapsed_secs1 = double(end1 - begin1) / CLOCKS_PER_SEC;
    // double elapsed_secs2 = double(end2 - begin2) / CLOCKS_PER_SEC;

    return List::create(Named("sigma") = sigma,
                        Named("mu") = mu,
                        Named("ci") = ciold,
                        Named("chi") = chiold);
                        // Named("time1") = elapsed_secs1,
                        // Named("time2") = elapsed_secs2);
}

// [[Rcpp::export]]
List est_param2(mat x, vec kappa, vec m, int n, int p, vec ciold, double phi, vec chiold, vec lambda2, bool intercept) {
  
  // initialize c, chi, dsigma and mu
  vec ci(n);
  vec chi(p);
  vec dsigma(p);
  vec mu(p);
  
  // used both with and without intercept
  mat xt = x.t();
  vec hinv = 1/(lambda2 + lambda2%sqrt(phi/chiold)); 
  vec w = 0.5*m/ciold%tanh(ciold/2);
  double trominv;
  double trominvsq;
  double sumkap;
  mat Om = diagmat(w);
  if(intercept) {
    trominv = 1/sum(w);
    trominvsq = pow(trominv, 2.0);
    Om = Om - trominv*(w*w.t());
    sumkap = sum(kappa);
  }
  mat C = x.each_row() % hinv.t();
  mat Ct = C.t();
  mat A = C*xt;
  mat B = A*Om;
  B.diag() += 1;
  B = Om*B.i();
  mat AB = A*B;
  vec Ak = A*kappa;
  mat CtB = Ct*B;
  vec ABAk = AB*Ak;
  mat ABC = AB*C;
  vec Ctk = Ct*kappa;
  vec CtBAk = CtB*Ak;
  mat had1 = C - ABC;
  
  // only used when intercept included
  rowvec wtA(n);
  mat ABA(n,n);
  vec ABAw(n);
  rowvec wtC(p);
  
  if(intercept) {
    wtA = w.t()*A;
    ABA = AB*A;
    ABAw = AB*wtA.t();
    wtC = w.t()*C;
    had1 = (had1.each_row() + (trominv*w.t()*ABC - trominv*wtC))%x;
    double val1 = conv_to<double>::from(trominv + trominvsq*wtA*w - trominvsq*w.t()*ABAw);
    double val2 = conv_to<double>::from(sumkap*trominv + sumkap*trominvsq*wtA*w - 
                                        sumkap*trominvsq*w.t()*ABAw - 
                                        trominv*wtA*kappa + trominv*w.t()*ABAk);
    ci = sqrt(val1 - trominv*wtA.t() + trominv*ABAw + sum(had1,1) + 
      square(A*kappa - ABAk - sumkap*trominv*wtA.t() + sumkap*trominv*ABAw + val2));
    mu = Ctk - CtBAk + sumkap*trominv*CtB*wtA.t() - sumkap*trominv*wtC.t();
  } else {
    ci = sqrt(sum(had1%x,1) + square(Ak - ABAk));
    mu = Ctk - CtBAk;
  }
  
  dsigma = hinv - sum(CtB%Ct,1);
  chi = lambda2%(dsigma + square(mu));
  
  return List::create(Named("ci") = ci,
                      Named("chi") = chi,
                      Named("dsigma") = dsigma,
                      Named("mu") = mu);
}



// [[Rcpp::export]]
List est_param3(mat xr, mat xu, vec kappa, vec m, int n, int p, vec ciold, double phi, vec chiold,
                double lambda2, vec lambdag, vec lambdagold, bool intercept, bool unpen, bool posterior, 
                bool elbo, bool start) {
  
  vec lambda2vec = lambda2*lambdag;
  int u = 0;
  int r = xr.n_cols;
  int nvars = xr.n_cols;
  mat x = xr;
  if(intercept==true && unpen==false) {
    u = 1;
    nvars += 1;
    vec evec(n,fill::ones);
    x = join_rows(evec, xr);
  } else if(unpen==true){
    u = xu.n_cols;
    nvars += xu.n_cols;
    x = join_rows(xu,xr);
  }
  
  // initialize c, chi, dsigma and mu
  double elboval;
  vec ci(n);
  vec chi(r);
  vec dsigma(nvars);
  mat sigma(nvars,nvars);
  vec mu(nvars);

  // used in all situations
  mat xrt = xr.t();
  vec h(r);
  vec om(n);
  if(start) {
    // if we are calculating starting values we need the diagonal weight 
    // matrix with wi = pi/(1-pi) and h=2*lambda
    om = ciold/(1 - ciold);
    h = 2*lambda2vec;
  } else {
    om = 0.5*m/ciold%tanh(ciold/2);
    h = lambda2vec + lambda2vec%sqrt(phi/chiold);
  }
  vec hinv = 1/h; 
  mat C = xr.each_row() % hinv.t();
  mat Ct = C.t();
  mat D = C*xrt;
  mat E(n,n), F(r,n), N(nvars,n);

  double ldetsig;
    
  if(p <= n) {
    sigma = (x.each_col() % om).t()*x;
    vec happend = h;
    happend.insert_rows(0, nvars - r, true);
    sigma.diag() += happend;
    if(elbo) {
      ldetsig = - real(log_det(sigma));
    }
    sigma = sigma.i();
    dsigma = sigma.diag();
    N = sigma*x.t();
  } else {
    if(intercept==false && unpen==false) { // checked and same as est_param and est_param2
      E = D;
      E.diag() += 1/om;
      E = E.i();
      F = Ct*E;
      N = Ct - F*D;
      if(posterior) {
        sigma = diagmat(hinv) - F*C;
        dsigma = sigma.diag();
      } else if(start==true){
        mat NW = N.each_row() % om.t();
        dsigma = sum(NW % N,1);
        N = NW*(N.t()*x.t());
      } else {
        dsigma = hinv - sum(F%Ct,1);
      }
      if(elbo) {
        ldetsig = -real(log_det(diagmat(1/om) + D)) - sum(log(om)) - sum(log(h));
      }
    } else if(intercept==true && unpen==false){ // this one is correct (compared with est_param2, est_param and R)
      double A = sum(om);
      double Ainv = 1/A;
      mat Om = diagmat(om);
      Om = Om - Ainv*om*om.t();
      E = D*Om;
      E.diag() += 1;
      E = E.i();
      F = Ct*Om*E;
      mat G = Ainv*om.t()*xr;
      mat J = C*G.t();
      mat K = G*F;
      mat M = G.each_row() % hinv.t();
      N.submat(0,0,0,n-1).fill(Ainv + conv_to<double>::from(M*G.t() - K*J));
      N.submat(0,0,0,n-1) += (G*F*D - G*Ct);
      N.submat(1,0,nvars-1,n-1).each_col() = F*J - M.t();
      N.submat(1,0,nvars-1,n-1) += Ct - F*D;
      if(posterior) {
        sigma.submat(0,0,0,0) = Ainv + conv_to<double>::from(M*G.t() - K*J);
        sigma.submat(0,1,0,nvars-1) = K*C - M;
        sigma.submat(1,0,nvars-1,0) = F*J - M.t();
        sigma.submat(1,1,nvars-1,nvars-1) = diagmat(hinv) - F*C;
        dsigma = sigma.diag();
      } else if(start==true){
        mat NW = N.each_row() % om.t();
        dsigma = sum(NW % N,1);
        N = NW*(N.t()*x.t());
      } else {
        dsigma(0) = Ainv + conv_to<double>::from(M*G.t() - K*J);
        dsigma.subvec(1,nvars-1) = hinv - sum(F % Ct,1);
      }
      if(elbo) {
        vec Dom = D*om;
        mat OmD = diagmat(1/om) + D;
        ldetsig = -real(log_det(OmD)) - sum(log(om)) - sum(log(h)) -
          log(A - conv_to<double>::from(om.t()*Dom + Dom.t()*OmD.i()*Dom));
      }
    } else { // checked (compared with calculations in R)
      mat xut = xu.t();
      mat Om = diagmat(om);
      mat xutom = (xut.each_row() % om.t());
      mat A = xutom*xu;
      mat Ainv = A.i();
      Om = Om - xutom.t()*Ainv*xutom;
      mat B = xutom*xr;
      E = D*Om;
      E.diag() += 1;
      E = E.i();
      F = Ct*Om*E;
      mat G = Ainv*B;
      mat J = G.t()*xut;
      mat K = G*F;
      mat L = K*C;
      mat M = G.each_row() % hinv.t(); // problem
      N.submat(0,0,u-1,n-1) = Ainv*xut + M*J - L*J - G*Ct + K*D;
      N.submat(u,0,nvars-1,n-1) = F*C*J - J.each_col() % hinv + Ct - F*D;
      if(posterior) {
        mat O = F*C;
        sigma.submat(0,0,u-1,u-1) = Ainv + (M - L)*G.t();
        sigma.submat(0,u,u-1,nvars-1) = -M + K*C;
        sigma.submat(u,0,nvars-1,u-1) = -M.t() + O*G.t();
        sigma.submat(u,u,nvars-1,nvars-1) = diagmat(hinv) - O;
        dsigma = sigma.diag();
      } else if(start==true){
        mat NW = N.each_row() % om.t();
        dsigma = sum(NW % N,1);
        N = NW*(N.t()*x.t());
      } else {
        dsigma.subvec(0,u-1) = Ainv.diag() + sum((M - L) % G, 1);
        dsigma.subvec(u,u+r-1) = hinv - sum(F % Ct, 1);
      }
      if(elbo) {
        mat BCt = B*Ct;
        mat OmD = diagmat(1/om) + D;
        ldetsig = -real(log_det(OmD)) - sum(log(om)) - sum(log(h)) -
          real(log_det(A - (B.each_row() % hinv.t())*B.t() + BCt*OmD.i()*BCt.t()));
      }
    }
  }

  mu = N*kappa;
  ci = sqrt(sum(x.t() % N,0).t() + square(x*mu));
  chi = lambda2vec%(dsigma.subvec(nvars-r,nvars-1) + square(mu.subvec(nvars-r,nvars-1)));

  if(elbo) {
    elboval = conv_to<double>::from(kappa.t()*x*mu + m.t()*(0.5*ci - log(exp(ci) + 1))) +
      0.5*sum(log(lambdag)) - 0.5*sqrt(phi)*sum(sqrt(chi)) -
      0.5*lambda2*sum(lambdag%(1 + sqrt(phi/chi))%(dsigma.subvec(u,u+r-1) +
      square(mu.subvec(u,u+r-1)))) + 0.5*ldetsig;
  }

  if(posterior && elbo) {
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("sigma") = sigma,
                        Named("mu") = mu,
                        Named("elbo") = elboval);
  } else if(elbo){
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("dsigma") = dsigma,
                        Named("mu") = mu,
                        Named("elbo") = elboval);
  } else if(posterior) {
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("sigma") = sigma,
                        Named("mu") = mu);
  } else {
    return List::create(Named("ci") = ci,
                        Named("chi") = chi,
                        Named("dsigma") = dsigma,
                        Named("mu") = mu);
  }

}



