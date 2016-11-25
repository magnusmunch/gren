/******************************* Preamble **************************
* code belonging to ENVB2.pdf					   *						
* version: 01							   *
* author: Magnus Münch						   *
* created: 25-11-2016						   *
* last edited: 25-11-2016					   *
*******************************************************************/

/******************************* Notes *****************************
*                                                                  *
*******************************************************************/



#include <math.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List est_param(mat x, vec kappa, vec m, int n, int p, vec ciold, vec phi, vec chiold, vec lambda2, bool intercept) {
    vec evec(n,fill::ones);
    mat xaug = join_rows(evec,x);
    mat trx = x.t();
    vec h = (1 + phi/chiold)%lambda2;  
    vec om = 0.5*(m/ciold) % tanh(ciold/2);
    vec omsq = sqrt(om);
    vec hinv = 1/h;

    double invtrom;
    mat Omadj(n,n);
    mat xhinv(n,p);
    mat xtrxom(n,n);
    mat Ainv(p,p);
    vec ainvxom(p);
    mat sigma(p + 1,p + 1);
    vec sigmapart(p);
    vec mu(p + 1);
    vec varmean;
    
    if (intercept) {
      invtrom = 1/sum(om);
      Omadj = diagmat(om) - invtrom*(om*om.t());
      xhinv = x.each_row() % hinv.t();
      xtrxom = x*xhinv.t()*Omadj;
      xtrxom.diag() += 1;
      Ainv= diagmat(hinv) - xhinv.t()*Omadj*xtrxom.i()*xhinv;
      ainvxom = Ainv*x.t()*om;
      sigmapart = -invtrom*ainvxom;
      sigma(0,0) = as_scalar(invtrom + pow(invtrom,2.0)*(om.t()*x*ainvxom));
      sigma.submat(1,0,p,0) = sigmapart;
      sigma.submat(0,1,0,p) = sigmapart.t();
      sigma.submat(1,1,p,p) = Ainv;
      mu = sigma*xaug.t()*kappa;
      ciold = pow(diagvec(xaug*sigma*xaug.t()) + pow(xaug*mu,2.0),0.5);
      varmean = diagvec(sigma) + pow(mu,2.0);
      chiold = lambda2%varmean.rows(1,p);
    } else {
      xaug = x.each_col() % omsq;
      xhinv = x.each_row() % hinv.t();
      xtrxom = xhinv * xaug.t();
      xtrxom.diag() += 1;
      sigma = diagmat(hinv) - xhinv.t()*xtrxom.i()*xhinv;
      mu = sigma*x.t()*kappa;
      ciold = pow(diagvec(x*sigma*x.t()) + pow(x*mu,2.0),0.5);
      varmean = diagvec(sigma) + pow(mu,2.0);
      chiold = lambda2%varmean;
    }

    return List::create(Named("sigma") = sigma,
                        Named("mu") = mu,
                        Named("ci") = ciold,
                        Named("chi") = chiold);
}

