#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
// MEM function with multiple arms to get draws:
// input is the following: 
// R_ctrl which is Response rate in control arm; N_ctrl which is sample size of control arm
// R_trt which is Response rate in treated arm; N_trt which is sample size of treated arm


// [[Rcpp::export]]
NumericVector cpprbinom(int n, double size, NumericVector prob) { 
  NumericVector v(n);            
  for (int i=0; i<n; i++) {v[i] = as<double>(rbinom(1, size, prob[i]));} 
  return(v); }
// [[Rcpp::export]]
int RJMCMC_c(int model_cur, NumericVector weights, 
             double prior){
  int model_prop = 1-model_cur;
  
  double num = 1.0; 
  double den = 1.0;

  if(model_prop == 0){
    num = weights[0];
    den = weights[1];
  }
  else if(model_prop == 1){
    num = weights[1];
    den = weights[0];
  }
  // Rcout << "model prop: " << model_prop << "\n";
  //Rcout << "weights 1: " << weights[0] << "\n";
  // Rcout << "weights 2: " << weights[1] << "\n";
  //Rcout << "v: " << num/den << "\n"; 
  // Rcout << "weights1/weights2: " << weights[0]/weights[2] << "\n";
  double v = num/den;
  //Rcout << "v: " << v << "\n";
  double A = std::min(1.0, v) ;
  // Rcout << "A: " << A << "\n";
  double u = runif(1)[0]; 
  // Rcout << "u: " << u << "\n";
  //bool FLAG = u <= A ;
  int ret = 0;
  if(u <= A){ret = model_prop; }
  else if(u > A) { ret = model_cur;}
  ///cout << "u <= A: " << FLAG << "\n";
  // Rcout << "ret: " << ret << "\n";
  return ret;
  
}
// [[Rcpp::export]]
NumericVector logMarg_Like_C(double Y_c, double Y_o, double N_c, double N_o, 
                   double pip, 
                   double a_c, double a_o, 
                   double b_c, double b_o){
  // Rcout << "Yc = " << Y_c << "\n";
  // Rcout << "Nc = " << N_c << "\n";
  // Rcout << "Yo = " << Y_o << "\n";
  // Rcout << "No = " << N_o << "\n";
  double a1_c = a_c+Y_c; 
  double b1_c = b_c+N_c-Y_c;
  double a1_o = a_o+Y_o;
  double b1_o = b_o+N_o-Y_o;
  double p1 = R::lbeta(a1_c, b1_c)-R::lbeta(a_c, b_c)+R::lbeta(a1_o, b1_o)-R::lbeta(a_o,b_o);
  double a2_c = a_c+Y_c+Y_o; 
  double b2_c = b_c+N_c-Y_c+N_o-Y_o;
  double p2 = R::lbeta(a2_c, b2_c)-R::lbeta(a_c, b_c);
  double norm = (1-pip)*exp(p1)+pip*exp(p2);
  double w1 = 0.0;
  double w2 = 0.0;
  // Rcout << "p1 = " << p1 << "\n";
  // Rcout << "p2 = " << p2 << "\n";
  double bic1 = -p1; 
  double bic2 = -p2;
  double min;
  if(bic1 < bic2){min = bic1;}
  else{min = bic2;}
  double p1_2 = (1-pip)*exp(-(bic1-min));
  double p2_2 = pip*exp(-(bic2-min));
  double norm_1 = p1_2+p2_2;
  double w1_2 = p1_2/norm_1;
  double w2_2= p2_2/norm_1;
  // Rcout << "norm = " << norm << "\n";
  // Rcout << "w1 (BIC)_p1 = " << w1_2 << "\n";
  // Rcout << "w2 (BIC)_p2 = " << w2_2 << "\n";
  // Rcout << "p1 (BIC)_p1 = " << p1_2 << "\n";
  // Rcout << "p2 (BIC)_p2 = " << p2_2 << "\n";
if(norm == 0){
 //model unexchangeable: the BIC value:
  //Rcout << "min = " << min << "\n";
  // Rcout << "ll2 = " << ll_2 << "\n";
  w1 = p1_2/norm_1;
  w2= p2_2/norm_1;
  // Rcout << "norm = " << norm << "\n";
  // Rcout << "w1 (BIC)_p1 = " << w1 << "\n";
  // Rcout << "w2 (BIC)_p2 = " << w2 << "\n";
  // Rcout << "p1 (BIC)_p1 = " << p1_2 << "\n";
  // Rcout << "p2 (BIC)_p2 = " << p2_2 << "\n";
  
  // double p1_3 = (1-pip)*exp(p1+100);
  // double p2_3 = pip*exp(p2+100);
  // double norm_2 = p1_3+p2_3;
  // w1 = p1_3/norm_2;
  // w2= p2_3/norm_2;
  // // Rcout << "norm = " << norm << "\n";
  // Rcout << "w1 (/10) = " << w1 << "\n";
  // Rcout << "w2 (/10) = " << w2 << "\n";
  // Rcout << "p1 (/10) = " << p1_2 << "\n";
  // Rcout << "p2 (/10) = " << p2_2 << "\n";
}
 else{
    w1 = ((1-pip)*exp(p1))/norm; 
    w2 = ((pip)*exp(p2) )/norm; 
    // Rcout << "w1 (RJMCMC) = " << w1 << "\n";
    // Rcout << "w2 (RJMCMC) = " << w2 << "\n";
    // Rcout << "p1 (RJMCMC) = " << p1 << "\n";
    // Rcout << "p2 (RJMCMC) = " << p2 << "\n";
  }
 
 // Rcout << "w1 = " << w1 << "\n";
 // Rcout << "w2 = " << w2 << "\n";
  
  NumericVector results = {w1,w2};

  return results;
}
// [[Rcpp::export]]
double log_odds_ratio(double pt, double pc)
                      {
  double results = log(pt)-log(1-pt)+log(1-pc)-log(pc);
  
  return results;
}

// [[Rcpp::export]]
List MEM_binary(double r_ctrl, double n_ctrl, 
                  double r_trt, double n_trt, 
                  double ext_r_ctrl, double ext_N_ctrl, 
                  double ext_r_trt, double ext_N_trt, 
                  int nburn, int M, 
                  double pip_ctrl, double pip_trt,
                  double a_c, double b_c, 
                  double a_t, double b_t, 
                  double a_oc, double b_oc,
                  double a_ot, double b_ot) {
  //set up all draws: 
  int Mall = M+nburn;
  // Initialize storage for draws:
  NumericVector draws_pc(Mall);
  NumericVector draws_pt(Mall);
  NumericVector draws_logOR(Mall);
  NumericVector draws_mem_IND_ctrl(Mall); 
  NumericVector draws_mem_IND_trt(Mall); 

  // Starting values for the draws/components:
  draws_pc[0] = 0.5;
  draws_pt[0] = 0.5;
  draws_logOR[0] = 1; 
  draws_mem_IND_ctrl[0] = 0; 
  draws_mem_IND_trt[0] = 0;
  
  //compute the weights from the log marginal likelihood for the control external data:
  NumericVector log_marg_ctrl_cur = logMarg_Like_C(r_ctrl, ext_r_ctrl, n_ctrl, ext_N_ctrl, 
                                                   pip_ctrl, 
                                                   a_c, a_oc, 
                                                   b_c, b_oc); 
  //compute the weights from the log marginal likelihood for the treated external data:
  NumericVector log_marg_trt_cur = logMarg_Like_C(r_trt, ext_r_trt, n_trt, ext_N_trt, 
                                                  pip_trt, 
                                                  a_t, a_ot, 
                                                  b_t, b_ot); 
  // Gibbs sampler -- big loop
  for (int s = 1; s < Mall; s++) {
    //current exchangeability indicator for control:
    double mem_ctrl = draws_mem_IND_ctrl[s-1];
    
    //current exchangeability indicator for treated:
    double mem_trt = draws_mem_IND_trt[s-1];
    
    // sampling current values control arm: 
    double a_ctrl = a_c+r_ctrl+ext_r_ctrl*mem_ctrl; 
    double b_ctrl = b_c+(n_ctrl - r_ctrl)+(ext_N_ctrl - ext_r_ctrl)*mem_ctrl; 
    //Rcout << "a_ctrl: " << a_ctrl << "\n";
    //Rcout << "b_ctrl: " << b_ctrl << "\n";
    draws_pc[s] = rbeta(1, a_ctrl, b_ctrl)[0];
    
    // sampling current values treated arm: 
    double a_trt = a_t+r_trt+ext_r_trt*mem_trt; 
    double b_trt = a_t+(n_trt - r_trt)+(ext_N_trt - ext_r_trt)*mem_trt; 
    draws_pt[s] = rbeta(1, a_trt, b_trt)[0];
    
    //sampling exchangeability indicator for control: 
  
    //draws_weights_ctrl[s, _] = log_marg_ctrl_cur; 
    
    //sample from the exchangeability indicator using RJMCMC: 
    int mem_ctrl_update = RJMCMC_c(mem_ctrl, log_marg_ctrl_cur, pip_ctrl); 
    draws_mem_IND_ctrl[s] = mem_ctrl_update; 
    
    //sampling exchangeability indicator for treated: 
    //draws_weights_trt[s, _] = log_marg_trt_cur; 
    
    //sample from the exchangeability indicator using RJMCMC: 
    int mem_trt_update = RJMCMC_c(mem_trt, log_marg_trt_cur, pip_trt); 
    draws_mem_IND_trt[s] = mem_trt_update; 
    
    draws_logOR[s] = log_odds_ratio(draws_pt[s], draws_pc[s]); 

    
  }

  
  
  List result; 
  result["pc"] = draws_pc ; 
  result["pt"] = draws_pt;
  result["logOR"] = draws_logOR;
  result["mem_ctrl"] = draws_mem_IND_ctrl;
  result["mem_trt"] = draws_mem_IND_trt;
  result["weights_ctrl"] = log_marg_ctrl_cur;
  result["weights_trt"] = log_marg_trt_cur;
  
  //Rcout << "pc value" << pc << "/n";
  return result;
}
