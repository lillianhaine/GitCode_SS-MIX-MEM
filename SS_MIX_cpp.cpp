#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector cpprbinom(int n, double size, NumericVector prob) { 
  NumericVector v(n);            
  for (int i=0; i<n; i++) {v[i] = as<double>(rbinom(1, size, prob[i]));} 
  return(v); }


// [[Rcpp::export]]
List SS_MIX(NumericVector ps_trial,  
                    NumericVector ps_ext, NumericVector y_ext,
                    int Mall,
              double sd_muct,
              double a_tauct,
              double b_tauct,
              double sd_muu,
              double a_tauu,
              double b_tauu,
              double a_w, 
              double b_w) {
  //inputs are: the x vectors (data), 
  //the number of draws (Mall),
  //AND the hyperparameters: 
  //the mean/sd for muct or component1
  //the mean/sd for muu or component2
  //the shape parameters for the beta w
  //the rate/scale parameters for tauct
  //the rate/scale parameters for tauu
  
  // values that will help for future comps:
  double nct = ps_trial.size();
  double n_ext = y_ext.size();
  double mean_ct = mean(ps_trial); 
  //double var_ct = var(x_trial);
  //double var_2 = var(x_obs);
  
  // Initialize storage for draws:
  NumericVector draws_muct(Mall);
  NumericVector draws_muu(Mall);
  NumericVector draws_tauct(Mall);
  NumericVector draws_tauu(Mall);
  NumericVector draws_w(Mall);
  NumericVector draws_nobsIN(Mall); 
  NumericVector draws_yobsIN(Mall);
  NumericMatrix draws_Zs(Mall, n_ext);
  
  // Starting values for the draws/components:
  draws_muct[0] = 0.0;
  draws_muu[0] = -3.0;
  draws_tauct[0] = 1.0;
  draws_tauu[0] = 1.0;
  draws_w[0] = 0.5;
  draws_nobsIN[0] = n_ext/2.0;
  draws_yobsIN[0] = sum(y_ext)/2.0;
  draws_Zs(0, _) = rep(0, n_ext); 
  //Rcout << "rows of Z: " << draws_Zs.nrow() << "\n";
  //Rcout << "cols of Z: " << draws_Zs.ncol() << "\n";
  // Gibbs sampler -- big loop
  // ATTENTION: We will still have Mall iterations, with s going from 0 to S - 1
  for (int s = 1; s < Mall; s++) {
    // previous iterations values (gibbs so used in future sampling):
    double muct_cur = draws_muct[s - 1];
    double tauct_cur = draws_tauct[s - 1];
    double muu_cur = draws_muu[s - 1];
    double tauu_cur = draws_tauu[s - 1];
    double w_cur = draws_w[s - 1];
    double sdct_cur = sqrt(1.0/tauct_cur);
    double sdu_cur = sqrt(1.0/tauu_cur);
    //Rcout << "The value of muu: " << muu_cur << "\n";
    //Rcout << "The value of mct: " << muct_cur << "\n";
    
    // w vector and new w value for current draw, then Z_vec made:
    NumericVector wi1 = w_cur*dnorm(ps_ext, muct_cur, sdct_cur);
    NumericVector wi2 = (1-w_cur)*dnorm(ps_ext, muu_cur, sdu_cur);
    NumericVector norm = wi1+wi2;
    NumericVector w_vec = wi1/norm; 
    NumericVector Z_vec = cpprbinom(n_ext, 1, w_vec);
    draws_Zs(s, _) = Z_vec;
    //Rcout << "s: " << s <<  "\n";
    // Rcout << "current w: " << w_cur <<  "\n";
    // Rcout << "current muct: " << muct_cur <<  "\n";
    // Rcout << "current sdct: " << sdct_cur <<  "\n";
    // Rcout << "current muu: " << muu_cur <<  "\n";
    // Rcout << "current sdu: " << sdu_cur <<  "\n";
    // Rcout << "wi1: " << wi1 <<  "\n";
    //int NA_vals = sum(is.na(Z_vec));
    //get nobsIN and nobsOUT from Z_vec
    double nobsIN = sum(Z_vec);
    double nobsOUT = Z_vec.size()-nobsIN;
    double yobsIN = sum(Z_vec*y_ext);
    draws_nobsIN[s] = nobsIN;
    draws_yobsIN[s] = yobsIN;
    // Rcout << "nobsIN: " << nobsIN << "\n";
    // Rcout << "YobsIN: " << yobsIN << "\n";
    
    //Rcout << "NA in Z vec: " << NA_vals << "/n";
    NumericVector xobsIN = ps_ext[Z_vec == 1];
    NumericVector xobsOUT = ps_ext[Z_vec == 0];
    //sampling and saving for for w:
    double a_w_new = a_w+nobsIN;
    double b_w_new = b_w+nobsOUT;
    double w_new = rbeta(1, a_w_new, b_w_new)[0];
    draws_w[s] = w_new;

    //Rcout << "a_w: " << a_w_new << "\n";
    //Rcout << "b_w: " << b_w_new << "\n";
    //if( R_IsNA(nobsIN) ){ nobsIN == 0}
    
    // sampling from the mean and the precision for the trial (note that the rgamma in cpp takes the scale not the rate)
    double shape_ct = ((nobsIN+nct)/2.0) + a_tauct;
    double rate_ct =  b_tauct + (0.5 * ( sum(pow( (xobsIN-muct_cur), 2) )+sum(pow( (ps_trial-muct_cur), 2)) ) );
    //double val_ct_check = sum(pow( (xobsIN-muct_cur), 2) )+sum(pow( (x_trial-muct_cur), 2));  
    double tauct_new = rgamma(1, shape_ct, 1.0 / rate_ct)[0];
    draws_tauct[s] = tauct_new;

    // Rcout << "nct: " << nct << "\n";
    // Rcout << "a_tauct: " << a_tauct << "\n";
    // Rcout << "val check: " << val_ct_check<< "\n";
    // Rcout << "shape_ct: " << shape_ct << "\n";
    // Rcout << "scale_ct: " << 1.0/rate_ct << "\n";
    // Rcout << "rate_ct: " << rate_ct << "\n";
    // Rcout << "tau_ct: " << tauct_new << "\n";
    // Generate and save muct. Note that rnorm takes SD, not variance
    double mean_xobsIN = mean(xobsIN);
    double mean_xobsOUT = mean(xobsOUT);
    bool NA_val1 = NumericVector::is_na(mean_xobsIN);
    if( NA_val1){mean_xobsIN = 0; }
    bool NA_val2 = NumericVector::is_na(mean_xobsOUT);
    if( NA_val2){mean_xobsOUT = 0; }
    // Rcout << "mean_XobsIN: " << mean_xobsIN << "\n";
    // Rcout << "size_XobsIN: " << xobsIN.size() << "\n";
    // Rcout << "mean_XobsOUT: " << mean(xobsOUT) << "\n";
    // Rcout << "size_XobsOUT: " << xobsOUT.size() << "\n";
    double mean_muct = tauct_new*(nct*mean_ct+nobsIN*mean_xobsIN)/(sd_muct+tauct_new*(nct+nobsIN));
    double prec_muct = sd_muct+ ( tauct_new*(nct+nobsIN) );
    double var_muct = 1.0/prec_muct;
    NumericVector muct_new = rnorm(1, mean_muct, sqrt( var_muct) );
    draws_muct[s] = muct_new[0];
    
    //Rcout << "mean_ct: " << mean_ct << "\n";
    // Rcout << "mean_XobsIN: " << mean(xobsIN) << "\n";
    //Rcout << "sd_muct: " << sd_muct << "\n";
    //Rcout << "mean_muct: " << mean_muct << "\n";
    //Rcout << "prec_muct: " << prec_muct<< "\n";

    
    // sampling from the mean and the precision for component2 (muu) (note that the rgamma in cpp takes the scale not the rate)
    double shape_u = ( ((nobsOUT)/2.0) + a_tauu );
    double rate_u =  b_tauu + 0.5 * (sum(pow( ( xobsOUT-muu_cur ), 2) ) );
    double tauu_new = rgamma(1, shape_u, 1.0 / rate_u)[0];
    draws_tauu[s] = tauu_new;
      
    // Generate and save muct. Note that rnorm takes SD, not variance
    double mean_muu = tauu_new*(nobsOUT*mean_xobsOUT)/(sd_muu+tauu_new*(nobsOUT));
    double prec_muu = sd_muu+( tauu_new*(nobsOUT) );
    double var_muu = 1.0/prec_muu;
    NumericVector muu_new = rnorm(1, mean_muu, sqrt(var_muu)  );
    draws_muu[s] = muu_new[0];
    // 

    // Rcout << "w_new: " << w_new << ", w_old: " <<w_cur << "\n";
    // Rcout << "shape_tauct: " << shape_ct << ", rate_ct: " << rate_ct << "\n";
    // Rcout << "scale_tauct: " << 1.0/rate_ct << "\n";
    // Rcout << "tauct_new: " << tauct_new << ", tauct_old: " << tauct_cur << "\n";
    // Rcout << "muu_new: " << muu_new << ", muu_old: " <<muu_cur << "\n";
    // Rcout << "tauu_new: " << tauu_new << ", tauu_old: " <<tauu_cur << "\n";
    //Rcout << "The value of mct: " << muct_cur << "\n";

  }
  
  List result; 
  result["muct"] = draws_muct ; 
  result["muu"] = draws_muu ;
  result["tauct"] = draws_tauct;
  result["tau2"] = draws_tauu; 
  result["w"] = draws_w; 
  result["yobsIN"] = draws_yobsIN;
  result["nobsIN"] = draws_nobsIN;
  result["Zvec"] = draws_Zs; 
  return result;
}


