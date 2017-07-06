
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_SPECIAL_FUNCTION_HH
#define SHTA_SPECIAL_FUNCTION_HH 1

namespace shta {
  
  double Skellam(double x, double p, double n){
    /*
      ridge of Bassel, strainforward sum of poisson, x = x_p - x_n 
      in area of x_n : [p - 2*\qsrt{p}, p + 2*\qsrt{p}],  [n - 2*\qsrt{n}, n + 2*\qsrt{n}] only
    */
    double result = 0.;
    double x_n = 0;
    double dx_n = 0.05;
    double sigma_p = 2*sqrt(p);
    double sigma_n = 2*sqrt(n);
    double x_n_min = std::min(p-sqrt(p), n-sqrt(n));
    double x_n_max = std::max(p+sqrt(p), n+sqrt(n));
    if(x_n_min < 0) x_n_min = 0;
    for(double x_n = x_n_min; x_n < x_n_max; x_n += dx_n){
      result += TMath::Poisson(x + x_n, p) * TMath::Poisson(x_n, n);
    }

    return result;
  }

  double Heaviside(double x){
    if(x >= 0.) return 1.0;
    return 0.;
  }
  
  double Heaviside(double x, double y){
    if(x >= y) return 1.0;
    return 0.;
  }
}

#endif



