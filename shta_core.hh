
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 


#ifndef SHTA_CORE_HH
#define SHTA_CORE_HH 1

// std
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>
#include <map>
#include <forward_list>
#include <dlfcn.h>

// ROOT
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"
#include "TColor.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TArrow.h"
#include "TTree.h"
#include "TFile.h"
#include "THStack.h"

#include "TApplication.h"

// shta random

TRandom3 _random_engine;

// is_eq_one
template<typename T, typename B> bool is_eq_one(T t, B f){
  return t == f;
}

template<typename T, typename B, typename... Args> bool is_eq_one(T t, B f, Args... args){
  if(t == f) return true;
  return is_eq_one(t, args...);
}

// is_eq_all

std::string get_shta_src_path(){
  return std::string( SHTA_SOURCE_DIR );
}

std::string get_shta_exe_path(){
  return std::string( SHTA_BINARY_DIR );
}

std::string get_shta_userfunctions_path(){
  return get_shta_src_path() + "/build-in-functions";
}

// TODO remove from core -------------------- root

namespace root {
  static int name_id_counter = 0;

  std::string get_date(){
    std::time_t result = std::time(nullptr);
    std::string name = std::string(std::asctime(std::localtime(&result)));
    name.pop_back();
    return name;
  }

  std::string get_name(){
    return "#" + std::to_string(name_id_counter++) + " " + get_date();
  }

  TCanvas * get_canvas (int width=640, int height=480){
    std::string name = get_name();
    return new TCanvas(name.c_str(), name.c_str(), width, height);
  }

  double get_th1d_limit_upper(TH1D * hist, double frac){
    int bin, nbins = hist->GetNbinsX();
    double subsum = 0;
    double sum = hist->Integral();
    for(bin = nbins; bin > 0; bin--){
      subsum += hist->GetBinContent( bin );
      if(subsum >= sum*frac) break;
    }
    return hist->GetBinLowEdge(bin);
  }

  void get_th2d_neyman_pallete(){
    const Int_t Number = 3;
    Double_t Red[Number]    = { 1.0, 0.05, 0.30};
    Double_t Green[Number]  = { 1.0, 0.10, 0.60};
    Double_t Blue[Number]   = { 1.0, 0.60, 1.00};
    Double_t Length[Number] = { 0.00, 0.50, 1.00};
    Int_t nb=50;
    TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
  }

  class TGraphShade{
    public:
      TGraphShade (unsigned long long int npts) : npoints(npts) {
        l = new TGraph(npoints);
        r = new TGraph(npoints);
        sh = new TGraph(2 * npoints);
      }

      ~TGraphShade(){
        delete l;
        delete r;
        delete sh;
      }

      void SetPoint(int indx, double x1, double y1, double x2, double y2){
        l->SetPoint(indx, x1, y1);
        r->SetPoint(indx, x2, y2);
        sh->SetPoint(indx, x1, y1);
        sh->SetPoint(2*npoints - indx, x2, y2);
      }

      unsigned long long int npoints;
      TGraph * l, * r, * sh;
  };
};

#include "shta_msg.hh"
namespace global_option{
  int def_verbose_lvl = shta::verbose::INFO;
}

class BaseClass{
  public:
    BaseClass(){
      verbose_lvl = global_option::def_verbose_lvl;
    }
    int verbose_lvl;
};

#include "shta_special_functions.hh"
#include "shta_function.hh"
#include "shta_mcmc.hh"
#include "shta_neyman.hh"
#include "shta_CLs.hh"
#include "shta_hep_model.hh"

#endif


