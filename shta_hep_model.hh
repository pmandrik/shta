
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_HEP_MODEL_HH
#define SHTA_HEP_MODEL_HH 1

namespace shta {

  struct HepTemplate{
    HepTemplate(TH1D * h) : hist(h) {}

    void AddMultFactor(Variable * var){ simple_factors.push_back( var ); }
    void AddMultFactor(Function * var){ complex_factors.push_back( var );}

    TH1D * hist;
    std::vector<Variable*> simple_factors;
    std::vector<Function*> complex_factors;
  };

  class ModelBuilder : public BaseClass {
  };

  class HepModelPoisson : public ModelBuilder {
    public:
      HepModelPoisson(std::string model_name){
        name = model_name;
        hist_data = nullptr;
      }

      void SetData(TH1D* hist){
        nbins = hist->GetXaxis()->GetNbins();
        hist_data = hist;
      }

      void AddTemplate(HepTemplate * temp){
        templates.push_back(temp);
      }

      void AddConstrainTerm(Function * f){
        constrains.push_back( f );
      }

      std::string GetTemplateBinContribution(HepTemplate * temp, int bin){
        std::string vars_mult = std::to_string( temp->hist->GetBinContent(bin+1) ); // FIXME

        for(auto var : temp->simple_factors)
          vars_mult += " * " + var->GetName();

        for(auto func : temp->complex_factors)
          vars_mult += " * " + func->GetDefinition();

        return vars_mult;
      }

      Function * GetLikelihood(){
        if(!hist_data){
          return nullptr;
        }

        std::string sb_string_mu = "";
        for(int b = 0; b < nbins; b++){
          double bin_i_data = hist_data->GetBinContent(b+1);
          std::string st_string_mu = "";
          for(int t = 0; t < templates.size(); t++){
            HepTemplate * temp = templates.at(t);
            std::string mu_name = "_mu_"+std::to_string(t)+"_"+std::to_string(b);
            std::string part = GetTemplateBinContribution(temp, b);
            fbuilder.DefineFunction(mu_name, part);
            if(t) st_string_mu += "+";
            st_string_mu += mu_name;
          }

          std::string mu_name = "_mu_"+std::to_string(b);
          fbuilder.DefineFunction(mu_name, st_string_mu);

          if(b) sb_string_mu +="*";
          sb_string_mu += "TMath::Poisson(" + std::to_string(bin_i_data) + ", " + mu_name + ")";
        }

        for(auto constrain : constrains){
          sb_string_mu+= " * " + constrain->GetName();
        }
        //std::string sb_string_tot = "TMath::Log( " + sb_string_mu + " )";
        std::string sb_string_tot = sb_string_mu;
        auto lhood = fbuilder.DefineFunction(name + "_nll", sb_string_tot);

        fbuilder.BuildFunctions(name + "_nll_def");
        return lhood;
      }

      Function * GetNLLikelihood(bool use_data=true){
        if(!hist_data and use_data){
          MSG_WARNING("shta::HepModelPoisson : GetLikelihood() : data histogram is not defined, return");
          return nullptr;
        }

        std::string sb_string_mu = "";
        for(int b = 0; b < nbins; b++){
          std::string st_string_mu = "";
          for(int t = 0; t < templates.size(); t++){
            HepTemplate * temp = templates.at(t);
            std::string mu_name = "_mu_"+std::to_string(t)+"_"+std::to_string(b);
            std::string part = GetTemplateBinContribution(temp, b);
            fbuilder.DefineFunction(mu_name, part);
            if(t) st_string_mu += "+";
            st_string_mu += mu_name;
          }
          std::string mu_name = "_mu_"+std::to_string(b);
          fbuilder.DefineFunction(mu_name, st_string_mu);

          if(use_data) sb_string_mu += " - (" + mu_name + ") + " + std::to_string( hist_data->GetBinContent(b+1) ) + " * TMath::Log(" + mu_name + ")";
          else         sb_string_mu += " - (" + mu_name + ") + " + "__data_bin_" + std::to_string(b)               + " * TMath::Log(" + mu_name + ")";
        }

        if(not use_data)
          for(int b = 0; b < nbins; b++)
            fbuilder.DefineVariable( "__data_bin_" + std::to_string(b), 1, 0., 100., shta::vartype::VARIABLE); // FIXME plz ...

        std::string sb_string_tot = sb_string_mu;
        for(auto constrain : constrains){
          sb_string_tot += " + TMath::Log(" + constrain->GetName() + ")";
        }

        auto lhood = fbuilder.DefineFunction(name + "_nll", sb_string_tot);
        fbuilder.BuildFunctions(name + "_nll_def");
        return lhood;
      }

      shta::FunctionBuilder fbuilder;

      TH1D* hist_data;
      std::vector<HepTemplate*> templates;
      std::vector<Function*> constrains;

      int nbins;
      std::string name;
  };
};

#endif








