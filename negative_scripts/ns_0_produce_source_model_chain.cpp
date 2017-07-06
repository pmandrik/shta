
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#include "shta_core.hh"

using namespace shta;
using namespace std;

std::vector<shta::MarkovChain *> *  GenerateMC(Function * pdf, int n_chains, int lenght){
  shta::ProposalGauss * prop = new shta::ProposalGauss();
  prop->step_factor = 0.08;
  shta::MetropolisHastings mh = shta::MetropolisHastings(pdf, dynamic_cast<shta::ProposalFunction *>(prop) );
  std::vector<shta::MarkovChain *> * chains_from_true = new std::vector<shta::MarkovChain *>;

  for(int i = 0; i < n_chains; i++){
    shta::MarkovChain * chain = mh.ConstructChain(lenght);
    chains_from_true->push_back( chain );
  }
  return chains_from_true;
}

void WriteMC(string file_name, vector<shta::MarkovChain *> * chains, Function * pdf, bool force_rewrite){
  TFile * out_file = new TFile(file_name.c_str(), "READ");
  if(out_file->IsOpen() and not force_rewrite){
    msg("file already exist",  file_name,"; return");
    return;
  }
  delete out_file;
  out_file = new TFile(file_name.c_str(), "RECREATE");

  auto vars = pdf->used_variables;
  for(int i = 0; i < chains->size(); i++){
    chains->at(i)->Write("base_chain_" + std::to_string(i));

    for(auto var : vars){
      TH1D * hist = chains->at(i)->DrawVariable(var, 500);
      hist->Write();
    }

    out_file->cd();
  }

  auto pars = pdf->used_parameters;
  TTree * generation_tree = new TTree("gen_description", "gen_description");
  vector<string> names; 
  vector<double> values;
  for(auto var : pars){
    names.push_back(var->name);
    values.push_back(var->value);
  }

  generation_tree->Branch("gen_par_names", &names);
  generation_tree->Branch("gen_par_values", &values);

  int size = chains->size();
  chains->at(0)->FindChainWeightedLenght();
  int lenght = chains->at(0)->GetChainWeightedLenght();
  msg("n_chains, lenght = ", size, lenght);
  generation_tree->Branch("n_chains", &size);
  generation_tree->Branch("lenght", &lenght);

  generation_tree->Fill();
  generation_tree->Write();

  out_file->Close();
}

int main(int argc, char** argv){
  /*
    This script will generate MCMC data from simple single dimension model to study negative weight treatment.
    We will suppouse that T_m is a effective number of weighted events with negative weight factors, and T_p with positive, so
    raw prediction is Poisson(x, T_p-T_m). Than we introduce into model a parameter P and barlow-beeston priors Poisson(t_p,T_p), Poisson(t_m,T_m).
    Rather than \mu = P*(T_p-T_m) we will use \mu = P*(1 + T_p-T_m) to avoid long posterior probability tail for P araised from T_p-T_m ~ 0.
  */

  shta::FunctionBuilder fbuilder;
  auto x   = fbuilder.DefineVariable( "x"  ,  5, 0., 100,  shta::vartype::VARIABLE );
  auto t_p = fbuilder.DefineVariable( "t_p",  5, 0., 100,  shta::vartype::VARIABLE );
  auto t_m = fbuilder.DefineVariable( "t_m",  5, 0., 100,  shta::vartype::VARIABLE );
  auto P   = fbuilder.DefineVariable( "P",    3, 0., 100,  shta::vartype::PARAMETER );
  auto T_p = fbuilder.DefineVariable( "T_p",  9, 0., 100,  shta::vartype::PARAMETER );
  auto T_m = fbuilder.DefineVariable( "T_m",  6, 0., 100,  shta::vartype::PARAMETER );
  auto PDF_BASE = fbuilder.DefineFunction("pdf_true",  "TMath::Poisson(x, P*(4 + T_p-T_m)) * TMath::Poisson(t_p,T_p) * TMath::Poisson(t_m,T_m)");

  msg("build model ... ");
  fbuilder.BuildFunctions("negative_scripts");

  msg("generate and write mc data ... ");
  int number_of_chains = 100;
  int chains_lenght    = 100000;

  // region A : T_p ~ 3 * T_m
  T_p->value = 12;
  T_m->value = 4;
  auto chains_A = GenerateMC(PDF_BASE, number_of_chains, chains_lenght);
  for(auto chain : *chains_A) chain->Print();
  WriteMC("ns_0_PDF_BASE_a.root", chains_A, PDF_BASE, true);

  // region B : T_p ~ T_m
  T_p->value = 9;
  T_m->value = 7;
  auto chains_B = GenerateMC(PDF_BASE, number_of_chains, chains_lenght);
  for(auto chain : *chains_B) chain->Print();
  WriteMC("ns_0_PDF_BASE_b.root", chains_B, PDF_BASE, true);

  return 0;
}









