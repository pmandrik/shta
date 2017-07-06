
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#include "shta_core.hh"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

using namespace shta;
using namespace std;

struct Limit{
  public:
  Limit(std::string n, TTree * t, bool add_or_read = true){
    name = n;
    tree = t;
    if(add_or_read) AddToTree();
    else            ReadFromTree();
  }

  void GetLimits(TH1D * hist){
    double CL_2sigma = 1. - 0.0455;
    l_2s = root::get_th1d_limit_upper(hist, CL_2sigma + (1.-CL_2sigma)*0.5);
    c_2s = root::get_th1d_limit_upper(hist, 0.5);
    u_2s = root::get_th1d_limit_upper(hist, (1.-CL_2sigma)*0.5);
    double CL_1sigma = 1. - 0.3173;
    l_1s = root::get_th1d_limit_upper(hist, CL_1sigma + (1.-CL_1sigma)*0.5);
    c_1s = root::get_th1d_limit_upper(hist, 0.5);
    u_1s = root::get_th1d_limit_upper(hist, (1.-CL_1sigma)*0.5);
  }

  void AddToTree(){
    tree->Branch( (name + "_l_1s").c_str(),   &l_1s  );
    tree->Branch( (name + "_c_1s").c_str(),   &c_1s  );
    tree->Branch( (name + "_u_1s").c_str(),   &u_1s  );
    tree->Branch( (name + "_l_2s").c_str(),   &l_2s  );
    tree->Branch( (name + "_c_2s").c_str(),   &c_2s  );
    tree->Branch( (name + "_u_2s").c_str(),   &u_2s  );
  }

  void ReadFromTree(){
    tree->SetBranchAddress( (name + "_l_1s").c_str(),   &l_1s  );
    tree->SetBranchAddress( (name + "_c_1s").c_str(),   &c_1s  );
    tree->SetBranchAddress( (name + "_u_1s").c_str(),   &u_1s  );
    tree->SetBranchAddress( (name + "_l_2s").c_str(),   &l_2s  );
    tree->SetBranchAddress( (name + "_c_2s").c_str(),   &c_2s  );
    tree->SetBranchAddress( (name + "_u_2s").c_str(),   &u_2s  );
  }

  double GetValue(string key){
    if( key == name + "_l_1s" ) return l_1s;
    if( key == name + "_c_1s" ) return c_1s;
    if( key == name + "_u_1s" ) return u_1s;
    if( key == name + "_l_2s" ) return l_2s;
    if( key == name + "_c_2s" ) return c_2s;
    if( key == name + "_u_2s" ) return u_2s;
  }

  TTree * tree;
  std::string name;
  double l_1s, c_1s, u_1s, l_2s, c_2s, u_2s, def_val;
};

std::vector<shta::MarkovChain *> * GetBaseMC(string file_name){
  TFile * in_file = new TFile(file_name.c_str(), "READ");
  if(not in_file->IsOpen()) {
    msg("cant open file", file_name);
    return nullptr;
  }

  TTree * generation_tree = (TTree*)in_file->Get("gen_description");
  int size, lenght;
  generation_tree->SetBranchAddress("n_chains", &size);
  generation_tree->SetBranchAddress("lenght", &lenght);
  generation_tree->GetEntry(0);
  msg("size, lenght = ", size, lenght);

  std::vector<shta::MarkovChain *> * chains_from_true = new std::vector<shta::MarkovChain *>;
  for(int i = 0; i < size; i++){
    shta::MarkovChain * chain = new shta::MarkovChain();
    chain->Read("base_chain_" + std::to_string(i));
    chains_from_true->push_back( chain );
    in_file->cd();
  }

  in_file->Close();
  delete in_file;
  return chains_from_true;
}

void FitAndWrite(std::vector<shta::MarkovChain *> * chains_from_true, shta::Function * model, string fname, int n_to_burn, int save_every_divisor, int model_chain_lenght, int max_number_of_points_to_read, int max_number_of_chains_to_read, shta::ProposalFunction * prop, int limit_hist_bins){
  shta::MetropolisHastings mh_model = shta::MetropolisHastings(model, prop);
  mh_model.verbose_lvl = shta::verbose::ERROR;

  TFile * out_file = new TFile(fname.c_str(), "RECREATE");
  TTree * output_data = new TTree("ttree","ttree");

  auto fchain = chains_from_true->at(0);
  auto parameters = &(model->used_parameters);
  auto variables  = &(model->used_variables);
  int n_pars = parameters->size();
  int n_vars = variables->size();

  int weight_index = chains_from_true->at(0)->GetLabelIndex("MCMC_weight");
  std::vector<int> indexes;
  for(auto par : *parameters)
    indexes.push_back( fchain->GetLabelIndex( par->name ) );

  std::vector<Limit*> limits;
  for(auto var : *variables)
    limits.push_back( new Limit(var->name, output_data) );

  double weight;
  double chain_id;
  output_data->Branch("chain_id",     &chain_id    );
  output_data->Branch("weight",       &weight      );

  shta::msg("main loop ... ");
  shta::MarkovChain * chain_from_model;
  int start_number = 25;
  for(chain_id = start_number; chain_id < chains_from_true->size(); chain_id++){
    if(chain_id >= start_number+max_number_of_chains_to_read) break;
    shta::msg("begin base chain", chain_id, "... ");
    std::vector <std::vector<double> * > * points = & (chains_from_true->at(chain_id)->chain);
    int chain_size = points->size();

    // 1. burn-in
    chains_from_true->at(chain_id)->Print(10);
    int index = chains_from_true->at(chain_id)->GetIndexAfterBurn(n_to_burn, weight);

    // 2. produce chain for every model
    while(index < chain_size and index < max_number_of_points_to_read){
      bool save_hist = false;
      if( (1+index) % save_every_divisor == 0) save_hist = true;
      std::vector<double> * point = points->at(index);
      weight = point->at(weight_index);

      for(int i = n_pars-1; i >= 0 ; i--)
        parameters->at(i)->value = point->at( indexes[i] );

      chain_from_model = mh_model.ConstructChain(model_chain_lenght);
      for(int i = n_vars-1; i >= 0 ; i--){
        Limit * limit = limits.at(i);
        Variable * var = variables->at(i);
        TH1D * hist = var->GetTH1D(limit_hist_bins);
        chain_from_model->FillHist( var->name, hist );
        limit->GetLimits(hist);

        if(save_hist){
          hist->SetName( (var->name + "_" + to_string(index) + "_" + to_string((int)chain_id)).c_str() );
          hist->Write();
        }
        delete hist;
      }
      output_data->Fill();
      delete chain_from_model;
      index++;
    }
  }

  TTree * generation_tree = new TTree("gen_description", "gen_description");
  vector<string> vnames; 
  vector<string> pnames; 
  for(auto var : *variables)  vnames.push_back(var->name);
  for(auto par : *parameters) pnames.push_back(par->name);
  generation_tree->Branch("gen_var_names", &vnames);
  generation_tree->Branch("gen_par_names", &pnames);
  generation_tree->Fill();
  out_file->Write();
  out_file->Close();
  delete out_file;
}


static double GL_x, GL_t_m, GL_t_p, GL_P, GL_C;
double PDF_ANAL_minimizer2d(const double *xx){
  const Double_t T_p = xx[0];
  const Double_t T_m = xx[1];
  if( T_p < T_m ) return 1000000;
  return -1*TMath::Poisson(GL_t_p,T_p) * TMath::Poisson(GL_t_m,T_m) * TMath::Poisson(GL_x, GL_P*(GL_C + T_p-T_m));
}

double PDF_ANAL_minimizer3d(const double *xx){
  const Double_t T_p = xx[0];
  const Double_t T_m = xx[1];
  const Double_t P   = xx[2];
  if( T_p < T_m ) return 1000000;
  return -1*TMath::Poisson(GL_t_p,T_p) * TMath::Poisson(GL_t_m,T_m) * TMath::Poisson(GL_x, P*(GL_C + T_p-T_m));
}

void FitAndWriteANALITIC(std::vector<shta::MarkovChain *> * chains_from_true, shta::Function * model, string fname, int n_to_burn, int save_every_divisor, int model_chain_lenght, int max_number_of_points_to_read, int max_number_of_chains_to_read, shta::ProposalFunction * prop, int limit_hist_bins, ROOT::Math::Minimizer* minimum){
  shta::MetropolisHastings mh_model = shta::MetropolisHastings(model, prop);
  mh_model.verbose_lvl = shta::verbose::ERROR;

  TFile * out_file = new TFile(fname.c_str(), "RECREATE");
  TTree * output_data = new TTree("ttree","ttree");

  double weight;
  double chain_id;
  double minuit_status, T_m_best, T_p_best;
  output_data->Branch("chain_id",     &chain_id    );
  output_data->Branch("weight",       &weight      );

  output_data->Branch("minuit_status",  &minuit_status );
  output_data->Branch("T_m_best",       &T_m_best      );
  output_data->Branch("T_p_best",       &T_p_best      );

  auto fchain = chains_from_true->at(0);
  auto model_parameters = &(model->used_parameters);
  auto model_variables  = &(model->used_variables);
  int n_model_pars = model_parameters->size();
  int n_model_vars = model_variables->size();
  int n_vars = model_variables->size();

  std::vector<Limit*> limits;
  for(auto var : *model_variables)
    limits.push_back( new Limit(var->name, output_data) );

  int weight_index = chains_from_true->at(0)->GetLabelIndex("MCMC_weight");

  // unfortunatly we will extract parameters by hands so fix this in case of no_0_* changes
  int data_index_t_p = fchain->GetLabelIndex( "t_p" );
  int data_index_t_m = fchain->GetLabelIndex( "t_m" );
  int data_index_x   = fchain->GetLabelIndex( "x" );
  Variable * model_x       = model->GetVariable( "x" );
  Variable * model_dT_best = model->GetVariable( "dT_best" );
  double t_p, t_m;

  // MINUIT PART
  minimum->SetMaxFunctionCalls(2000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(100000);  // for GSL
  minimum->SetTolerance(0.00001);
  minimum->SetPrintLevel(0);

  //ROOT::Math::Functor f(&PDF_ANAL_minimizer3d,3); // TODO generation of function for ROOT MINUIT minimization in SHTA
  bool MINIMIZE_P = true;

  ROOT::Math::Functor f2d(&PDF_ANAL_minimizer2d,2);
  ROOT::Math::Functor f3d(&PDF_ANAL_minimizer3d,3);
  if(MINIMIZE_P) minimum->SetFunction( f3d );
  else           minimum->SetFunction( f2d );

  shta::msg("main loop ... ");
  shta::MarkovChain * chain_from_model;
  for(chain_id = 0; chain_id < chains_from_true->size(); chain_id++){
    if(chain_id >= max_number_of_chains_to_read) break;
    shta::msg("begin base chain", chain_id, "... ");
    std::vector <std::vector<double> * > * points = & (chains_from_true->at(chain_id)->chain);
    int chain_size = points->size();

    // 1. burn-in
    chains_from_true->at(chain_id)->Print(10);
    int index = chains_from_true->at(chain_id)->GetIndexAfterBurn(n_to_burn, weight);

    // 2. produce chain for every model
    while(index < chain_size and index < max_number_of_points_to_read){
      bool save_hist = false;
      if( (1+index) % save_every_divisor == 0) save_hist = true;
      std::vector<double> * point = points->at(index);
      weight = point->at(weight_index);

      model_x->value = point->at( data_index_x );
      t_p = point->at( data_index_t_p );
      t_m = point->at( data_index_t_m );
      /*  
          For ANALITIC pdf we performe minimization with ROOT TMinuit at every base chain point to get dT_best
          for basic function:
          -TMath::Poisson(x, P*(4 + T_p-T_m)) * TMath::Poisson(t_p,T_p) * TMath::Poisson(t_m,T_m) * shta::Heaviside(T_p, T_m)
          we can performe minimization over (1) P,T_p,T_m or over (2) T_p,T_m with some prior know value of P
          (2) is a case when we try to set a better limits on already very know values
          and better suite to our study, so we will use (2) strategy and even with fails of minimization it will lead to better P estimation
      */

      // FIXME find better way to pass parameters
      GL_x   = model_x->value;
      GL_t_m = t_m;
      GL_t_p = t_p;
      GL_P   = 3;
      GL_C   = 4;

      minimum->SetVariable(0, "T_p", t_p, 0.01);
      minimum->SetVariable(1, "T_m", t_m, 0.01);
      if(MINIMIZE_P) minimum->SetVariable(2, "P", 3., 0.01);

      minimum->Minimize();

      const double *xs = minimum->X();
      //std::cout << "  x = " << GL_x << ", t_p = " << t_p << ", t_m = " << t_m << endl;
      //std::cout << "T_p = " << xs[0] << ", T_m = " << xs[1] << ", P = " << xs[2] << ", status = " << minimum->Status() <<std::endl;

      minuit_status = minimum->Status();
      T_p_best      = xs[0];
      T_m_best      = xs[1];
      model_dT_best->value = xs[0] - xs[1];

      out_file->cd();
      chain_from_model = mh_model.ConstructChain(model_chain_lenght);
      for(int i = n_vars-1; i >= 0 ; i--){
        Limit * limit = limits.at(i);
        Variable * var = model_variables->at(i);
        TH1D * hist = var->GetTH1D(limit_hist_bins);
        chain_from_model->FillHist( var->name, hist );
        limit->GetLimits(hist);

        if(save_hist){
          hist->SetName( (var->name + "_" + to_string(index) + "_" + to_string((int)chain_id)).c_str() );
          hist->Write();
        }
        delete hist;
      }
      output_data->Fill();
      delete chain_from_model;
      index++;
    }
  }

  TTree * generation_tree = new TTree("gen_description", "gen_description");
  vector<string> vnames; 
  vector<string> pnames; 
  for(auto var : *model_variables)  vnames.push_back(var->name);
  for(auto par : *model_parameters) pnames.push_back(par->name);
  generation_tree->Branch("gen_var_names", &vnames);
  generation_tree->Branch("gen_par_names", &pnames);
  generation_tree->Fill();

  out_file->Write();
  out_file->Close();
  delete out_file;
}

int main(int argc, char** argv){
  /*
    We will fit here every MC point generated by MCMC
    As long as number of resulting chains are too big to be saved so we will store only limits numbers.
    But in addition we will store some chains histograms for validation.
  */

  //ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent");
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit", "Minimize");
  
  string mode = argv[1];
  msg("will use mode = ", mode);

  shta::FunctionBuilder fbuilder;
  auto x   = fbuilder.DefineVariable( "x"  ,  5, 0., 100,  shta::vartype::PARAMETER );
  auto t_p = fbuilder.DefineVariable( "t_p",  5, 0., 100,  shta::vartype::PARAMETER );
  auto t_m = fbuilder.DefineVariable( "t_m",  5, 0., 100,  shta::vartype::PARAMETER );
  auto T_p = fbuilder.DefineVariable( "T_p",  8, 0., 100,  shta::vartype::VARIABLE );
  auto T_m = fbuilder.DefineVariable( "T_m",  6, 0., 100,  shta::vartype::VARIABLE );
  auto P   = fbuilder.DefineVariable( "P",    2, 0., 100,  shta::vartype::VARIABLE );
  auto T   = fbuilder.DefineVariable( "T",    2, 0., 100,  shta::vartype::VARIABLE );
  auto dT_best = fbuilder.DefineVariable( "dT_best", 2, 0., 100, shta::vartype::PARAMETER );
  auto PDF_NAIVE = fbuilder.DefineFunction("pdf_naive", "TMath::Poisson(x, P*(4 + t_p-t_m))");
  auto PDF_TWOPOISSON = fbuilder.DefineFunction("pdf_twopoisson", "TMath::Poisson(x, P*(4 + T_p-T_m)) * TMath::Poisson(t_p,T_p) * TMath::Poisson(t_m,T_m) * shta::Heaviside(T_p, T_m)");
  auto PDF_GAUS = fbuilder.DefineFunction("pdf_gaus", "TMath::Poisson(x, P*(4 + T)) * TMath::Gaus(t_p-t_m, T, sqrt(t_p+t_m)) * shta::Heaviside(T)");
  auto PDF_ANALITIC = fbuilder.DefineFunction("pdf_anal", "TMath::Poisson(x, P*(4 + dT_best))");

  shta::ProposalFunction * prop;

  msg(mode, ":", "build model ... ");
  fbuilder.BuildFunctions("negative_scripts_fit" + mode);

  msg(mode, ":", "read mc data ... ");
  auto chains_A = GetBaseMC("ns_0_PDF_BASE_a.root");
  auto chains_B = GetBaseMC("ns_0_PDF_BASE_b.root");

  msg(mode, ":", "fit mc data and save resutls ... ");


  int model_chain_lenght = 100000;
  int max_number_of_points_to_read = 10; // 50000
  int max_number_of_chains_to_read = 2; // 50000
  int burn_number = 0;
  int save_hist_every = 10;
  int limit_hist_bins = 100000; // = 100 / sigma_err

  model_chain_lenght = 150000;
  max_number_of_points_to_read = 10500; // 50000
  max_number_of_chains_to_read = 25; // 50000
  burn_number = 500;
  save_hist_every = 1500;
  limit_hist_bins = 200000;

  msg(mode, ":", "process model ``PDF_NAIVE`` ... ");
  shta::ProposalGaussDiffusion * prop_naive = new shta::ProposalGaussDiffusion();
  prop_naive->widths.push_back( 1. );
  prop = dynamic_cast<shta::ProposalFunction *>(prop_naive);
  if(mode == "na" or mode == "all") FitAndWrite(chains_A, PDF_NAIVE, "ns_1_PDF_NAIVE_a.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins);
  if(mode == "nb" or mode == "all") FitAndWrite(chains_B, PDF_NAIVE, "ns_1_PDF_NAIVE_b.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins);

  msg(mode, ":", "process model ``PDF_POISSON`` ... ");
  shta::ProposalGaussDiffusion * prop_poisson = new shta::ProposalGaussDiffusion();
  prop_poisson->widths.push_back( 1.0 );
  prop_poisson->widths.push_back( 1.2 );
  prop_poisson->widths.push_back( 0.9 );
  prop = dynamic_cast<shta::ProposalFunction *>(prop_poisson);
  if(mode == "pa" or mode == "all") FitAndWrite(chains_A, PDF_TWOPOISSON, "ns_1_PDF_TWOPOISSON_a.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins);
  if(mode == "pb" or mode == "all") FitAndWrite(chains_B, PDF_TWOPOISSON, "ns_1_PDF_TWOPOISSON_b.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins);

  msg(mode, ":", "process model ``PDF_GAUSS`` ... ");
  shta::ProposalGaussDiffusion * prop_gauss = new shta::ProposalGaussDiffusion();
  prop_gauss->widths.push_back( 1.0 );
  prop_gauss->widths.push_back( 1.5 );
  prop = dynamic_cast<shta::ProposalFunction *>(prop_gauss);
  if(mode == "ga" or mode == "all") FitAndWrite(chains_A, PDF_GAUS, "ns_1_PDF_GAUS_a.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins);
  if(mode == "gb" or mode == "all") FitAndWrite(chains_B, PDF_GAUS, "ns_1_PDF_GAUS_b.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins);

  msg(mode, ":", "process model ``PDF_ANALITIC`` ... ");

  shta::ProposalGaussDiffusion * prop_anal = new shta::ProposalGaussDiffusion();
  prop_anal->widths.push_back( 1.0 );
  prop = dynamic_cast<shta::ProposalFunction *>(prop_anal);
  if(mode == "aa" or mode == "all") FitAndWriteANALITIC(chains_A, PDF_ANALITIC, "ns_1_PDF_ANALITIC_a.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins, minimum);
  if(mode == "ab" or mode == "all") FitAndWriteANALITIC(chains_B, PDF_ANALITIC, "ns_1_PDF_ANALITIC_b.root", burn_number, save_hist_every, model_chain_lenght, max_number_of_points_to_read, max_number_of_chains_to_read, prop, limit_hist_bins, minimum);

  msg(mode, ":", "done ... ");

  delete prop_naive;
  delete prop_poisson;
  delete prop_gauss;
  delete prop_anal;
  delete minimum;

  msg(mode, ": happy exit ... ");

  return 0;
}









