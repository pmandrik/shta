
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#include "shta_core.hh"

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

void FitAndWrite(std::vector<shta::MarkovChain *> * chains_from_true, shta::Function * model, string fname, int model_chain_lenght, int max_number_of_points_to_read, int number_of_chains, shta::ProposalFunction * prop, int limit_hist_bins, vector<int> data_indexes){
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
  int data_index, try_points;
  output_data->Branch("chain_id",   &chain_id   );
  output_data->Branch("weight",     &weight     );
  output_data->Branch("data_index", &data_index );
  output_data->Branch("try_points", &try_points );

  shta::msg("main loop ... ");
  msg_progress( 0. );
  shta::MarkovChain * chain_from_model;
  bool save_hist = true;
  for(chain_id = 0; chain_id < chains_from_true->size(); chain_id++){
    if(chain_id) break;
    shta::msg("begin base chain", chain_id, "... ");
    std::vector <std::vector<double> * > * points = & (chains_from_true->at(chain_id)->chain);
    int chain_size = points->size();

    // 2. produce chain for every model
    try_points = 0;
    while(try_points < max_number_of_points_to_read){
      msg_progress( double(try_points)/max_number_of_points_to_read );
      data_index = data_indexes.at( try_points );

      std::vector<double> * point = points->at(data_index);
      weight = point->at(weight_index);

      for(int i = n_pars-1; i >= 0 ; i--)
        parameters->at(i)->value = point->at( indexes[i] );

      save_hist = true;
      for(int mc_index = 0; mc_index < number_of_chains; mc_index++){
        chain_from_model = mh_model.ConstructChain(model_chain_lenght);
        for(int i = n_vars-1; i >= 0 ; i--){
          Limit * limit = limits.at(i);
          Variable * var = variables->at(i);
          TH1D * hist = var->GetTH1D(limit_hist_bins);
          chain_from_model->FillHist( var->name, hist );
          limit->GetLimits(hist);

          if(save_hist){
            hist->SetName( (var->name + "_" + to_string(data_index) + "_" + to_string((int)chain_id)).c_str() );
            hist->Write();
            save_hist = false;
          }

          delete hist;
        }
        output_data->Fill();
        delete chain_from_model;
      }
      try_points++;
    }
  }
  msg_progress( 1. );
  TTree * generation_tree = new TTree("gen_description", "gen_description");
  vector<string> vnames; 
  vector<string> pnames; 
  for(auto var : *variables)  vnames.push_back(var->name);
  for(auto par : *parameters) pnames.push_back(par->name);
  generation_tree->Branch("gen_var_names", &vnames);
  generation_tree->Branch("gen_par_names", &pnames);
  generation_tree->Fill();

  out_file->Write();
}


int main(int argc, char** argv){
  /*
    To check error in pdf estimation due to MCMC method we will produce limits at the same data point
    with different random generator initial states for different number of chain lenght.
    As long as for ns_1_fit twopoison model has most number of parameters and
    worsest expected chain quality we will check only this model.
  */

  shta::FunctionBuilder fbuilder;
  auto x   = fbuilder.DefineVariable( "x"  ,  5, 0., 100,  shta::vartype::PARAMETER );
  auto t_p = fbuilder.DefineVariable( "t_p",  5, 0., 100,  shta::vartype::PARAMETER );
  auto t_m = fbuilder.DefineVariable( "t_m",  5, 0., 100,  shta::vartype::PARAMETER );
  auto T_p = fbuilder.DefineVariable( "T_p",  8, 0., 100,  shta::vartype::VARIABLE );
  auto T_m = fbuilder.DefineVariable( "T_m",  6, 0., 100,  shta::vartype::VARIABLE );
  auto P   = fbuilder.DefineVariable( "P",    2, 0., 100,  shta::vartype::VARIABLE );
  auto T   = fbuilder.DefineVariable( "T",    2, 0., 100,  shta::vartype::VARIABLE );
  auto PDF_TWOPOISSON = fbuilder.DefineFunction("pdf_twopoisson", "TMath::Poisson(x, P*(4 + T_p-T_m)) * TMath::Poisson(t_p,T_p) * TMath::Poisson(t_m,T_m) * shta::Heaviside(T_p, T_m)");

  shta::ProposalFunction * prop;

  msg("build model ... ");
  fbuilder.BuildFunctions("fir_mcmc_clq");

  msg("read mc data ... ");
  auto chains_A = GetBaseMC("ns_0_PDF_BASE_a.root");
  auto chains_B = GetBaseMC("ns_0_PDF_BASE_b.root");

  msg("fit mc data and save resutls ... ");
  std::vector<int> model_chain_lenghts = {10, 25, 50, 75, 100, 150, 250, 350, 500}; // * 1000
  int max_number_of_points_to_read = 50; // 50000
  int number_of_chains = 100; // 50000
  int limit_hist_bins = 100000;
  vector<int> data_indexes;
  for(int i = 0; i < max_number_of_points_to_read; i++) data_indexes.push_back( _random_engine.Integer( chains_A->at(0)->chain.size() ) );

  msg("process model ``PDF_POISSON`` ... ");
  shta::ProposalGaussDiffusion * prop_poisson = new shta::ProposalGaussDiffusion();
  prop_poisson->widths.push_back( 1.0 );
  prop_poisson->widths.push_back( 1.2 );
  prop_poisson->widths.push_back( 0.9 );
  prop = dynamic_cast<shta::ProposalFunction *>(prop_poisson);

  for(int mcl : model_chain_lenghts){
    int model_chain_lenght = mcl * 1000;
    msg("begin chain ... ", model_chain_lenght);
    FitAndWrite(chains_A, PDF_TWOPOISSON, "ns_3_PDF_TWOPOISSON" + to_string(model_chain_lenght) + "_a.root", 
                model_chain_lenght, max_number_of_points_to_read, number_of_chains, prop, limit_hist_bins, data_indexes);
    //FitAndWrite(chains_B, PDF_TWOPOISSON, "ns_1_PDF_TWOPOISSON" + to_string(model_chain_lenght) + "_b.root",
    //          model_chain_lenght, max_number_of_points_to_read, number_of_chains, prop, limit_hist_bins, data_indexes);
  }

  msg("done ... ");
  return 0;
}









