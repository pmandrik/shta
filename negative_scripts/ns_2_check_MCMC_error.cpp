
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

struct gen_data{
  map<string, double> parvals;
  int n_chains, lenght;
};

gen_data GetBaseMCparameters(string file_name){
  TFile * in_file = new TFile(file_name.c_str(), "READ");
  if(not in_file->IsOpen()) {
    msg("cant open file", file_name);
    return gen_data();
  }

  TTree * generation_tree = (TTree*)in_file->Get("gen_description");
  std::vector<std::string> names; 
  std::vector<double> values;
  auto * names_ptr  = &names;
  auto * values_ptr = &values;

  gen_data gdata;

  generation_tree->SetBranchAddress("gen_par_names", &names_ptr);
  generation_tree->SetBranchAddress("gen_par_values", &values_ptr);
  generation_tree->SetBranchAddress("n_chains", &(gdata.n_chains));
  generation_tree->SetBranchAddress("lenght",   &(gdata.lenght));

  generation_tree->GetEntry(0);

  msg("read data model with (", gdata.n_chains, "x", gdata.lenght, ") points and parameters:");
  for(int i = 0; i < names.size(); i++){
    msg(names.at(i), values.at(i));
    gdata.parvals[ names.at(i) ] = values.at(i);
  }
  in_file->Close();
  return gdata;
};

double get_diff(double val, double c){
  return (val - c) / c;
}

void GetPlots(gen_data gdata, const char * inpname, const char * outname, vector <double> shifts){
  TFile * in_file = new TFile(inpname, "READ");
  if(not in_file->IsOpen()) {
    msg("cant open file", inpname);
    return;
  }
  msg("process file ", inpname, "...");

  int n_chains = gdata.n_chains;
  vector<string> postfixes = {"_l_1s", "_c_1s", "_u_1s", "_l_2s", "_c_2s", "_u_2s"};

  // read "gen_description" -----------------------
  TTree * generation_tree = (TTree*)in_file->Get("gen_description");
  std::vector<std::string> pnames; 
  std::vector<std::string> vnames; 
  auto * pnames_ptr  = &pnames;
  auto * vnames_ptr  = &vnames;

  generation_tree->SetBranchAddress("gen_par_names", &pnames_ptr);
  generation_tree->SetBranchAddress("gen_var_names", &vnames_ptr);
  generation_tree->GetEntry(0);

  // setup read "ttree" -----------------------
  TTree * mcmc_tree = (TTree*)in_file->Get("ttree");
  //mcmc_tree->Print();

  std::vector<Limit*> limits;
  //msg( "used paramters : " );
  for(auto name : vnames){
    //msg( name );
    Limit * lim = new Limit(name, mcmc_tree, false);
    lim->def_val = (*(gdata.parvals.find( name ))).second;
    limits.push_back( lim );
  }
  double weight;
  double chain_id;
  mcmc_tree->SetBranchAddress("chain_id",     &chain_id    );
  mcmc_tree->SetBranchAddress("weight",       &weight      );
  Long64_t nentries = mcmc_tree->GetEntries();

  // produce images part -----------------------
  TFile * out_hists = new TFile( outname, "RECREATE");

  // distributions of limits
  map <string, THStack *> dlimits;
  map <string, TH1D *> dhists;
  for(auto name : vnames)
    for(auto postfix : postfixes){
      dlimits[name+postfix] = new THStack( (name+postfix).c_str(), (name+postfix).c_str() );
      for(int n = 0; n < n_chains; n++){
        auto hists = new TH1D(root::get_name().c_str(), "", 100, -10, 10);
        hists->SetFillColor(kOrange+n);
        hists->SetLineColor(kOrange+n);
        dlimits[name+postfix]->Add( hists );
        dhists[name+postfix + to_string(n)] = hists;
      }
    }

  for(Long64_t i=0; i<nentries; i++){
    mcmc_tree->GetEntry(i);
    for(auto limit : limits){
      if(TMath::Abs(limit->l_2s - limit->u_2s) < 0.0001 ) continue;
      for(auto postfix : postfixes)
        (*(dhists.find( limit->name + postfix + to_string((int)chain_id)   ))).second->Fill( get_diff(limit->GetValue(limit->name + postfix), limit->def_val), weight);
    }
  }
  for(auto pair : dlimits){
    root::get_canvas();
    pair.second->Draw("hist f");
    pair.second->Write();
  }

  double shift_l_1s = shifts[0];
  double shift_u_1s = shifts[1];
  double shift_l_2s = shifts[2];
  double shift_u_2s = shifts[3];

  // distributions 
  for(auto limit : limits){
    vector <double> find_n_1s = vector<double>(n_chains, 0.);
    vector <double> find_n_2s = vector<double>(n_chains, 0.);
    vector <double> find_n_1s_shift = vector<double>(n_chains, 0.);
    vector <double> find_n_2s_shift = vector<double>(n_chains, 0.);
    vector <double> tot_n = vector<double>(n_chains, 0.);
    double entries_sum = 0.;
    for(Long64_t i=0; i<nentries; i++){
      mcmc_tree->GetEntry(i);

      tot_n[ (int)chain_id ] += weight;
      if(TMath::Abs(limit->l_2s - limit->u_2s) < 0.0001 ) continue;

      if(limit->l_1s <= limit->def_val and limit->def_val <= limit->u_1s)
        find_n_1s[ (int)chain_id ] += weight;
      if(limit->l_2s <= limit->def_val and limit->def_val <= limit->u_2s)
        find_n_2s[ (int)chain_id ] += weight;
      if(limit->l_1s + shift_l_1s <= limit->def_val and limit->def_val <= limit->u_1s - shift_u_1s)
        find_n_1s_shift[ (int)chain_id ] += weight;
      if(limit->l_2s + shift_l_2s <= limit->def_val and limit->def_val <= limit->u_2s - shift_u_2s)
        find_n_2s_shift[ (int)chain_id ] += weight;
    }

    TH1D * frac_1s = new TH1D( (limit->name+" in frac 1#sigma").c_str() , (limit->name+" in frac 1#sigma").c_str(), 20000, 0.0, 1.0);
    TH1D * frac_2s = new TH1D( (limit->name+" in frac 2#sigma").c_str() , (limit->name+" in frac 2#sigma").c_str(), 20000, 0.0, 1.0);
    TH1D * frac_1s_shift = new TH1D( (limit->name+" in frac 1#sigma shifted limits").c_str() , (limit->name+" in frac 1#sigma shifted limits").c_str(), 20000, 0.0, 1.0);
    TH1D * frac_2s_shift = new TH1D( (limit->name+" in frac 2#sigma shifted limits").c_str() , (limit->name+" in frac 2#sigma shifted limits").c_str(), 20000, 0.0, 1.0);
    if(limit->name == "P") msg( "Fracs for ", limit->name, " : " );
    for(int i = 0; i < (int)chain_id; i++){
      frac_1s->Fill( find_n_1s[ i ] / tot_n[ i ] );
      frac_2s->Fill( find_n_2s[ i ] / tot_n[ i ] );
      frac_1s_shift->Fill( find_n_1s_shift[ i ] / tot_n[ i ] );
      frac_2s_shift->Fill( find_n_2s_shift[ i ] / tot_n[ i ] );
      if(limit->name == "P" and false) msg("nominal = ", find_n_1s[ i ] / tot_n[ i ], find_n_2s[ i ] / tot_n[ i ], "| shifted = ", find_n_1s_shift[ i ] / tot_n[ i ], find_n_2s_shift[ i ] / tot_n[ i ] );
    }
    if(limit->name == "P"){
      msg( "nominal 1s mean/rms = ", frac_1s->GetMean(), " / ", frac_1s->GetRMS() );
      msg( "nominal 2s mean/rms = ", frac_2s->GetMean(), " / ", frac_2s->GetRMS() );
      msg( "shifted 1s mean/rms = ", frac_1s_shift->GetMean(), " / ", frac_1s_shift->GetRMS() );
      msg( "shifted 2s mean/rms = ", frac_2s_shift->GetMean(), " / ", frac_2s_shift->GetRMS() );
      double err_1s_a = TMath::Abs(frac_1s->GetMean() - frac_1s_shift->GetMean());
      double err_1s_b = frac_1s->GetRMS();
      msg("fnal 1s = ", frac_1s->GetMean(), "+/-", err_1s_a, "+/-", err_1s_b);
      double err_2s_a = TMath::Abs(frac_2s->GetMean() - frac_2s_shift->GetMean());
      double err_2s_b = frac_2s->GetRMS();
      msg("fnal 2s = ", frac_2s->GetMean(), "+/-", err_2s_a, "+/-", err_2s_b);
    }
  
    frac_1s->Write();
    frac_2s->Write();
  }

  out_hists->Close();
};

int main(int argc, char** argv){
  /* 
    In this script we will read limits and compare with true ones.
    As a measure ... TODO
  */
  //TApplication theApp("App",&argc, argv);

  auto gdata_a = GetBaseMCparameters("ns_0_PDF_BASE_a.root");
  auto gdata_b = GetBaseMCparameters("ns_0_PDF_BASE_b.root");

                          /*l1s   u1s    l2s   u2s*/
  vector <double> shifts = {0.01, 0.02, 0.01,  0.04};

  GetPlots(gdata_a, "ns_1_PDF_NAIVE_a.root", "ns_2_PDF_NAIVE_a.root",           shifts );
  GetPlots(gdata_a, "ns_1_PDF_TWOPOISSON_a.root", "ns_2_PDF_TWOPOISSON_a.root", shifts );
  GetPlots(gdata_a, "ns_1_PDF_GAUS_a.root", "ns_2_PDF_GAUS_a.root",             shifts );
  GetPlots(gdata_a, "ns_1_PDF_ANALITIC_a.root", "ns_2_PDF_ANALITIC_a.root",     shifts );

  GetPlots(gdata_b, "ns_1_PDF_NAIVE_b.root", "ns_2_PDF_NAIVE_b.root",           shifts );
  GetPlots(gdata_b, "ns_1_PDF_TWOPOISSON_b.root", "ns_2_PDF_TWOPOISSON_b.root", shifts );
  GetPlots(gdata_b, "ns_1_PDF_GAUS_b.root", "ns_2_PDF_GAUS_b.root",             shifts );
  //GetPlots(gdata_b, "ns_1_PDF_ANALITIC_b.root", "ns_2_PDF_ANALITIC_b.root",     shifts );

  //shta::msg("Press Ctrl+C to exit root TApplication ...");
  //theApp.Run();

  return 0;
}









