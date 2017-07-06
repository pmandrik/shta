
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_MCMC_HH
#define SHTA_MCMC_HH 1

namespace shta {
  class ProposalFunction{
    public:
      ProposalFunction(){}
      virtual void Eval( std::vector<double> & point_prev, std::vector<double> & point_next ) = 0;
      std::vector <double> var_range_mins;
      std::vector <double> var_range_maxs;
      std::vector <double> var_range_deltas;
      int n_variables, n_parameters;
  };

  class ProposalUniform : public ProposalFunction {
    public:
      void Eval( std::vector<double> & point_prev, std::vector<double> & point_next ){
        for(int i = 0; i < n_variables; i++){
          point_next[i] = _random_engine.Uniform( var_range_mins[i], var_range_maxs[i] );
          //msg( var_range_mins[i], var_range_maxs[i], point_next[i] );
        }
      }
  };

  class ProposalGauss : public ProposalFunction {
    public:
      void Eval( std::vector<double> & point_prev, std::vector<double> & point_next ){
        double point, rmax, rmin, delta;
        for(int i = 0; i < n_variables; i++){
          point = point_prev[i];
          rmin = var_range_mins[i];
          rmax = var_range_maxs[i];
          delta = var_range_deltas[i];
          point += delta * step_factor * _random_engine.Gaus();
          while(point < rmin) point += delta;
          while(point > rmax) point -= delta;
          point_next[i] = point;
        }
      }
      double step_factor;
  };

  class ProposalGaussDiffusion : public ProposalFunction {
    public:
      void Eval( std::vector<double> & point_prev, std::vector<double> & point_next ){
        for(int i = 0; i < n_variables; i++)
          point_next[i] = point_prev[i] + widths[i] * _random_engine.Gaus();
      }
      
      double TuneWidth(double width, int n_pars){
        return 2.38 * width / TMath::Sqrt( n_pars );
      }

      std::vector <double> widths;
  };

  class MarkovChain : public BaseClass {
    public:
      MarkovChain(){
        n_variables = 0;
        weighted_lenght = -1;
      }

      MarkovChain(int n_vars){
        n_variables = n_vars;
        labels.reserve( n_variables );
      }

      ~MarkovChain(){
        for(std::vector <std::vector<double> * >::iterator iter = chain.begin(); iter != chain.end(); ++iter) 
          delete *iter;
      }

      void Write(std::string folder_name) {
        TDirectory * mc_dir = gDirectory->mkdir(folder_name.c_str());
        mc_dir->cd();

        TTree * labels_tree = new TTree("labels", "labels");
        labels_tree->Branch("labels", &labels);
        labels_tree->Fill();
        labels_tree->Write();

        TTree * chain_tree = new TTree("chain", "chain");
        std::vector<double> * mc_vec;
        chain_tree->Branch("chain", &mc_vec);
        for(auto point : chain){
          mc_vec = point;
          chain_tree->Fill();
        }
        chain_tree->Write();

        delete labels_tree;
        delete chain_tree;
      }

      bool Read(std::string folder_name) {
        if(gDirectory->cd(folder_name.c_str()) != kTRUE){
          msg_err("WARNING : shta : class MarkovChain : Read() : can't find Markov Chain folder \"", folder_name, "\", return false");
          return false;
        }

        TTree * labels_tree = (TTree*)gDirectory->Get("labels");
        TTree * chain_tree  = (TTree*)gDirectory->Get("chain");
        if(not labels_tree or not chain_tree){
          msg_err("WARNING : shta : class MarkovChain : Read() : can't find \"label\" or \"chain\" Tree in folder, return false");
          return false;
        }
        if(not labels_tree->GetBranch("labels") or not chain_tree->GetBranch("chain")){
          msg_err("WARNING : shta : class MarkovChain : Read() : can't find \"label\" or \"chain\" branch in Tree, return false");
          return false;
        }
        auto * labels_ptr = &labels;
        labels_tree->SetBranchAddress("labels", &labels_ptr);
        labels_tree->GetEntry(0);
        n_variables = labels.size();

        std::vector<double> * mc_vec = 0;
        chain_tree->SetBranchAddress("chain", &mc_vec);
        Long64_t nentries = chain_tree->GetEntries();
        for (Long64_t i=0; i<nentries; i++) {
          chain_tree->GetEntry(i);
          std::vector<double> * new_vec = new std::vector<double>( *mc_vec );
          chain.emplace_back( new_vec );
        }
        delete labels_tree;
        delete chain_tree;
        return true;
      }

      int GetLabelIndex(std::string label) const {
        unsigned int index = 0;
        for(index = 0; index < labels.size(); index++) if(labels[index] == label) break;
        if (index == labels.size()){
          MSG_WARNING("WARNING : shta : class MarkovChain : GetLabelIndex() : cant find label \"", label, "\", return -1");
          return -1;
        }
        return index;
      }

      MarkovChain * GetSlice(std::string label, double min_val, double max_val) const {
        int index = GetLabelIndex( label );
        if(index == -1) return nullptr;

        MarkovChain * slice = new MarkovChain( n_variables );
        slice->labels = labels;
        double val;

        for( auto point : chain ) {
          val = (*point)[index];
          if(val < min_val or val > max_val) continue;
          slice->chain.push_back( point );
        }
        return slice;
      }

      MarkovChain * GetSlice( const std::vector<unsigned int> & indexes, const std::vector<double> & min_vals, const std::vector<double> & max_vals ) const {
        MarkovChain * slice = new MarkovChain( n_variables );
        slice->labels = labels;
        double val;
        for( auto point : chain ) {
          for(int i = indexes.size() - 1; i >= 0; i--){
            val = (*point)[indexes[i]];
            if(val < min_vals[i] or val > max_vals[i]) goto sKIP_THIS_POINT;
          }
          slice->chain.push_back( point );
          sKIP_THIS_POINT:;
        }
        return slice;
      }

      MarkovChain * GetSlice( const std::vector<std::string> & labels, const std::vector<double> & min_vals, const std::vector<double> & max_vals ) const {
        std::vector<unsigned int> indexes;
        for(auto label : labels) indexes.push_back( GetLabelIndex( label ) );
        return GetSlice(indexes, min_vals, max_vals);
      }

      double GetChainWeightedLenght() const { return weighted_lenght; }
      double FindChainWeightedLenght() { 
        weighted_lenght = 0;
        unsigned int chain_index_weight = labels.size()-2;
        for(auto point : chain) weighted_lenght += (*point)[chain_index_weight];
      }

      void FillHist(std::string label, TH1D * hist) const {
        int index = GetLabelIndex( label );
        if(index == -1) return;

        unsigned int chain_index_weight = labels.size()-2;
        for(auto point : chain){
          hist->Fill( (*point)[index], (*point)[chain_index_weight] );
        }
      }

      void Fill2DHist(std::string label_x, std::string label_y, TH2D * hist) const {
        int index_x = GetLabelIndex( label_x );
        int index_y = GetLabelIndex( label_y );
        if(index_x == -1 or index_y == -1) return;

        unsigned int chain_index_weight = labels.size()-2;
        for(auto point : chain){
          hist->Fill( (*point)[index_x], (*point)[index_y], (*point)[chain_index_weight] );
        }
      }

      long long int GetIndexAfterBurn(int n_to_burn, double & last_weight){
        int weight_index = GetLabelIndex( "MCMC_weight" );
        int sum_weight = 0;
        long long int index = 0;
        long long int size = chain.size();
        if(n_to_burn > size) return size-1;
        while(index < size){
          sum_weight += chain.at(index)->at(weight_index);
          if(sum_weight > n_to_burn){
            last_weight = sum_weight - n_to_burn;
            break;
          }
          index++;
        }
        return index;
      }

      TH1D * DrawVariable(const Variable * var, int nbins){
        TH1D * hist = var->GetTH1D(nbins);
        FillHist(var->name, hist);
        root::get_canvas();
        hist->Draw();
        return hist;
      }

      void Print(unsigned long long int deep = 0){
        std::cout << "shta::MarkovChain [";
        for(std::string label : labels) std::cout << label << ", ";
        std::cout << "] " << chain.size()  << std::endl;

        deep = std::min(deep, (unsigned long long int )chain.size());
        for(unsigned long long int i = 0; i < deep; i++){
          std::vector<double> * point = chain[i];
          std::cout << i << " | ";
          for(double val : *point) std::cout << val << ", ";
          std::cout << std::endl;
        }
      }

      int n_variables;
      std::vector <std::string> labels;
      std::vector <std::vector<double> * > chain;
      double weighted_lenght;
  };

  class MarkovChainMonteCarlo : public BaseClass {
    public:
      MarkovChainMonteCarlo(Function * func, ProposalFunction * prop) : pdf(func), proposal(prop) {
        n_variables  = func->used_variables.size();
        n_parameters = func->used_parameters.size();

        proposal->var_range_mins.reserve( n_variables );
        proposal->var_range_maxs.reserve( n_variables );
        proposal->var_range_deltas.reserve( n_variables );

        proposal->n_variables  = n_variables;
        proposal->n_parameters = n_parameters;

        pdf_eval = pdf->eval;

        for( auto var : func->used_variables){
          proposal->var_range_mins.push_back( var->range_min );
          proposal->var_range_maxs.push_back( var->range_max );
          proposal->var_range_deltas.push_back( var->range_max-var->range_min );
        }
      }

      MarkovChain * PrepareNewChain(){
        MarkovChain * chain = new MarkovChain( n_variables+2 );
        for(const Variable * var : pdf->used_variables ) chain->labels.push_back( var->name );
        chain->labels.push_back( "MCMC_weight" );
        chain->labels.push_back( "MCMC_pdf" );
        chain->verbose_lvl = verbose_lvl;
        return chain;
      }

      virtual MarkovChain * ConstructChain(int n_iterations) = 0;

      int n_variables, n_parameters;
      Function * pdf;
      ProposalFunction * proposal;
      double (*pdf_eval)(std::vector<double> & x, std::vector<double> & p);
  };

  class MetropolisHastings : public MarkovChainMonteCarlo{
    public:
      MetropolisHastings(Function * func, ProposalFunction * prop) : MarkovChainMonteCarlo(func, prop) {}

      MarkovChain * ConstructChain(int n_iterations){
        MarkovChain * chain = PrepareNewChain();
        std::vector < std::vector<double> * > * chain_v = & chain->chain;
        int chain_width = chain->n_variables;
        int chain_index_weight = chain_width-2, chain_index_pdf = chain_width-1;

        std::vector<double> * point_prev = new std::vector<double>();
        std::vector<double> * point_next = new std::vector<double>( chain_width,  0.0 );
        std::vector<double> * parameters = new std::vector<double>();
        parameters->reserve( n_parameters );
        for(const Variable * variable : pdf->used_variables) point_prev->push_back( variable->value );
        for(const Variable * parameter : pdf->used_parameters) parameters->push_back( parameter->value );

        //for(const Variable * parameter : pdf->used_parameters) msg( parameter->value );

        //proposal->Eval( *point_prev, *point_next );
        double T, P_new, P_prev = pdf_eval( *point_prev, *parameters );
        point_prev->push_back( 0 );      // "MCMC_weight"
        point_prev->push_back( P_prev ); // "MCMC_pdf"
        chain_v->push_back( point_prev );

        MSG_INFO("shta::MetropolisHastings begin chain construction ... ");
        int barWidth = 50;
        int progress_point = n_iterations - n_iterations / barWidth;
        if(verbose_lvl >= shta::verbose::INFO) msg_progress(0.);
        for(int i = n_iterations; i; i--){
          proposal->Eval( *point_prev, *point_next );
          P_new = pdf_eval( *point_next, *parameters );
          T = TMath::Min(1., P_new / P_prev);
          if( _random_engine.Uniform(0., 1) < T ){
            // add new point
            (*point_next)[chain_index_weight] = 1.;
            (*point_next)[chain_index_pdf]    = P_prev;
            chain_v->push_back( point_next );
            point_prev = point_next;
            P_prev     = P_new;
            point_next = new std::vector<double>( chain_width, 0.0 );

            // print progress sometimes
            if(verbose_lvl >= shta::verbose::INFO and i < progress_point){
              msg_progress( 1. - double(i)/n_iterations );
              progress_point = std::max(progress_point - n_iterations / barWidth, 1);
            }
            continue;
          }
          // add same point ~> increase weight
          (*point_prev)[chain_index_weight] += 1.;
        }
        if(verbose_lvl >= shta::verbose::INFO) msg_progress(1.);
        MSG_INFO("\nshta::MetropolisHastings end chain construction ... ");
        return chain;
      }

      MarkovChain * ConstructChainNLL(int n_iterations){
        MarkovChain * chain = PrepareNewChain();
        std::vector < std::vector<double> * > * chain_v = & chain->chain;
        int chain_width = chain->n_variables;
        int chain_index_weight = chain_width-2, chain_index_pdf = chain_width-1;

        std::vector<double> * point_prev = new std::vector<double>();
        std::vector<double> * point_next = new std::vector<double>( chain_width,  0.0 );
        std::vector<double> * parameters = new std::vector<double>();
        parameters->reserve( n_parameters );
        for(const Variable * variable : pdf->used_variables) point_prev->push_back( variable->value );
        for(const Variable * parameter : pdf->used_parameters) parameters->push_back( parameter->value );

        double T, P_new, P_prev = pdf_eval( *point_prev, *parameters );
  
        (*point_prev)[chain_index_weight] = 1.;
        (*point_prev)[chain_index_pdf] = P_prev;
        chain_v->push_back( point_prev );

        MSG_INFO("shta::MetropolisHastings begin chain construction ... ");
        int barWidth = 50;
        int progress_point = n_iterations - n_iterations / barWidth;
        if(verbose_lvl >= shta::verbose::INFO) msg_progress(0.);
        else progress_point = -1;
        for(int i = n_iterations; i; i--){
          proposal->Eval( *point_prev, *point_next );
          P_new = pdf_eval( *point_next, *parameters );
          if( P_new <= P_prev or _random_engine.Uniform(0., 1) < TMath::Exp(P_new - P_prev) ){
            // add new point
            (*point_next)[chain_index_weight] = 1.;
            (*point_next)[chain_index_pdf]    = P_prev;
            chain_v->push_back( point_next );
            point_prev = point_next;
            P_prev     = P_new;
            point_next = new std::vector<double>( chain_width, 0.0 );

            // print progress sometimes
            if(i < progress_point){
              msg_progress( 1. - double(i)/n_iterations );
              progress_point = std::max(progress_point - n_iterations / barWidth, 1);
            }
            continue;
          }
          // add same point ~> increase weight
          (*point_prev)[chain_index_weight] += 1.;
        }
        if(verbose_lvl >= shta::verbose::INFO) msg_progress(1.);
        MSG_INFO("\nshta::MetropolisHastings end chain construction ... ");
        return chain;
      }
  };
};

#endif 




