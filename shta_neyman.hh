
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_NEYMAN_HH
#define SHTA_NEYMAN_HH 1

/*
  Classes for Neyman Construction of confedences intervals
*/

namespace shta {

  class Interval {
    public:
      bool ContainDummy(){
        msg_err("WARNING : shta :: class Interval : ContainDummy() : Contain is not implemented for this Interval class, return false");
        return false;
      }
      virtual bool Contain( const double & point ){ return ContainDummy(); };
      virtual bool Contain( const std::vector <double> & point ){ return ContainDummy(); };
  };

  class Interval1D : public Interval {
    public:
      Interval1D(long long int lndx, long long int undx) : lower_index(lndx), upper_index(undx) {};
      using Interval::Contain;
      bool Contain(const double & point){
        if(point < lower_abs) return false;
        if(point > upper_abs) return false;
        return true;
      }

      long long int lower_index, upper_index;
      double lower_abs, upper_abs;
  };
  
  class ElipticInterval : public Interval {
    public:
      using Interval::Contain;
      bool Contain( const std::vector <double> & point ){
        int size = point.size();
        for(int i = 0; i < size; i++)
          if( abs(point[i] - center[i]) > area_ranges[i] ) return false;
        return true;
      }

      std::vector <double> center;
      std::vector <double> area_ranges;
  };

  class OrderingFunction{
    public:
      OrderingFunction() : ConfidenceLevel(0.95) {};

      virtual Interval1D * GetInterval( const std::vector <double> & slice ){};
      virtual Interval * FindArea(const MarkovChain * toys){};
      double ConfidenceLevel;
  };

  class OrderingFunction1D : public OrderingFunction {
    public:
      OrderingFunction1D() : intervalType(1), OrderingFunction() {
      }

      Interval1D * GetInterval(const std::vector <double> & slice){
        // central
        double CL_left  = (1.-ConfidenceLevel) * 0.5;
        double CL_right = (1.-ConfidenceLevel) * 0.5;
        if(intervalType == 2){
          CL_left = 0;
          CL_right = 1.-ConfidenceLevel;
        } else 
        if(intervalType == 0) {
          CL_left = 1.-ConfidenceLevel;
          CL_right = 0;
        }

        // sum
        double sum = 0;
        for(auto item : slice) sum += item;
        // lover 
        unsigned long long int l_index = 0;
        for(double lpoint = 0; l_index < slice.size(); l_index++){
          lpoint += slice[l_index];
          if(lpoint > CL_left * sum) break;
        }
        // upper
        unsigned long long int u_index = slice.size()-1;
        for(double rpoint = 0.; u_index > 0; u_index--){
          rpoint += slice[u_index];
          if(rpoint > CL_right * sum) break;
        }
        return new Interval1D(--l_index, ++u_index);
      }
      
      void SetTypeLower(){   intervalType = 0; }
      void SetTypeCentral(){ intervalType = 1; }
      void SetTypeUpper(){   intervalType = 2; }

      int intervalType;
  };
  
  class ElipticZeroOrderingFunction : public OrderingFunction {
    public:
      bool LessDistance(const std::pair<std::vector<double> *, double> & s1, const std::pair<std::vector<double> *, double> & s2){
        return s1.second < s2.second;
      }

      Interval * FindArea(const MarkovChain * toys){
        /* 
          This ordering function used to find minimal eliptical area around zero point \vec{x}=\vec{0}:
            1. Find ranges of \vec{x}  most of points ~ RMS from zero
            2. Sort points based on distance to \vec{0} (normalized on ranges from punkt 1.)
            3. Add points till we had CL points
        */

        // 1. Find ranges of \vec{x}  most of points ~ RMS from zero
        const std::vector <std::vector<double> * > * chain = & toys->chain;
        // if toys do not cover this area, then no information like uniform, return 
        if(chain->size() == 0){
          msg("WARNING : shta :: class ElipticZeroOrderingFunction : FindArea() : empty chain, return");
          return nullptr;
        }
        
        double integral = 0, sum = 0;
        unsigned int vec_x_dimension = toys->n_variables-2; // exclude weight and pdf
        unsigned int chain_index_weight = toys->n_variables-2; // weight index
        double chain_w_size = toys->GetChainWeightedLenght();
        std::vector<double> x_points_rms    = std::vector<double>( vec_x_dimension,  0.0 );

        for(auto point : *chain){
          integral += (*point)[chain_index_weight];
          for(unsigned int i = 0; i < vec_x_dimension; ++i)
            x_points_rms[i]  += ((*point)[i]) * ((*point)[i]); // rms from 0
        }
        for(unsigned int i = 0; i < vec_x_dimension; ++i) 
          x_points_rms[i] /= chain_w_size;

        // 2. Sort points based on distance to \vec{0} (normalized on ranges from punkt 1.)
        std::vector <std::pair<std::vector<double> *, double> > sorted_copy;
        sorted_copy.reserve( chain->size() );
        double distance;
        for(auto point : *chain){
          distance = 0;
          for(unsigned int i = 0; i < vec_x_dimension; ++i) distance += pow((*point)[i], 2) / x_points_rms[i] ;
          //distance += pow((*point)[i], 2);
          sorted_copy.push_back( std::make_pair(point, TMath::Sqrt( distance )) );
        }

        struct {
          bool operator()(const std::pair<std::vector<double> *, double> & s1, const std::pair<std::vector<double> *, double> & s2){   
            return s1.second < s2.second;
          }   
        } LessDistance;
        std::sort(sorted_copy.begin(), sorted_copy.end(), LessDistance);

        // 3. Add points till we had CL points
        unsigned long long int chain_size = sorted_copy.size();
        unsigned long long int index;
        integral *= this->ConfidenceLevel;
        for(index = 0; index < chain_size; ++index ){
          std::vector<double> * point = sorted_copy[index].first;
          sum += (*(sorted_copy[index].first))[chain_index_weight];
          if(sum >= integral) break;
        }

        // slice from 0 to index is a accepted area defined as point \vec{0} and ranges
        ElipticInterval * interval_A = new ElipticInterval();
        interval_A->center = std::vector<double>( vec_x_dimension,  0.0 );
        std::vector<double> * edge_point = sorted_copy[index].first;
        for(unsigned int i = 0; i < vec_x_dimension; ++i){
          interval_A->area_ranges.push_back( std::abs((*edge_point)[i]) );
          //msg("area #", i, std::abs((*edge_point)[i]));
        }

        return dynamic_cast<shta::Interval *>(interval_A);
      }
  };

  class NeymanConstruction{
    public:
      NeymanConstruction(){
        ordering_func = nullptr;
        mc_toys = nullptr;
      };

      std::vector< std::vector<double>* > * ConstructFromMarkovChain(MarkovChain * chain, const std::vector<double> & par_range_steps){
        if( ordering_func == nullptr  ){
          msg_err("WARNING : shta :: class NeymanConstruction : ConstructFromMarkovChain() : ordering_func == nullptr, return");
          return nullptr;
        }
        if( parameters.size() != par_range_steps.size() ){
          msg_err("WARNING : shta :: class NeymanConstruction : ConstructFromMarkovChain() : parameters.size() != par_range_steps.size(), return");
          return nullptr;
        }
        /*
          Neyman construction I(\vec{x}_{exp}) is a set of points \vec{p}: 
          \vec{x}_{exp} \in A(\vec{p}) => \vec{p} \in I(\vec{x}_{exp})

          So, we should know:
            1. which variables in markov chain is a parameter for \vec{p} and witch is a variable for \vec{x} <- parameters, variables
            2. steps in parameter space
            3. ordering function to construct A(\vec{p})
        */

        mc_toys = chain;
        std::forward_list <std::vector<double> * > chain_original_fl;
        for(auto point : mc_toys->chain) chain_original_fl.push_front( point );

        // map paramters
        unsigned int n_variables_original = mc_toys->n_variables;

        unsigned int N_parameters = parameters.size();
        unsigned int N_variables  = variables.size();

        std::vector<unsigned int> parameters_id;
        std::vector<unsigned int> variables_id;
        parameters_id.reserve( N_parameters );
        variables_id.reserve( N_variables );
        for(auto var : parameters) parameters_id.push_back( chain->GetLabelIndex(var->name) );
        for(auto var : variables)  variables_id.push_back( chain->GetLabelIndex(var->name) );

        // ranges of parameters
        std::vector<double> min_vals;
        std::vector<double> max_vals;
        std::vector<double> par_range_mins;
        std::vector<double> par_range_maxs;
        min_vals.reserve(  N_parameters );
        max_vals.reserve(  N_parameters );
        par_range_mins.reserve(  N_parameters );
        par_range_maxs.reserve(  N_parameters );

        for(unsigned int i = 0; i < N_parameters; i++){
          const Variable * par = parameters.at(i);
          par_range_mins.push_back( par->range_min );
          par_range_maxs.push_back( par->range_max );
          min_vals.push_back( par->range_min );
          max_vals.push_back( par->range_min + par_range_steps[i] );
        }

        // variables_exp
        std::vector<double> variables_exp_values;
        variables_exp_values.reserve( variables.size() );
        for(auto var : variables) variables_exp_values.push_back(var->value);

        // answer is a set_I - set of accepted \vec{p} values
        set_I = std::vector<std::vector<double>* >();
        unsigned long long int set_I_size = 1;
        for(unsigned int i = 0; i < N_parameters; i++) set_I_size *= (par_range_maxs[i] - par_range_mins[i])/par_range_steps[i];
        set_I.reserve( set_I_size );
        std::vector<double> * parameter_point;

        MarkovChain * slice = new MarkovChain( N_variables+2 );
        auto slice_chain = & slice->chain;
        std::vector<double> * new_point;
        for(auto var : variables) slice->labels.push_back(var->name);
        slice->labels.push_back( "MCMC_weight" );
        slice->labels.push_back( "MCMC_pdf" );

        Interval * area_A;
        unsigned int i = 0;

        // iterate over parameters ranges
        msg("shta::NeymanConstruction begin Neyman construction ... ");
        unsigned long long int n_points = 0;
        int barWidth = 50;
        unsigned long long int progress_point = set_I_size / barWidth;
        msg_progress(0.);

        double val;
        std::vector <std::vector<double> * >::iterator iter_slice;
        std::forward_list <std::vector<double> * >::iterator iter_before;
        std::forward_list <std::vector<double> * >::iterator iter;
        
        while(true){
          // for this \vec{p}' construct set  A(\vec{p}') of \vec{x}
          //slice = chain->GetSlice( parameters_id, min_vals, max_vals ); // <- this is a P(\vec{x}|\vec{p}')
          //delete slice;
          for(iter_slice = slice_chain->begin(); iter_slice != slice_chain->end(); ++iter_slice)
            delete *iter_slice;
          slice_chain->clear();

          for( iter_before = chain_original_fl.before_begin(); ; ++iter_before ) {
            iter = std::next( iter_before );
            if(iter == chain_original_fl.end()) break;
            for(int i = parameters_id.size() - 1; i >= 0; i--){
              val = (**iter)[parameters_id[i]];
              if(val < min_vals[i] or val > max_vals[i]) goto sKIP_THIS_POINT;
            }
            new_point = new std::vector<double>();
            for(int i = N_variables - 1; i >= 0; --i){
              new_point->push_back( (**iter)[variables_id[i]] );
            }
            new_point->push_back( (**iter)[n_variables_original-2] );
            new_point->push_back( (**iter)[n_variables_original-1] );
            slice_chain->push_back( new_point );
            // remove this point
            chain_original_fl.erase_after( iter_before );
            if( std::next( iter_before )== chain_original_fl.end() ) break;
            sKIP_THIS_POINT:;
          }
          
          area_A = ordering_func->FindArea( slice ); // <- this is a A(\vec{p}')
          
          // check if \vec{x}_{exp} \in A(\vec{p}')
          if( area_A == nullptr or area_A->Contain(variables_exp_values) ){
            parameter_point = new std::vector<double>();
            parameter_point->reserve(N_parameters);
            for(i = 0; i < N_parameters; i++) parameter_point->push_back( (max_vals[i] + min_vals[i])*0.5 );
            set_I.push_back( parameter_point ); // add \vec{p}' to I(\vec{x}_{exp})
          }
          if(area_A != nullptr) areas_A.push_back(area_A);

          n_points++;
          // change parameters edges
          for(i = 0; i < N_parameters; i++){
            min_vals[i] = max_vals[i];
            max_vals[i] += par_range_steps[i];
            if( max_vals[i] <= par_range_maxs[i] ) break;
            min_vals[i] = par_range_mins[i];
            max_vals[i] = min_vals[i] + par_range_steps[i];
            if(i == N_parameters - 1) goto lOOP_EXIT;
            if(n_points > progress_point){
              msg_progress( double(progress_point)/set_I_size );
              progress_point = std::min(n_points + set_I_size / barWidth, set_I_size);
            }
            continue;
          }

          // debug printf
          /*
          std::cout << "------------" << std::endl;
          for(i = 0; i < N_parameters; i++){
            std::cout << parameters[i]->name << " " << min_vals[i] << " " << max_vals[i] << std::endl;
          }
          msg( area_A );
          msg( "slice", slice->chain.size() );
          if(area_A) msg( "contain = ", area_A->Contain(variables_exp_values) );
          */
        }
        lOOP_EXIT:
        msg_progress(1.);
        msg("\nshta::NeymanConstruction end Neyman construction ... ");
        return & set_I;
      }

      TCanvas * GetProjectionPlot(std::string label_x, std::string label_y, double step_x, double step_y){
        if(label_x == label_y){
          msg_err("WARNING : shta :: class NeymanConstruction : GetProjectionPlot() : label_x = label_y =", label_y, "return nullptr");
          return nullptr;
        }

        root::get_th2d_neyman_pallete();
        TCanvas * canv = root::get_canvas();

        const Variable * x = GetVariableOrParameter(label_x);
        const Variable * y = GetVariableOrParameter(label_y);
        const int par_x_id = GetParameterIndex( label_x );
        const int par_y_id = GetParameterIndex( label_y );
        const int var_x_id = GetVariableIndex( label_x );
        const int var_y_id = GetVariableIndex( label_y );

        // There is a three areas, first of all it is a Toys themselfs (Markov Chains)
        std::string hist_name = label_x + " & " + label_y;
        TH2D * neyman_plot = new TH2D(hist_name.c_str(), "", 
                               x->GetRangeLenght() / step_x, x->range_min, x->range_max,
                               y->GetRangeLenght() / step_y, y->range_min, y->range_max);
        neyman_plot->GetXaxis()->SetTitle( label_x.c_str() );
        neyman_plot->GetYaxis()->SetTitle( label_y.c_str() );

        mc_toys->Fill2DHist(label_x, label_y, neyman_plot);
        neyman_plot->Draw("colz");
        neyman_plot->SetStats(0);

        //Then it is a CL% areas A(\vec{p}') accepted by Neyman Construction
        for(auto interval : areas_A){
          
        }

        // Then it is a \vec{x}_{exp} line or point
        if(var_x_id != -1 or var_x_id != -1){
          double pos_x_min = x->range_min, pos_x_max = x->range_max, osmlx = (x->range_max - x->range_min)*0.001;
          double pos_y_min = y->range_min, pos_y_max = y->range_max, osmly = (y->range_max - y->range_min)*0.001;
          if(var_x_id != -1) { pos_x_min = x->value,  pos_x_max = x->value+osmly; }
          if(var_y_id != -1) { pos_y_min = y->value,  pos_y_max = y->value+osmly; }
          TLine *line = new TLine(pos_x_min, pos_y_min, pos_x_max, pos_y_max);
          line->SetLineColor(kOrange+9);
          line->SetLineWidth(4);
          line->SetLineStyle(2);
          line->Draw();
        }

        // And Finally it is a set I(\vec{x}_{exp}) in space \vec{p}
        // empty if projection on non-parameters plot
        int n_points = set_I.size();
        if(par_x_id == -1 and par_y_id == -1) return canv;
        else if(par_x_id != -1 and par_y_id != -1){
          TH2D * neyman_plot_area = new TH2D("set_I", "", 
                               x->GetRangeLenght() / step_x, x->range_min, x->range_max,
                               y->GetRangeLenght() / step_y, y->range_min, y->range_max);
          for(n_points--;n_points>=0;n_points--){
            std::vector<double> * point = set_I[n_points];
            neyman_plot_area->Fill((*point)[par_x_id], (*point)[par_y_id]);
          }
          neyman_plot_area->Draw("same box");
        }
        else{
          TGraph * points_on_line = new TGraph(n_points);
          double pos_x = x->value, pos_y = y->value;
          for(n_points--;n_points>=0;n_points--){
            std::vector<double> * point = set_I[n_points];
            if(par_x_id != -1) pos_x = (*point)[par_x_id];
            else               pos_y = (*point)[par_y_id];
            points_on_line->SetPoint(n_points, pos_x, pos_y);
          }
          points_on_line->Draw("same P");
          points_on_line->SetMarkerSize(1.);
          points_on_line->SetMarkerStyle(21);
          points_on_line->SetMarkerColor(kRed);
        }

        return canv;
      }

      void AddParameter(Variable * var){ parameters.push_back( var ); }
      void AddVariable(Variable * var) { variables.push_back( var );  }

      int GetParameterIndex(std::string name) const {
        for(int index = 0; index < parameters.size(); index++)
          if(parameters[index]->name == name) return index;
        return -1;
      }

      int GetVariableIndex(std::string name) const {
        for(int index = 0; index < variables.size(); index++)
          if(variables[index]->name == name) return index;
        return -1;
      }

      const Variable * GetVariableOrParameter(std::string name) const {
        for(const Variable * var : parameters)
          if(var->name == name) return var;
        for(const Variable * var : variables) 
          if(var->name == name) return var;
        msg_err("WARNING : shta :: class NeymanConstruction : GetVariableOrParameter() : cant find name =", name, "return nullptr");
        return nullptr;
      }
      
      std::vector<const Variable * > parameters;
      std::vector<const Variable * > variables;

      MarkovChain * mc_toys;

      OrderingFunction * ordering_func;

      std::vector<const Interval*> areas_A;
      std::vector< std::vector<double> * > set_I;
  };

  class NeymanConstruction1D{
    public:
      NeymanConstruction1D(Function * func, OrderingFunction * ord, std::string label_var, std::string label_par) : pdf_func(func), ordering_func(ord) {
        pdf_eval = pdf_func->eval;
        neyman_plot = nullptr;
        intervals_area = nullptr;
        parameter_area = nullptr;
        variable_label = label_var;
        parameter_label = label_par;
      }

      ~NeymanConstruction1D(){
        delete neyman_plot;
        delete intervals_area;
        delete parameter_area;
      }

      TCanvas * GetPlot(){
        if( is_eq_one(nullptr, neyman_plot, intervals_area, parameter_area) ){
          msg_err("WARNING : shta :: class NeymanConstruction1D : GetPlot() : neyman_plot, intervals_area, parameter_area = ", neyman_plot, intervals_area, parameter_area, " return NULL");
          return nullptr;
        }

        root::get_th2d_neyman_pallete();
        TCanvas * canv = root::get_canvas();
        neyman_plot->Draw("colz");
        neyman_plot->SetStats(0);
        msg("5");
        /*
        intervals_area->sh->Draw("same f");
        intervals_area->sh->SetFillColor(30);
        intervals_area->sh->SetFillStyle(3350);
        */

        const Variable * var = pdf_func->GetVariable( variable_label );
        const Variable * par = pdf_func->GetVariable( parameter_label );
        if( is_eq_one(nullptr, var, par) ){
          msg_err("WARNING : shta :: class NeymanConstruction1D : GetPlot() : cant find variables with labels:", variable_label, parameter_label, " in this function");
          return nullptr;
        }

        TLine *line = new TLine(var->value, par->range_min, var->value, par->range_max);
        line->SetLineColor(kOrange+9);
        line->SetLineWidth(4);
        line->SetLineStyle(2);
        line->Draw();

        intervals_area->l->Draw("same L");
        intervals_area->l->SetLineColor(1);
        intervals_area->l->SetLineWidth(3);
        intervals_area->l->SetLineStyle(9);

        intervals_area->r->Draw("same L");
        intervals_area->r->SetLineColor(1);
        intervals_area->r->SetLineWidth(3);
        intervals_area->r->SetLineStyle(9);

        parameter_area->Draw("same P");
        parameter_area->SetMarkerColor(kOrange-3);
        parameter_area->SetMarkerStyle(21);
        parameter_area->SetMarkerSize(1.0);

        long long int size = parameter_area->GetN();
        if(not size) return canv;
        double max_y = *std::max_element(parameter_area->GetY(), parameter_area->GetY() + size);
        double min_y = *std::min_element(parameter_area->GetY(), parameter_area->GetY() + size);

        TArrow *ar1 = new TArrow(var->range_min, max_y, var->value, max_y, 0.02, "<|");
        ar1->SetLineColor(kOrange+9);
        ar1->SetFillColor(kOrange+9);
        ar1->SetLineWidth(3);
        ar1->Draw();

        TArrow *ar2 = new TArrow(var->range_min, min_y, var->value, min_y, 0.02, "<|");
        ar2->SetLineColor(kOrange+9);
        ar2->SetFillColor(kOrange+9);
        ar2->SetLineWidth(3);
        ar2->Draw();

        TLegend * legend = new TLegend(0.45,0.15,0.88,0.4);
        //legend->AddEntry(neyman_plot,       ("P("+var->GetName()+"|"+par->GetName()+")").c_str());
        legend->AddEntry(line,              (var->GetName() + " = " + std::to_string(var->value)).c_str(), "l");
        legend->AddEntry(intervals_area->l, "Neyman construction area",  "l");
        legend->AddEntry(parameter_area,    ("Neyman confedence interval for " + par->GetName()).c_str(),  "p");
        legend->Draw();

        return canv;
      }

      std::vector<double> * ConstructIntervals(double step_var, double step_par){
        if(not pdf_func or not ordering_func){
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervals() : do not set pdf_func, ordering_func = ", pdf_func, ordering_func);
          return nullptr;
        }
        if(pdf_func->GetVariablesNumber() + pdf_func->GetParametersNumber() < 2){
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervalsD() : wrong PDFunction dimension ", pdf_func->GetVariablesNumber(), " + ", pdf_func->GetParametersNumber(), "return;");
          return nullptr;
        }

        const Variable * var = pdf_func->GetVariable( variable_label );
        const Variable * par = pdf_func->GetVariable( parameter_label );
        if( is_eq_one(nullptr, var, par) ){
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervalsD() : cant find variables with labels:", variable_label, parameter_label, " in this function");
          return nullptr;
        }

        if(not var->CheckDefinition() or not par->CheckDefinition()){
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervalsD() : definition of variables is uncorrect: ", variable_label, parameter_label, ", return");
          return nullptr;
        }

        int var_bins = (var->range_max - var->range_min)/step_var;
        int par_bins = (par->range_max - par->range_min)/step_par;
        if( var_bins < 1 or par_bins < 1 ){
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervalsD() : step of integration is greater than parameter/variable range:", variable_label, parameter_label, ", return");
          return nullptr;
        }
        if( var_bins > 50000 or par_bins > 50000){
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervalsD() : number of bins for intergation is too big:", var_bins, par_bins, ", try to continue");
        }

        int index_var = pdf_func->GetVariableId( variable_label );
        int index_par = pdf_func->GetVariableId( parameter_label );

        int type_index_var = var->IsVariable() ? 0 : 1; // for now its only to possibility
        int type_index_par = par->IsVariable() ? 0 : 1;

        msg("INFO : shta :: class NeymanConstruction1D : ConstructIntervals() : begin construction ... ");
        std::vector <double> * par_set = new std::vector <double>();

        // for every possible value of <par> construct a 1-D pdf for <var>
        // we will use integration for this as this is a one dimension only and store result in TH2D
        // and in vector
        delete neyman_plot;
        delete intervals_area;
        std::string name = root::get_name();
        // msg("debug 0.5", var_bins, par_bins, var->value, var->range_min, var->range_max, par->range_min, par->range_max);
        neyman_plot = new TH2D(name.c_str(), "", 
                         var_bins, var->range_min, var->range_max, 
                         par_bins, par->range_min, par->range_max);
        neyman_plot->GetXaxis()->SetTitle( var->GetName().c_str() );
        neyman_plot->GetYaxis()->SetTitle( par->GetName().c_str() );

        intervals_area = new root::TGraphShade( par_bins - 1 );
        parameter_area = new TGraph( par_bins );

        // fill variable/parameter vector
        std::vector<double> * vars = pdf_func->GetVariablesVector();
        std::vector<double> * pars = pdf_func->GetParametersVector();
        std::vector<std::vector<double>* > pair_vector = {vars, pars};
        double pdf_val;
        std::vector <double> slice;
        slice.reserve( (var->range_max - var->range_min)/step_var );
        long long int index = 0;
        for(double par_ = par->range_min + 0.5 * step_par; par_ < par->range_max; par_ += step_par){
          slice.clear();
          for(double var_ = var->range_min + 0.5 * step_var; var_ < var->range_max; var_ += step_var){
            (*(pair_vector.at(type_index_var)))[index_var] = var_;
            (*(pair_vector.at(type_index_par)))[index_par] = par_;
            pdf_val = pdf_eval( *vars, *pars );
            neyman_plot->Fill( var_, par_, pdf_val );
            slice.push_back( pdf_val );
          }
          // now for every slice we would like to get an interval
          Interval1D * interval = ordering_func->GetInterval( slice );
          interval->lower_abs = var->range_min + 0.5 * step_var + interval->lower_index * step_var;
          interval->upper_abs = var->range_min + 0.5 * step_var + interval->upper_index * step_var;

          // check, if this interval contain a given experimental value of <var>
          if( interval->Contain( var->value ) ) par_set->push_back( par_ );

          // in addition we will store intervals edges and our finded <par>s in TF1s
          intervals_area->SetPoint( index, interval->lower_abs, par_, interval->upper_abs, par_ );
          
          index++;
        }

        msg("INFO : shta :: class NeymanConstruction1D : ConstructIntervals() : end construction ... ");
        if(not par_set->size())
          msg_err("WARNING : shta :: class NeymanConstruction1D : ConstructIntervalsD() : Neyman interval is empty ... ");

        parameter_area = new TGraph( par_set->size() );
        for(unsigned long long int i = 0; i < par_set->size(); i++) parameter_area->SetPoint(i, var->value, (*par_set)[i]);
        
        return par_set;
      }

      const Function * pdf_func;
      OrderingFunction * ordering_func;
      double (*pdf_eval)(std::vector<double> & x, std::vector<double> & p);
      std::string variable_label, parameter_label;

      TH2D * neyman_plot;
      TGraph * parameter_area;
      root::TGraphShade * intervals_area;
  };
};

#endif 
