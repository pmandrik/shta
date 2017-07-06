
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_CLS_HH
#define SHTA_CLS_HH 1 

namespace shta {

/*
  Classes for Modified Frequentist Approach for upper limit settings
*/

  class ModifiedFrequentistApproach1D {
    public:

      ModifiedFrequentistApproach1D(Function * func_sb, Function * func_b, std::string label_var, std::string label_par) : 
                                    pdf_sb_func(func_sb), pdf_b_func(func_b), variable_label(label_var), parameter_label(label_par) {
        pdf_sb_eval = pdf_sb_func->eval;
        pdf_b_eval = pdf_b_func->eval;

        cl_vs_p_graph = nullptr;
        cl_vs_p_graph_accepted = nullptr;
      }

      /*
        Modified Frequentist Approach or CLs methos is a tool for setting upper limits
        on parameters of interests. In this class we introduce a single parameter \vec{p} dimension realisation
        of CLs without nuisances and with only one dimension of ramdom variables \vec{x}.
        
        CLs+b = \int_{ Q(x) < Q(x_{obs}) } P_{s+b}(x|p_{up}) dx 
        CLb   = \int_{ Q(x) < Q(x_{obs}) } P_{b}  (x|p_{up}) dx 
        CLs   = CLs+b / CLb and usuallu is a 0.95 or greater =>
        p_{up} = f(CLs, x_{obs})
        where Q is a test statistics, originaly in CLs is Likelihood ratio => function of parameter too.

        Strainforvard realisation is a calcullate <CLs> as a function of <p> and then find a point where disired <CLs> intersect this function
      */

      std::vector<double> * FindUpperLimit(double CL, double step_var, double step_par){
        // preparation, similare to neyman, but double work because of two pdfs
        const Variable * var_sb = pdf_sb_func->GetVariable( variable_label );
        const Variable * par_sb = pdf_sb_func->GetVariable( parameter_label );
        const Variable * var_b = pdf_b_func->GetVariable( variable_label );
        const Variable * par_b = pdf_b_func->GetVariable( parameter_label );
        if( is_eq_one(nullptr, var_sb, par_sb, var_b) ){
          msg_err("WARNING : shta :: class ModifiedFrequentistApproach1D : FindUpperLimit() : cant find variables with labels:", variable_label, parameter_label, " in this function");
            return nullptr;
        }

        int index_var_sb = pdf_sb_func->GetVariableId( variable_label );
        int index_par_sb = pdf_sb_func->GetVariableId( parameter_label );
        int index_var_b = pdf_b_func->GetVariableId( variable_label );

        int type_index_var_sb = var_sb->IsVariable() ? 0 : 1; // for now its only to possibility
        int type_index_par_sb = par_sb->IsVariable() ? 0 : 1;
        int type_index_var_b = var_b->IsVariable() ? 0 : 1;

        int index_par_b, type_index_par_b;
        if(par_b != nullptr){
          index_par_b = pdf_b_func->GetVariableId( parameter_label );
          type_index_par_b = par_b->IsVariable() ? 0 : 1;
        }

        msg("INFO : shta :: class ModifiedFrequentistApproach1D : FindUpperLimit() : start ... ");

        // for every possible value of <par> find a <CL>
        std::vector<std::pair<double, double> > cl_vs_p;
        double Q_obs, Q, CLs;
        double CLsb, Psb_integral, Psb; 
        double CLb,  Pb_integral,  Pb;
        double var_obs = var_sb->value;
        CLs_max = 1.-CL;

        std::vector<double> * vars_sb = pdf_sb_func->GetVariablesVector();
        std::vector<double> * pars_sb = pdf_sb_func->GetParametersVector();
        std::vector<std::vector<double>* > pair_vector_sb = {vars_sb, pars_sb};

        std::vector<double> * vars_b = pdf_b_func->GetVariablesVector();
        std::vector<double> * pars_b = pdf_b_func->GetParametersVector();
        std::vector<std::vector<double>* > pair_vector_b = {vars_b, pars_b};

        std::vector <double> * par_set = new std::vector <double>();

        for(double par_ = par_sb->range_min + 0.5 * step_par; par_ < par_sb->range_max; par_ += step_par){
          CLb = 0.;
          CLsb = 0.;
          Pb_integral  = 0.;
          Psb_integral = 0.;

          (*(pair_vector_b.at(type_index_var_b)))[index_var_b] = var_obs;
          if(par_b != nullptr) (*(pair_vector_b.at(type_index_par_b)))[index_par_b] = par_;
          (*(pair_vector_sb.at(type_index_var_sb)))[index_var_sb] = var_obs;
          (*(pair_vector_sb.at(type_index_par_sb)))[index_par_sb] = par_;
          Pb = pdf_b_eval(  *vars_b, *pars_b );
          if(Pb <= 0) continue; // FIXME
          Q_obs = pdf_sb_eval( *vars_sb, *pars_sb ) / Pb; // FIXME not normalized :c

          for(double var_ = var_sb->range_min + 0.5 * step_var; var_ < var_sb->range_max; var_ += step_var){
            (*(pair_vector_b.at(type_index_var_b)))[index_var_b] = var_;
            (*(pair_vector_sb.at(type_index_var_sb)))[index_var_sb] = var_;
            Pb = pdf_b_eval(  *vars_b, *pars_b );
            Psb = pdf_sb_eval( *vars_sb, *pars_sb );
            if(Pb <= 0) continue; // FIXME
            Q = Psb / Pb;
            Psb_integral += Psb;
            Pb_integral  += Pb;
            if(Q > Q_obs) continue;
            CLsb += Psb;
            CLb  += Pb;
          }
          if(CLb <= 0 or Psb_integral <= 0 or Pb_integral <= 0) continue;
          CLs = CLsb / CLb * Pb_integral / Psb_integral; // CLs is >= than frequency of wrong exclusion
          cl_vs_p.push_back( std::make_pair(CLs, par_) );
          if(CLs < CLs_max) par_set->push_back(par_);
          //msg( par_, CLs, CLsb, CLb, Pb_integral, Psb_integral, 1. - CLs, CL );
        }

        // fill a graph
        delete cl_vs_p_graph;
        delete cl_vs_p_graph_accepted;
        cl_vs_p_graph          = new TGraph( cl_vs_p.size() );
        cl_vs_p_graph_accepted = new TGraph( cl_vs_p.size() - par_set->size() );
        cl_vs_p_graph->GetXaxis()->SetTitle( variable_label.c_str() );
        cl_vs_p_graph->GetYaxis()->SetTitle( "CLs" );
        int index = 0;
        for(int i = 0; i < cl_vs_p.size(); i++){
          cl_vs_p_graph->SetPoint(i, (cl_vs_p.at(i)).second, (cl_vs_p.at(i)).first );
          if((cl_vs_p.at(i)).first > CLs_max){
            cl_vs_p_graph_accepted->SetPoint(index, (cl_vs_p.at(i)).second, (cl_vs_p.at(i)).first );
            index++;
          }
        }

        return par_set;
      }

      TCanvas * GetPlot(){
        TCanvas * canv = root::get_canvas();
        canv->SetLogy();

        cl_vs_p_graph->SetLineColor(1);
        cl_vs_p_graph->SetMarkerColor(kRed);
        cl_vs_p_graph->SetMarkerStyle(21);
        cl_vs_p_graph->SetMarkerSize(0.5);
        cl_vs_p_graph->Draw("APL");

        cl_vs_p_graph_accepted->SetLineColor(1);
        cl_vs_p_graph_accepted->SetMarkerColor(kGreen);
        cl_vs_p_graph_accepted->SetMarkerStyle(21);
        cl_vs_p_graph_accepted->SetMarkerSize(0.5);
        cl_vs_p_graph_accepted->Draw("PL same");

        const Variable * var = pdf_sb_func->GetVariable( variable_label );
        const Variable * par = pdf_sb_func->GetVariable( parameter_label );

        TLine *line = new TLine(var->value, 0, var->value, 1.);
        line->SetLineColor(kBlack);
        line->SetLineWidth(4);
        line->SetLineStyle(2);
        line->Draw();

        TLine *line_cls = new TLine(par->range_min, CLs_max, par->range_max, CLs_max);
        line_cls->SetLineColor(kGreen-5);
        line_cls->SetLineWidth(4);
        line_cls->SetLineStyle(2);
        line_cls->Draw();

        TLegend * legend = new TLegend(0.45,0.15,0.88,0.4);
        legend->AddEntry(line,              (var->GetName() + " = " + std::to_string(var->value)).c_str(), "l");
        legend->AddEntry(line_cls,          ("CLs = " + std::to_string(CLs_max)).c_str(),  "l");
        legend->AddEntry(cl_vs_p_graph_accepted,   ("Accepted " + par->GetName()).c_str(),  "pl");
        legend->AddEntry(cl_vs_p_graph,            ("Excluded " + par->GetName()).c_str(),  "pl");
        legend->Draw();

        return canv;
      }

      const Function * pdf_sb_func, * pdf_b_func;
      double (*pdf_sb_eval)(std::vector<double> & x, std::vector<double> & p);
      double (*pdf_b_eval) (std::vector<double> & x, std::vector<double> & p);

      double CLs_max;
      std::string variable_label, parameter_label;

      TGraph * cl_vs_p_graph;
      TGraph * cl_vs_p_graph_accepted;
  };
}

#endif 
