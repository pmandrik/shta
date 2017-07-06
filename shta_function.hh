
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_FUNCTION_HH
#define SHTA_FUNCTION_HH 1

namespace shta {

  // SetName, GetName, Print, operators +, -, *, /

  static const std::string _functionSeparators = " `~!@#$%^&*()+|-={}\'\";:<>?/,";

  enum vartype{
    VARIABLE = 1,
    PARAMETER
  };

  class Variable {
    public:
      Variable(std::string n, double val, double rmin, double rmax, int typo=shta::vartype::VARIABLE) : name(n), value(val), range_min(rmin), range_max(rmax), type(typo) {};

      std::string GetName() const { return name; };

      void Print() const {
        msg( name, "=", value, "[", range_min, range_max, " ]" );
      }

      bool CheckDefinition() const {
        if( range_min > range_max ) return false;
        if( value < range_min or value > range_max ) return false;
        return true;
      }

      bool IsParameter() const { return (type == shta::vartype::PARAMETER ? true : false); }
      bool IsVariable() const {  return (type == shta::vartype::VARIABLE  ? true : false); }
      double GetRangeLenght() const { return range_max - range_min; }

      TH1D * GetTH1D(int bins) const {
        return new TH1D(root::get_name().c_str(), name.c_str(), bins, range_min, range_max);;
      }

      std::string name;
      double value, range_min, range_max;
      int type;

      std::map<std::string, double> extra_variables;
  };

  class Function {
    public:
      Function(std::string name) : name(name) { }
      Function(std::string name, std::string def) : name(name), definition(def) {
        initial_definition = def;
        eval = NULL;
      }

      void SetDefinition(std::string def){ definition = def; initial_definition = def; }
      void SetName(std::string def){ definition = def; }

      std::string GetName() const { return name; }
      std::string GetDefinition() const { return definition; }
      int GetDepSize() const { return dependences.size(); }
      int GetVariablesNumber()  const { return used_variables.size(); }
      int GetParametersNumber() const { return used_parameters.size(); }

      void Print() const {
        msg( name, definition );
      }

      void ParceDefinition(const std::vector <Function*> & functions){
        substrings.clear();
        unsigned int sub_i = 0;
        for ( unsigned int i = 0 ; i < definition.length(); i++){
          if( _functionSeparators.find( definition[i] ) == std::string::npos ) continue;
          if(i-sub_i) substrings.push_back( definition.substr(sub_i, i-sub_i) );
          substrings.push_back( definition.substr(i, 1) );
          sub_i = i+1;
        }
        substrings.push_back( definition.substr(sub_i, definition.length()) );

        dependences.clear();
        for(auto item : substrings){
          auto it = find_if( functions.begin(), functions.end(), [&item](const Function * f) {return f->GetName() == item;} );
          if(it == functions.end()) continue;
          dependences.push_back( *it );
        }
      }

      void ResetDefinition(){
        definition = "";
        for(auto item : substrings) definition += item;
      };

      bool ExpandFunction(const Function* func){
        for(auto it = dependences.begin(); it != dependences.end(); ){
          if( func->GetName() == (*it)->GetName() ) it = dependences.erase( it );
          else ++it;
        }

        bool answer = false;
        for(auto it = substrings.begin(); it != substrings.end(); ){
          if(func->GetName() == *it) {
            *it = "";
            substrings.insert( it, func->substrings.begin(), func->substrings.end() );
            answer = true;
            it = substrings.begin();
            continue;
          }
          ++it; 
        }
        return answer;
      }

      void ExpandVariables(const std::vector <Variable*> & variables){
        variables_to_array_map.clear();
        std::string prefix;
        int index;
        for(auto it = substrings.begin(); it != substrings.end(); it++){
          auto finded = find_if( variables.begin(), variables.end(), [it](const Variable * v) {return v->GetName() == *it;} );
          if(finded == variables.end()) continue;
          if( (*finded)->IsVariable() ) prefix = "x";
          else                          prefix = "p";

          auto find_id = variables_to_array_map.find( (*finded)->GetName() );
          if( find_id == variables_to_array_map.end() ){
            if( (*finded)->IsVariable() ){ 
              index = used_variables.size();
              used_variables.push_back( (*finded) ); 
            }
            else{
              index = used_parameters.size();
              used_parameters.push_back( (*finded) ); 
            }
            variables_to_array_map[ (*finded)->GetName() ] = std::make_pair(index, (*finded));
          }
          else index = find_id->second.first;
          *it = prefix + "[" + std::to_string( index ) + "]";
        }
      }

      std::string GetSource () const {
        std::string src = "// " + name + ", \'" + initial_definition + "\'\n";
        src += "extern \"C\" double " + name + "(std::vector<double> & x, std::vector<double> & p){\n";
        src += "  return ( /* -> */ " + definition + " /* <- */ );\n";
        src += "}\n";
        return src;
      }

      int GetVariableId (std::string name) const {
        auto find_id = variables_to_array_map.find( name );
        if( find_id == variables_to_array_map.end() ) return -1;
        return find_id->second.first;
      }

      Variable * GetVariable (std::string name) const {
        auto find_id = variables_to_array_map.find( name );
        if( find_id == variables_to_array_map.end() ) return nullptr;
        return find_id->second.second;
      }

      void MapArray(double * x, double * p) const {
      }

      std::vector<double> * GetVariablesVector() const {
        std::vector<double> * answer = new std::vector<double>();
        answer->reserve( used_variables.size() );
        for(auto var : used_variables) answer->push_back( var->value );
        return answer;
      }

      std::vector<double> * GetParametersVector() const {
        std::vector<double> * answer = new std::vector<double>();
        answer->reserve( used_parameters.size() );
        for(auto var : used_parameters) answer->push_back( var->value );
        return answer;
      }

      double EvalTest(){
        std::vector <double> vars;
        std::vector <double> pars;
        for(const Variable * var : used_variables) vars.push_back( var->value );
        for(const Variable * par : used_variables) pars.push_back( par->value );
        return (this->eval(vars, pars));
      }

      std::string name, definition, initial_definition;
      std::vector <std::string> substrings;
      std::vector <Function*> dependences;
      std::vector <Variable*> used_variables;
      std::vector <Variable*> used_parameters;
      std::map <std::string, std::pair<int, Variable*> > variables_to_array_map;
      double (*eval)(std::vector<double> & x, std::vector<double> & p);
  };

  class FunctionBuilder {
    public:
      FunctionBuilder(){
        preambule += "/*  ";
        preambule += "This src file was generated automatically by shta::FunctionBuilder";
        preambule += "  */ \n ";
        name_space = "test";
        name_out   = "test.hh";

        AddInclude("shta_core.hh");

        func_library_handle = nullptr;
      }

      ~FunctionBuilder(){
        for( Variable * var  : variables) delete var;
        for( Function * func : functions) delete func;
        if(func_library_handle != nullptr) dlclose(func_library_handle);
      }

      Variable * DefineVariable(std::string name, double val, double rmin, double rmax, int type=shta::vartype::VARIABLE){
        Variable * var = new Variable(name, val, rmin, rmax, type);
        variables.push_back( var );
        return var;
      }

      Function * DefineFunction(std::string name, std::string def){
        Function * func = new Function(name, def);
        functions.push_back( func );
        return func;
      }

      void AddInclude(std::string name){
        includes.push_back( name );
      }

      void Print(){
        msg("Variables : ");
        for(auto var : variables){
          var->Print();
        }
        msg("Functions : ");
        for(auto func : functions){
          func->Print();
        }
      }

      bool CollectDefinitions(){
        /*
          - iterate over function and check, if they use other functions
            construct a graph of functions dependence and check 
            if it is possible to expand them
          - expand functions, replace variables to array items and create dictionary name:order
        */

        // Parse
        for(Function * func : functions){
          func->ParceDefinition( functions );
        }

        // resolve function dependences (hard-coding, n^3)
        std::vector <Function*> functions_expanded;
        for(int resolves = 1; resolves; ){
          resolves=0;
          for(Function * func_to_expand : functions)
            if( func_to_expand->GetDepSize() == 0 )
              for(Function * func : functions)
                resolves += func->ExpandFunction( func_to_expand );
        }

        // resolve variables
        for(Function * func : functions){
          func->ResetDefinition();
          func->ExpandVariables( variables );
        }

        // check and exit
        for(Function * func : functions){
          if( func->GetDepSize() != 0 ){
            msg_err("WARNING : shta :: class FunctionBuilder : CollectDefinitions() : function ", func->GetName(), "has unresolved function dependences");
          }
          func->ResetDefinition();
        }

        return true;
      }

      std::string GetSource(){
        std::string src = preambule + "\n";
        for(auto include : includes)   src += "#include <" + include + ">\n";
        src += "namespace " + name_space + " { ";
        for(auto function : functions) src += function->GetSource() + "\n";
        src += "}; ";
        return src;
      }

      void CreateSrcFile(std::string name){
        write_string_to_file( get_shta_userfunctions_path() + "/" + name + ".cpp", GetSource());
      }

      void CompileFunctions(std::string name){
        std::string fname = "shta_uftd_" + name;
        system( ("mkdir "+fname+"; cd "+fname+"; cmake -DUSER_FUNCTION_NAME=" + get_shta_userfunctions_path() + "/" + name + ".cpp ../.. ; make; cp lib" + name + ".so ../.").c_str());
      }

      void LoadFunctionsLibrary(std::string name){
        func_library_handle = dlopen((get_shta_exe_path() + "/" + name).c_str(), RTLD_LAZY);
        if(!func_library_handle) msg_err("WARNING : shta :: class FunctionBuilder : LoadFunctionsLibrary() : Cannot open library: ", dlerror());
      }

      void LinkFunction(Function * func){
        if(!func_library_handle) msg_err("WARNING : shta :: class FunctionBuilder : LinkFunction() : function library is not loaded");
        dlerror();
        func->eval = (double (*)(std::vector<double> & x, std::vector<double> & p)) dlsym(func_library_handle, func->name.c_str());
        char *error;
        if ((error = dlerror()) != NULL){ 
          msg_err("WARNING : shta :: class FunctionBuilder : LinkFunction() : dlsym error for function ", func->name, " : ");
          fputs(error, stderr);
          func->eval = NULL;
        };
      }

      void LinkFunctions(){
        for(Function* func : functions) LinkFunction(func);
      }

      void BuildFunctions(std::string pname){
        CollectDefinitions();
        CreateSrcFile( pname );
        CompileFunctions( pname );
        LoadFunctionsLibrary("lib" + pname + ".so");
        LinkFunctions();
      }

      std::vector <Function*> functions;
      std::vector <Variable*> variables;

      // source code variables
      std::vector <std::string> includes;
      std::string name_out;
      std::string preambule;
      std::string name_space;
      void * func_library_handle;
  };

};

#endif


