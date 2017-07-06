
//  This is  SHTA package for statistical data analyses in high energy physics
//  SHTA available under GNU Lesser General Public License v3.0
//  more information in README file
// 
//  Mandrik P., IHEP, PROTVINO, 2017 

#ifndef SHTA_MSG_HH
#define SHTA_MSG_HH 1

namespace shta {

  template<typename T> void print_simple_vector(std::vector <T> v, std::string separator = "\n"){
    int size = v.size();
    for(int i = 0; i < size; ++i) std::cout << v[i] << separator;
    std::cout << std::endl;
  }

  // msg
  void msg(){ std::cout << "\n" << std::endl; };
  template<typename T> void msg(T t){
    std::cout << t << std::endl;
  }

  template<typename T, typename... Args> void msg(T t, Args... args){
    std::cout << t << " ";
    msg(args...);
  }

  // msg_err
  void msg_err(){ std::cerr << "\n" << std::endl; };
  template<typename T> void msg_err(T t){
    std::cerr << t << std::endl;
  }

  template<typename T, typename... Args> void msg_err(T t, Args... args){
    std::cerr << t << " ";
    msg_err(args...);
  }

  enum verbose{
    SILENCE=0,
    ERROR,
    WARNING,
    INFO,
    DEBUG,
    VERBOSE
  };

  #define MSG_ERROR(...) if(verbose_lvl >= shta::verbose::ERROR)   msg_err(__VA_ARGS__)
  #define MSG_WARNING(...) if(verbose_lvl >= shta::verbose::WARNING) msg_err(__VA_ARGS__)
  #define MSG_INFO(...) if(verbose_lvl >= shta::verbose::INFO)    msg(__VA_ARGS__)
  #define MSG_DEBUG(...) if(verbose_lvl >= shta::verbose::DEBUG)   msg(__VA_ARGS__)
  #define MSG_VERBOSE(...) if(verbose_lvl >= shta::verbose::VERBOSE) msg(__VA_ARGS__)

  void msg_progress(double progress, int barWidth = 50){
    if(progress > 1.) progress = 1.;
    if(progress < 0.) progress = 0.;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }

  void write_string_to_file(std::string const & file_name, std::string const & data){
  std::ofstream file( file_name );
  if(!file){
    msg_err("WARNING : shta : write_string_to_file() : can’t open output file \"", file_name, "\"");
  }

  file << data;
}
};

#endif
