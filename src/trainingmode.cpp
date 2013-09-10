#include "trainingmode.hpp"

namespace Training {
  std::istream& operator>>(std::istream& in, Method& method)
  {
    std::string token;
    in >> token;
    if(token == "none" or token == "fixed" or token == "fix")
      method = Method::none;
    else if(token == "reestimation" or token == "reestimate" or token == "em" or token == "generative")
      method = Method::reestimation;
    else if(token == "gradient" or token == "discriminative" or token == "disc")
      method = Method::gradient;
    else {
      std::cout << "Training method '" << token << "' not implemented. See -h or --help for help." << std::endl;
      exit(-1);
    }
    return(in);
  }

  Method measure2method(Measure measure)
  {
    if(Measures::is_generative(measure))
      return(Method::reestimation);
    else if(Measures::is_discriminative(measure))
      return(Method::gradient);
    else
      return(Method::none);
  }


  std::string method2string(Method method)
  {
    std::string s;
    switch(method) {
      case Method::none:
        s = "none";
        break;
      case Method::reestimation:
        s = "reestimation";
        break;
      case Method::gradient:
        s = "gradient";
        break;
    }
    return(s);
  }

  std::ostream& operator<<(std::ostream& out, Method method)
  {
    out << method2string(method);
    return(out);
  }


  Plasma::Objective corresponding_objective(const Objective &x, bool use_mi_to_seed) {
    Plasma::Objective y;
    y.motif_name = x.motif_name;
    y.series_expression = x.series_expression;
    if(use_mi_to_seed)
      y.measure = Measures::Discrete::Measure::MutualInformation;
    else
      y.measure = Measures::corresponding_measure(x.measure);
    return(y);
  }
  Plasma::Objectives corresponding_objectives(const Objectives &x, bool use_mi_to_seed) {
    Plasma::Objectives y;
    for(auto &z: x)
      y.push_back(corresponding_objective(z, use_mi_to_seed));
    return(y);
  }
}
