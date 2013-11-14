#include "trainingmode.hpp"

namespace Training {
  std::istream& operator>>(std::istream& in, Method& method)
  {
    std::string token;
    in >> token;
    if(token == "none" or token == "fixed" or token == "fix")
      method = Method::None;
    else if(token == "reestimation" or token == "reestimate" or token == "em" or token == "generative")
      method = Method::Reestimation;
    else if(token == "gradient" or token == "discriminative" or token == "disc")
      method = Method::Gradient;
    else {
      std::cout << "Training method '" << token << "' not implemented. See -h or --help for help." << std::endl;
      exit(-1);
    }
    return(in);
  }

  Method measure2method(Measure measure)
  {
    if(Measures::is_generative(measure))
      return(Method::Reestimation);
    else if(Measures::is_discriminative(measure))
      return(Method::Gradient);
    else
      return(Method::None);
  }


  std::string method2string(Method method)
  {
    std::string s;
    switch(method) {
      case Method::None:
        s = "none";
        break;
      case Method::Reestimation:
        s = "reestimation";
        break;
      case Method::Gradient:
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


  Seeding::Objective corresponding_objective(const Objective &x, bool use_mi_to_seed) {
    Seeding::Objective y;
    y.motif_name = x.motif_name;
    y.series_expression = x.series_expression;
    if(use_mi_to_seed)
      y.measure = Measures::Discrete::Measure::MutualInformation;
    else
      y.measure = Measures::corresponding_measure(x.measure);
    return(y);
  }
  Seeding::Objectives corresponding_objectives(const Objectives &x, bool use_mi_to_seed) {
    Seeding::Objectives y;
    for(auto &z: x)
      y.push_back(corresponding_objective(z, use_mi_to_seed));
    return(y);
  }
}
