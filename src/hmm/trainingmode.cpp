#include "trainingmode.hpp"

using namespace std;

namespace Training {
istream& operator>>(istream& in, Method& method) {
  string token;
  in >> token;
  if (token == "none" or token == "fixed" or token == "fix")
    method = Method::None;
  else if (token == "reestimation" or token == "reestimate" or token == "em"
           or token == "generative")
    method = Method::Reestimation;
  else if (token == "gradient" or token == "discriminative" or token == "disc")
    method = Method::Gradient;
  else {
    cout << "Training method '" << token
         << "' not implemented. See -h or --help for help." << endl;
    exit(-1);
  }
  return in;
}

Method measure2method(Measure measure) {
  if (Measures::is_generative(measure))
    return Method::Reestimation;
  else if (Measures::is_discriminative(measure))
    return Method::Gradient;
  else
    return Method::None;
}

string method2string(Method method) {
  string s;
  switch (method) {
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
  return s;
}

ostream& operator<<(ostream& out, Method method) {
  out << method2string(method);
  return out;
}

Seeding::Objective corresponding_objective(const Objective& x,
                                           bool use_mi_to_seed) {
  Seeding::Objective y;
  y.motif_name = x.motif_name;
  y.contrast_expression = x.contrast_expression;
  if (use_mi_to_seed)
    y.measure = Measures::Discrete::Measure::MutualInformation;
  else
    y.measure = Measures::corresponding_measure(x.measure);
  return y;
}
Seeding::Objectives corresponding_objectives(const Objectives& x,
                                             bool use_mi_to_seed) {
  Seeding::Objectives y;
  for (auto& z : x)
    y.push_back(corresponding_objective(z, use_mi_to_seed));
  return y;
}
}
