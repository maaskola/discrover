#include <iostream>
#include <cstdlib>
#include "logo.hpp"

using namespace std;

int main(int argc, char**argv) {
  if(argc < 2) {
    cout << "Please provide the paths to one or more .hmm files." << endl;
    return -1;
  }

  Logo::matrix_t matrix = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
    {0.5, 0.5, 0, 0},
    {0.5, 0, 0.5, 0},
    {0.5, 0, 0, 0.5},
    {1.0/3, 1.0/3, 1.0/3, 0},
    {1.0/3, 1.0/3, 0, 1.0/3},
    {1.0/3, 0, 1.0/3, 1.0/3},
    {0, 1.0/3, 1.0/3, 1.0/3}
  };

  Logo::Options options;

  for(int i = 1; i < argc; ++i) {
    auto paths = draw_logo(matrix, argv[i], options);
    if(paths.size() == 0) {
      cout << "Problem drawing logos for " << argv[i] << endl;
      return -1;
    }
  }

  return EXIT_SUCCESS;
}
