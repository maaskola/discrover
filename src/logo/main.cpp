#include <iostream>
#include <cstdlib>
#include "logo.hpp"

using namespace std;

int main(int argc, char**argv) {
  if(argc < 2) {
    cout << "Please provide the paths to one or more .hmm files." << endl;
    return -1;
  }

  for(int i = 1; i < argc; ++i)
    if(draw_logo(argv[i])) {
      cout << "Problem drawing logo for " << argv[i] << endl;
      return -1;
    }

  return EXIT_SUCCESS;
}
