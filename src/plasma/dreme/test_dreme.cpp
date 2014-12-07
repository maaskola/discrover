#include "dreme.hpp"
#include <iostream>

using namespace std;

int main(int argc, const char** argv) {
  if (argc < 2) {
    cout << "Please provide at least one path to a FASTA file." << endl;
    return -1;
  }
  if (argc > 3) {
    cout << "Please provide at most two paths to a FASTA file." << endl;
    return -2;
  }

  string path1 = argv[1];
  string path2 = "";

  if (argc == 3)
    path2 = argv[2];

  try {
    auto regexes = Dreme::run(path1, path2, 8, 8, false, 1);
    for (auto& x : regexes)
      cout << "Dreme found motif " << x.first << " with E-value " << x.second
           << endl;
  } catch (exception& e) {
    cout << "Exception: " << e.what() << endl;
  }
  return EXIT_SUCCESS;
}
