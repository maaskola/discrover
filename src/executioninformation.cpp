
#include <cstdio>
#include <ctime>
#include <boost/filesystem.hpp>
#include "executioninformation.hpp"

std::string cmdline(int argc, const char** argv)
{
  std::string cmd;
  for(int i = 0; i < argc; i++)
    cmd += (i != 0 ? " " : "") + std::string(argv[i]);
  return(cmd);
}

ExecutionInformation generate_exec_info(const std::string &name, const std::string &hmm_version, const std::string &cmdline)
{
  time_t rawtime;
  time(&rawtime);

  std::string datetime = ctime(&rawtime);
  datetime = datetime.substr(0,datetime.size() - 1);

  std::string dir = boost::filesystem::initial_path().string();

  ExecutionInformation exec_info = {name, hmm_version, cmdline, datetime, dir};
  return(exec_info);
}


