#ifndef OPTIONS_HPP__
#define OPTIONS_HPP__
#include <string>
#include <vector>
#include <map>

class Options {
public:
  static int help;

  static int get(std::string opts, std::string* arg);
  static int get(std::string opts) {
    return get(opts, 0);
  }
  

  static void init(int argc, const char** argv);

  static std::string str();
  
  static bool inited;

  static std::vector<std::string> argv;
private:
  
  static std::string input;
  static std::map<std::string, std::string> opts_map;
  
};

#endif
