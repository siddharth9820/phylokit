#include "Options.hpp"
#include <cassert>
#include <iostream>
#include <iterator>
#include <sstream>

bool Options::inited = false;

std::vector<std::string> Options::argv;
std::map<std::string, std::string> Options::opts_map;
std::string Options::input;

enum option_type {SHORT, LONG, ARG, END, EMPTY};

option_type get_option_type(std::string& arg) {
  if (arg.size() == 0) return EMPTY;
  if (arg[0] != '-') return ARG;
  if (arg.size() == 1) {
    std::cerr << "INVALID ARGUMENT " << arg << std::endl;
    exit(1);
  }
  if (arg[1] == '-') {
    if (arg.size() == 2) return END;
    arg = std::string(arg, 2); //remove starting --
    return LONG;
  }
  if (arg.size() > 2) {
    std::cerr << "INVALID ARGUMENT " << arg << std::endl;
    exit(1);
  }
  arg = std::string(arg, 1);
  return SHORT;
}

std::string Options::str() {
  return input;
}

void Options::init(int argc_, const char** argv_) {

  for (int i = 1; i < argc_; i++) {
    //    cerr << i << " " << argv_[i] << endl;
    argv.push_back(std::string(argv_[i]));
    input += std::string(argv_[i]) + " ";
  }

  argv.push_back("--");

  std::string last_option = "";
  
  for (std::string arg : argv) {
    option_type opttype = get_option_type(arg);
    //    DEBUG << arg << " " <<  opttype << endl;
    switch(opttype) {
      case SHORT:
      case LONG:
        if (last_option != "") {
          opts_map[last_option] = "";
        }
        last_option = arg;
        break;
      case ARG:
        if (last_option == ""){
          std::cerr << "ARGUMENT WITHOUT OPTION: " << arg << std::endl;
          exit(1);
        }
        opts_map[last_option] = arg;
        last_option = "";
        break;
      case END:
        if (last_option != "") {
          opts_map[last_option] = "";
        }
        break;
      case EMPTY:
        break;
    }

  }

  // DEBUG << "OPTIONS MAP:\n";
  // for (auto& kv : opts_map) {
  //   DEBUG << kv.first << " = " << kv.second << endl;
  // }
  
  inited = true;
}


int Options::get(std::string opts, std::string* arg) {
  assert(inited);

  std::stringstream ss(opts);
  std::istream_iterator<std::string> begin(ss);
  std::istream_iterator<std::string> end;
  std::vector<std::string> vopts(begin, end);
  
  for (auto& opt : vopts) {
    if (opts_map.count(opt)) {
      if (arg)
	*arg = opts_map[opt];

      return 1;
    }

  }
  return 0;
  
  // opterr = 0;
  
  // struct option opts[] = {{ longopt.c_str(), optional_argument, 0, 0},
  // 			  { 0, 0, 0, 0}
  // };
  // int c;
  // int option_index;
  // string optstring = "-:";
  // if (shortopt) {
  //   if (arg)
  //     optstring = ":" + string(&shortopt, 1) + ":";
  //   else
  //     optstring = ":" + string(&shortopt, 1) ;
  // }
  // while ((c = getopt_long(argc, argv, optstring.c_str(), &(opts[0]), &option_index)) != -1) {
  //   if (c == shortopt || c == 0 ) {
  //     if (arg && (size_t)optarg)
  // 	*arg = string(optarg);
  //     optind = 1;
  //     return 1;
  //     cout << shortopt << " " << longopt << " " << arg << endl;
  //   }
  // }
  // optind = 1;
  // return 0;  
}

