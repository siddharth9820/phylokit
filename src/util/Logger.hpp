#ifndef __LOGGER_HPP__
#define __LOGGER_HPP__

#include <unordered_set>
#include <iostream>

using namespace std;

class NullStream : public ostream {
public:
    void setFile() { /* no-op */ }
    template<typename TPrintable>
    NullStream& operator<<(TPrintable const&)
  { return *this;/* no-op */ }
};

class Logger {
public:
  static Logger& get();
  static ostream& log(string tag, string fname, int lineno);
  static void enable(string tag);
  static bool isEnabled(string tag);
  static void enable(string tag, Logger& logger);

  static void disable(string tag);
  static void setLevel(Logger& logger);
private:
  Logger();
  static unordered_set<string> enabled_tags;
  static unordered_set<string> enabled_files;
  static unordered_set<string> disabled_tags;
  static unordered_set<string> disabled_files;
  static int ilevel;
  static void getIlevel(string& level);
  static NullStream nstream;

};

#define LOG(tag) (Logger::isEnabled(tag)) && (Logger::log((tag), (__BASE_FILE__), (__LINE__)))

#ifndef DEBUG
#define DEBUG LOG("DEBUG")
#endif

#ifndef ERR
#define ERR LOG("ERROR")
#endif

#ifndef WARN
#define WARN LOG("WARNING")
#endif

#ifndef INFO
#define INFO LOG("INFO")
#endif

#ifndef PROGRESS
#define PROGRESS LOG("PROGRESS")
#endif

#endif
