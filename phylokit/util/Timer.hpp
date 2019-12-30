#ifndef __TIMER_HPP__
#define __TIMER_HPP__

#include <set>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

class Timer {
public:
  static void start(string tag);
  static void stop(string tag);
  static void writeAll(ostream& out);
  static void reset();
private:
  static vector<pair<string, milliseconds> > times; 
  static unordered_map<string, steady_clock::time_point> starts; 
};

#endif
