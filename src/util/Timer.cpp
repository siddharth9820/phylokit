#include "Timer.hpp"

#include <cassert>

vector<pair<string, milliseconds> > Timer::times;
unordered_map<string, steady_clock::time_point> Timer::starts;

void Timer::start(string tag) {
  assert(starts.count(tag) == 0);
  starts[tag] = steady_clock::now();
}

void Timer::stop(string tag) {
  assert(starts.count(tag) == 1);
  times.push_back(make_pair(tag,  duration_cast<milliseconds>(steady_clock::now() - starts[tag])));
}


void Timer::writeAll(ostream& out) {
  for (auto p : times ) {
    out << p.first << "\t" << p.second.count() << endl;
  }
}

void Timer::reset() {
  starts.clear();
  times.clear();
}
