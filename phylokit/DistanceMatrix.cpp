#include "DistanceMatrix.hpp"
#include <queue>
#include <sstream>
#include <glog/logging.h>
#include <vector>

#include <boost/tokenizer.hpp>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

DistanceMatrix::DistanceMatrix(const TaxonSet &ts) :
    ts(&ts) {
  d.resize(ts.size() * ts.size(), 0);
  mask_.resize(ts.size() * ts.size(), 0);
}

DistanceMatrix::DistanceMatrix(const TaxonSet &ts, std::string newick) :
    ts(&ts) {
  d.resize(ts.size() * ts.size(), 0);
  mask_.resize(ts.size() * ts.size(), 0);

  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(newick, sep);

  std::vector<double> dists(ts.size(), 0);
  std::vector<double> ops(ts.size(), 0);

  std::vector<Taxon> seen;
  std::string prevtok = "";

  for (auto tok : tokens) {

    if (tok == "(") {
      for (Taxon s : seen) {
        ops[s] += 1;
        dists[s] += 1;
      }
    } else if (tok == ")") {
      for (Taxon s : seen) {
        if (ops[s]) {
          dists[s] -= 1;
          ops[s] -= 1;
        } else {
          dists[s] += 1;
        }
      }
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        continue;
      }
      boost::algorithm::trim(tok);
      Taxon id = ts[tok];
      for (Taxon other : seen) {
        get(other, id) += dists[other] + 2;
        get(other, id, mask_) += 1;
      }
      seen.push_back(id);
    }
    prevtok = tok;
  }
}

std::string DistanceMatrix::str() {
  std::stringstream ss;
  for (Taxon t1 : *ts) {
    for (Taxon t2 : *ts) {
      ss << get(t1, t2) << ";" << masked(t1, t2) << "\t";
    }
    ss << std::endl;
  }
  return ss.str();
}


std::ostream& DistanceMatrix::writePhylip(std::ostream& out) {
  out << ts->size();
  out << std::endl;
  for (Taxon t1 : *ts) {
    for (Taxon t2 : *ts) {
      out << get(t1, t2) << " ";
    }
    out << std::endl;
  }
  return out;
}


double &DistanceMatrix::get(Taxon t1, Taxon t2, std::vector<double> &myD) {
  Taxon a = std::min(t1, t2);
  Taxon b = std::max(t1, t2);

  return myD[(b * (b + 1)) / 2 + a];
};

double DistanceMatrix::get(Taxon t1, Taxon t2, const std::vector<double> &myD) const {
  Taxon a = std::min(t1, t2);
  Taxon b = std::max(t1, t2);

  return myD[(b * (b + 1)) / 2 + a];
};

std::unordered_set<Clade> DistanceMatrix::upgma() {
  DLOG(INFO) << "Running UPGMA\n";
  std::vector<double> myD(d);
  std::vector<double> myMask(mask_);

  std::unordered_set<Clade> clades;
  DisjointSet sets(ts->size(), *ts);
  std::priority_queue<std::tuple<double, Taxon, Taxon, int, int>,
                      std::vector<std::tuple<double, Taxon, Taxon, int, int>>,
                      std::greater<std::tuple<double, Taxon, Taxon, int, int>>> pq;

  for (Taxon t1 : ts->taxa_bs) {
    Clade c(*ts);
    c.add(t1);
    clades.insert(c);
    for (Taxon t2 : ts->taxa_bs) {
      if (t1 < t2 && !isMasked(t1, t2))
        pq.push(std::make_tuple(get(t1, t2), t1, t2, 1, 1));
    }
  }

  double dist;
  Taxon t1, t2;
  int s1, s2;

  while (pq.size()) {
    std::tie(dist, t1, t2, s1, s2) = pq.top();
    pq.pop();

    if (sets.find(t1) != t1 || sets.find(t2) != t2)
      continue;
    if (sets.size[t1] != s1 || sets.size[t2] != s2)
      continue;

    sets.merge(t1, t2);

    clades.insert(sets.clade[sets.find(t1)]);

    for (Taxon t : *ts) {
      if (sets.find(t) == t && !masked(t, sets.find(t1))) {

        get(sets.find(t1), t, myD) = (get(t1, t, myMask) * get(t1, t, myD) + get(t2, t, myMask) * get(t2, t, myD))
            / ((get(t1, t, myMask) + get(t2, t, myMask)));
        masked(sets.find(t1), t) = 1;

        pq.push(std::make_tuple(get(sets.find(t1), t, myD), sets.find(t1), t, sets.size[sets.find(t1)], sets.size[t]));

      }
    }
  }

  return clades;
}
