#include "DistanceMatrix.hpp"
#include <queue>
#include "util/Logger.hpp"
#include <sstream>


#include <boost/tokenizer.hpp>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>


DistanceMatrix::DistanceMatrix(const TaxonSet& ts) :
  ts(&ts) {
    d.resize(ts.size() * ts.size(), 0);
    mask_.resize(ts.size() * ts.size(), 0);
}

DistanceMatrix::DistanceMatrix(const TaxonSet& ts, string newick) :
ts(&ts) {
  d.resize(ts.size() * ts.size(), 0);
  mask_.resize(ts.size() * ts.size(), 0);

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(newick, sep);

  vector<double> dists(ts.size(), 0);
  vector<double> ops(ts.size(), 0);

  vector<Taxon> seen;
  string prevtok = "";

  for (auto tok : tokens) {

    if (tok == "(") {
      for (Taxon s : seen) {
	       ops[s] += 1;
	        dists[s] += 1;
      }
    }
    else if (tok == ")") {
      for (Taxon s : seen) {
	if (ops[s]) {
	  dists[s] -= 1;
	  ops[s] -= 1;
	} else {
	  dists[s] += 1;
	}
      }
    }
    else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" or prevtok == ":" or (tok == " " and prevtok == ",")) {
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


string DistanceMatrix::str() {
  stringstream ss;
  for (Taxon t1 : *ts) {
    for (Taxon t2 : *ts) {
      ss << get(t1, t2) << ";" << masked(t1, t2) <<  "\t";
    }
    ss << endl;
  }
  return ss.str();
}

double& DistanceMatrix::get(Taxon t1, Taxon t2, vector<double>& myD) {
  Taxon a = min(t1, t2);
  Taxon b = max(t1, t2);

  return myD[(b * (b+1))/2 + a];
};


double DistanceMatrix::get(Taxon t1, Taxon t2, const vector<double>& myD) const {
  Taxon a = min(t1, t2);
  Taxon b = max(t1, t2);

  return myD[(b * (b+1))/2 + a];
};


unordered_set<Clade> DistanceMatrix::upgma() {
  PROGRESS << "Running UPGMA\n";
  vector<double> myD(d);
  vector<double> myMask(mask_);

  unordered_set<Clade> clades;
  DisjointSet sets(ts->size(), *ts);
  priority_queue<tuple<double, Taxon, Taxon, int, int>, vector<tuple<double, Taxon, Taxon, int, int>>, greater<tuple<double, Taxon, Taxon, int, int>>> pq;

  for (Taxon t1 : ts->taxa_bs) {
    Clade c(*ts);
    c.add(t1);
    clades.insert(c);
    for (Taxon t2 : ts->taxa_bs) {
      if (t1 < t2 && !isMasked(t1, t2))
        pq.push(make_tuple(get(t1, t2), t1, t2, 1, 1));
    }
  }

  double dist;
  Taxon t1, t2;
  int s1, s2;

  while(pq.size()) {
    tie(dist, t1, t2, s1, s2) = pq.top();
    pq.pop();

    if (sets.find(t1) != t1 || sets.find(t2) != t2)
      continue;
    if (sets.size[t1] != s1 || sets.size[t2] != s2)
      continue;

    sets.merge(t1, t2);

    clades.insert(sets.clade[sets.find(t1)]);

    for (Taxon t : *ts) {
      if (sets.find(t) == t && !masked(t, sets.find(t1)) ) {

        get(sets.find(t1), t, myD) = (get(t1, t, myMask) * get(t1, t, myD) + get(t2, t, myMask) * get(t2, t, myD))/((get(t1, t, myMask) + get(t2, t, myMask)));
        masked(sets.find(t1), t) = 1;

        pq.push(make_tuple(get(sets.find(t1), t, myD), sets.find(t1), t, sets.size[sets.find(t1)], sets.size[t]));

      }
    }
  }

  return clades;
}
