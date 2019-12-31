#ifndef PHYLONAUT_DISTANCEMATRIX_HPP__
#define PHYLONAUT_DISTANCEMATRIX_HPP__

#include "Clade.hpp"
#include "TaxonSet.hpp"
#include <unordered_set>
#include <vector>
#include <functional>

class DistanceMatrix {
private:
  const TaxonSet* ts;
  std::vector<double> d;
  std::vector<double> mask_;

public:
  DistanceMatrix(const TaxonSet& ts);
  DistanceMatrix(const TaxonSet& ts, std::string newick);
  DistanceMatrix(const DistanceMatrix& other) :
    ts(other.ts),
    d(other.d),
    mask_(other.mask_)
  {}

  double& get(Taxon t1, Taxon t2, std::vector<double>& myD);
  double get(Taxon t1, Taxon t2, const std::vector<double>& myD) const;

  double& operator()(Taxon t1, Taxon t2) {
    return get(t1, t2);
  }
  double& get(Taxon t1, Taxon t2) {
    return get(t1, t2, d);
  }

  double has(Taxon t1, Taxon t2) const {
    return get(t1, t2, mask_) != 0;
  }


  double isMasked(Taxon t1, Taxon t2) const {
    return get(t1, t2, mask_) == 0;
  }
  double& masked(Taxon t1, Taxon t2) {
    return get(t1, t2, mask_);
  }


  DistanceMatrix& operator=(const DistanceMatrix& other) {
    ts = other.ts;
    d = other.d;
    mask_ = other.mask_;
    return *this;
  }


  DistanceMatrix& operator+=(const DistanceMatrix& other) {
    for (size_t i= 0; i < d.size(); i++) {
      d[i] += other.d[i];
    }
    for (size_t i= 0; i < mask_.size(); i++) {
      mask_[i] += other.mask_[i];
    }
    return *this;
  }

  DistanceMatrix& operator*=(const double val) {
    for (size_t i= 0; i < d.size(); i++) {
      d[i] *= val;
    }
    for (size_t i= 0; i < mask_.size(); i++) {
      mask_[i] *= val;
    }
    return *this;
  }

  std::string str();


  std::unordered_set<Clade> upgma();

};


struct DisjointSet {
  std::vector<int> parent;
  std::vector<int> rank;
  std::vector<int> size;
  std::vector<Clade> clade;
  DisjointSet(int n, const TaxonSet& ts) :
  parent(n), rank(n, 1), size(n, 1), clade(n, ts) {
    for (int i = 0; i < n; i++) {
      parent[i] = i;

      clade[i].add(i);

    }
  }
  int find(int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  }

  void merge(int x, int y) {
    x = find(x);
    y = find(y);
    if (rank[x] > rank[y]) {
      parent[y] = x;
      size[x] += size[y];
    } else if (rank[y] > rank[x]) {
      parent[x] = y;
      size[y] = size[x] + size[y];

    } else {
      parent[x] = y;
      rank[x] ++;
      size[x] += size[y];

    };
    clade[x] = clade[x] + clade[y];
    clade[y] = clade[x];
  }
};

#endif
