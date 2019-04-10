#ifndef CLADE_HPP__
#define CLADE_HPP__

#include <string>
#include <vector>
#include <iostream>
#include <bitset>
#include <unordered_set>
#include <cassert>
#include <unordered_map>
#include <string.h>
#include "TaxonSet.hpp"


class Clade {
private:


  clade_bitset taxa;
  const TaxonSet* ts_;

  int sz;

public:
  
  Clade(TaxonSet& ts, string& str);
  Clade(const TaxonSet& ts, Taxon t);
  Clade(const TaxonSet& ts, const clade_bitset& taxa);
  Clade(const TaxonSet& ts, const unordered_set<Taxon>& taxa);
  Clade(const TaxonSet& ts);
  Clade(const Clade& other);

  Clade& operator=(const Clade& other);
  bool operator==(const Clade& other) const;


  string str() const;


  bool contains(const Clade& other) const;
  bool contains(const Taxon taxon) const;

  bool compatible(const Clade& other) const;
  bool compatible(const Clade& other, const Clade& restr) const;

  Clade overlap(const Clade& other) const;
  int overlap_size(const Clade& other) const ;

  static void test();

  void add(const Taxon taxon);
  void remove(const Taxon taxon);
  void add(const Clade& other);
  void remove(const Clade& other);
  Clade complement() const;
  Clade minus(const Clade& other) const;
  Clade plus(const Clade& other) const;
  Clade minus(const Taxon other) const;
  Clade plus(const Taxon other) const;

  Clade& operator-=(const Clade& other);
  Clade& operator+=(const Clade& other);
  Clade& operator-=(const Taxon other);
  Clade& operator+=(const Taxon other);


  Clade operator-(const Clade& other) const;
  Clade operator+(const Clade& other) const;
  Clade operator-(const Taxon other) const;
  Clade operator+(const Taxon other) const;


  const TaxonSet& ts() const {return *ts_;}
  int size() const;
  const clade_bitset& get_taxa() const {return taxa;}



  BVFIterator begin() const {
    return taxa.begin();
  }

  BVFIterator end() const {
    return taxa.end();
  }

  void do_swap(Clade& other);
  size_t hash() const { return taxa.hash(); }
};

ostream& operator<<(ostream& os, const Clade& c);

template<class c>
struct TripartitionG {
  c a1, a2, rest;
  TripartitionG(TaxonSet& ts, c& clade, c& subclade) :
    a1(clade.minus(subclade)),
    a2(subclade),
    rest(clade.complement()){ }

  string str() const  {
    assert(a1.overlap(rest).size() == 0);
    assert(a2.overlap(rest).size() == 0);
    assert(a2.overlap(a1).size() == 0);
    return "{" + a1.str() + "/" + a2.str() + "/" + rest.str() + "}";
  }
};

typedef TripartitionG<Clade> Tripartition;

struct Bipartition {
  Clade a1, a2;
  Bipartition(const Clade& clade1, const Clade& clade2) :
    a1(clade1),
    a2(clade2)
  {}
  size_t hash() const { return a1.hash() ^ a2.hash(); }
  bool operator==(const Bipartition& other) const {
    return ((a1 == other.a1) && (a2 == other.a2)) || ((a2 == other.a1) && (a1 == other.a2));
  }
  string str() const;
};


namespace std {
  template <> struct hash<Clade> {
    size_t operator()(const Clade& bvf) const {
      return bvf.hash() ;
    }
  };

  template <> struct hash<Bipartition> {
    size_t operator()(const Bipartition& bp) const {
      return bp.hash();
    }
  };
}

#endif // CLADE_HPP__
