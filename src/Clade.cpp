#include "Clade.hpp"
#include "util/Options.hpp"
#include "util/Logger.hpp"
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

Clade::Clade(TaxonSet& ts_, string& str) :
  taxa(ts_.size()),
  ts(ts_),
  sz(0)
{
  char* cladestr = &(str[1]);
  char* token;
  char* saveptr;

  while((token = strtok_r(cladestr, ",} ", &saveptr))) {
    cladestr = NULL;
    add(ts[string(&(token[0]))]);
  }
}

Clade::Clade(TaxonSet& ts_) :
  taxa(ts_.size()),
  ts(ts_),
  sz(0)
{}

Clade::Clade(TaxonSet& ts_, const clade_bitset& taxa) :
  taxa(taxa),
  ts(ts_),
  sz(taxa.popcount())
{
}

Clade::Clade(TaxonSet& ts_, const unordered_set<Taxon>& taxa) :
  taxa(ts_.size()),
  ts(ts_),
  sz(taxa.size())
{
  for (Taxon t : taxa) {
    add(t);
  }
}

Clade::Clade(const Clade& other) :
  taxa(other.taxa),
  ts(other.ts),
  sz(other.sz)
{
}



Clade& Clade::operator=(const Clade& other) {
   taxa = other.taxa;
   ts = other.ts;
   sz = other.sz;
   return *this;
}

bool Clade::operator==(const Clade& other) const {
   return taxa == other.taxa;
}





string Clade::str() const {
  stringstream ss;
  vector<string> strings;

  for (Taxon i : *this) {
    strings.push_back(ts[i]);
  }

  sort(strings.begin(), strings.end());
  
  ss << '{';
  for (string s : strings) {
    ss << s << ", ";
  }
  ss << '}';
  return ss.str();
}




void Clade::test() {
  string str = string("{tx1, tx8, tx3, tx2, tx4}");
  TaxonSet ts(str);
  cout << Clade(ts, str).str() << endl;
  cout << ts.str() << endl;
}

Clade Clade::overlap(const Clade& other) const {
  clade_bitset cb = other.taxa & taxa;
  return Clade(ts, cb);
}

int Clade::overlap_size(const Clade& other) const {
  return taxa.overlap_size(other.taxa);
}


bool Clade::contains(const Clade& other) const {
  if (other.size() > size())
    return false;
  
  bool status = 1;  
  for (size_t i = 0; i < taxa.cap; i++) {
    status &= ((other.taxa.data[i] & taxa.data[i]) == other.taxa.data[i]);
  }
  
  return status;
}
bool Clade::contains(const Taxon taxon) const {
  return taxa.get(taxon);
}

bool Clade::compatible(const Clade& other) const {
  return contains(other) || other.contains(*this) || overlap_size(other) == 0;
}

void Clade::add(const Taxon taxon) {
  taxa.set(taxon);
  sz++;
}

Clade Clade::complement() const {
  BitVectorFixed comp = ts.taxa_bs & (~taxa);
  Clade c(ts, comp);
  return c;
}

Clade Clade::minus(const Clade& other) const {
  BitVectorFixed m(taxa & (~other.taxa));
  Clade c(ts, m);
  return c;
}

int Clade::size() const {
  return sz;
}

void Clade::do_swap(Clade& other) {
  std::swap(taxa, other.taxa);
  std::swap(sz, other.sz);
}



string Bipartition::str() const {
  return "{" + a1.str() + " " + a2.str() + "}";
}

namespace std
{
    template<>
    void swap<Clade>(Clade& lhs, Clade& rhs)
    {
      lhs.do_swap(rhs);
    }
}
