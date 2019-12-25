#include "Clade.hpp"
#include "util/Logger.hpp"
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>


#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>


Clade::Clade(TaxonSet& ts_, std::string& str) :
  taxa(ts_.size()),
  ts_(&ts_),
  sz(0)
{

  typedef boost::tokenizer<boost::char_separator<char> >
    tokenizer;
  boost::char_separator<char> sep("{,}");
  tokenizer tokens(str, sep);

  for (auto& token : tokens) {
    add(ts_[token]);
  }
  
}

Clade::Clade(const TaxonSet& ts_) :
  taxa(ts_.size()),
  ts_(&ts_),
  sz(0)
{}


Clade::Clade(const TaxonSet& ts_, Taxon t) :
  taxa(ts_.size()),
  ts_(&ts_),
  sz(0)
{
  add(t);
}

Clade::Clade(const TaxonSet& ts_, const clade_bitset& taxa) :
  taxa(taxa),
  ts_(&ts_),
  sz(taxa.popcount())
{
}

Clade::Clade(const TaxonSet& ts_, const unordered_set<Taxon>& taxa) :
  taxa(ts_.size()),
  ts_(&ts_),
  sz(taxa.size())
{
  for (Taxon t : taxa) {
    add(t);
  }
}

Clade::Clade(const Clade& other) :
  taxa(other.taxa),
  ts_(&other.ts()),
  sz(other.sz)
{
}



Clade& Clade::operator=(const Clade& other) {
   taxa = other.taxa;
   ts_ = other.ts_;
   sz = other.sz;
   return *this;
}

bool Clade::operator==(const Clade& other) const {
   return taxa == other.taxa;
}





std::string Clade::str() const {
  std::stringstream ss;
  std::vector<std::string> strings;

  for (Taxon i : *this) {
    strings.push_back(ts()[i]);
  }

  std::sort(strings.begin(), strings.end());

  ss << '{';
  for (string s : strings) {
    ss << s << ", ";
  }
  ss << '}';
  return ss.str();
}




void Clade::test() {
  std::string str = std::string("{tx1, tx8, tx3, tx2, tx4}");
  TaxonSet ts(str);
  cout << Clade(ts, str).str() << endl;
  cout << ts.str() << endl;
}

Clade Clade::overlap(const Clade& other) const {
  clade_bitset cb = other.taxa & taxa;
  return Clade(ts(), cb);
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

bool Clade::compatible(const Clade& other, const Clade& restr) const {
  return contains(other) || other.contains(*this) || overlap_size(other) == 0;
}

void Clade::add(const Taxon taxon) {
  taxa.set(taxon);
  sz++;
}

void Clade::remove(const Taxon taxon) {
  taxa.unset(taxon);
  sz--;
}


void Clade::add(const Clade& other) {
  taxa |= other.taxa;
  sz = taxa.popcount();
}

void Clade::remove(const Clade& other) {
  taxa &= ~other.taxa;
  sz = taxa.popcount();

}


Clade Clade::complement() const {
  BitVectorFixed comp = ts().taxa_bs & (~taxa);
  Clade c(ts(), comp);
  return c;
}

Clade Clade::minus(const Clade& other) const {
  BitVectorFixed m(taxa & (~other.taxa));
  Clade c(ts(), m);
  return c;
}
Clade Clade::plus(const Clade& other) const {
  BitVectorFixed m(taxa | other.taxa);
  Clade c(ts(), m);
  return c;
}

Clade Clade::minus(const Taxon other) const {
  Clade c(ts(), taxa);
  c.remove(other);
  return c;
}
Clade Clade::plus(const Taxon other) const {
  Clade c(ts(), taxa);
  c.add(other);
  return c;
}


Clade Clade::operator+(const Clade& other) const {
  return plus(other);
}

Clade Clade::operator-(const Clade& other) const {
  return minus(other);
}



Clade Clade::operator+(const Taxon other) const {
  return plus(other);
}

Clade Clade::operator-(const Taxon other) const {
  return minus(other);
}


Clade& Clade::operator-=(const Clade& other){
  remove(other);
  return *this;
}
Clade& Clade::operator+=(const Clade& other){
  add(other);
  return *this;
}
Clade& Clade::operator-=(const Taxon other){
  remove(other);
  return *this;
}
Clade& Clade::operator+=(const Taxon other){
  add(other);
  return *this;
}


ostream& operator<<(ostream& os, const Clade& c) {
  os << c.str();
  return os;
}

int Clade::size() const {
  return sz;
}

void Clade::do_swap(Clade& other) {
  std::swap(taxa, other.taxa);
  std::swap(sz, other.sz);
}



std::string Bipartition::str() const {
  return "{" + a1.str() + " " + a2.str() + "}";
}



