#include "TaxonSet.hpp"
#include <iostream>
#include <cstring>
#include <cassert>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#define strtok_r strtok_s
TaxonSet::TaxonSet(int size):
  frozen(false),
  taxa_bs(size) {
}

TaxonSet::TaxonSet(string str):
  frozen(false),
  taxa_bs(resize_clades(str))
{
  taxa.reserve(taxa_set.size());

  for (const string& i : taxa_set) {
    add(i);
  }

}

TaxonSet::TaxonSet(const TaxonSet&& other) :
  taxa_set(other.taxa_set),
  taxa(other.taxa),
  index(other.index),
  frozen(other.frozen),
  taxa_bs(other.taxa_bs)
{}
TaxonSet& TaxonSet::operator=(const TaxonSet&& other) {
  if (this == &other) {
    return *this;
  }
  taxa_set = other.taxa_set;
  taxa = other.taxa;
  index = other.index;
  frozen = other.frozen;
  taxa_bs = other.taxa_bs;
  return *this;
}

void TaxonSet::freeze() {
  frozen = true;
}

int TaxonSet::resize_clades(string str) {
  stringstream stream(str);
  string s;

  while(!stream.eof()) {
    getline(stream, s);
    if (s.size() == 0)
      continue;
    add_clade_taxa(s, taxa_set);
  }
  return taxa_set.size();
}

void TaxonSet::add_clade_taxa(string str, unordered_set<string>& taxa_set) {

  typedef boost::tokenizer<boost::char_separator<char> >
    tokenizer;
  boost::char_separator<char> sep("{,}");
  tokenizer tokens(str, sep);
  
  for (auto& token : tokens) {
    taxa_set.insert(token);
  }

}

Taxon TaxonSet::add(const string& str) {
  if (index.count(str)) {
    return index[str];
  }

  if (frozen) {
    cerr << "Trying to add " << str << " to frozen taxon set\n";
    for (const string& i : taxa_set) {
      cerr << (i) << endl;
    }
    assert(false);
  }


  int i = taxa.size();
  taxa.push_back(str);
  index[str] = i;
  taxa_bs.set(i);
  return i;
}

size_t TaxonSet::size() const {
  return taxa.size();
}
