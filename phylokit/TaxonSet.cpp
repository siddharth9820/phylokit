#include "TaxonSet.hpp"
#include <cassert>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

#define strtok_r strtok_s
TaxonSet::TaxonSet(int size) : frozen(false), taxa_bs(size) {}

TaxonSet::TaxonSet(std::string str)
    : frozen(false), taxa_bs(resize_clades(str)) {
  taxa.reserve(taxa_set.size());

  for (const std::string &i : taxa_set) {
    add(i);
  }
}

TaxonSet::TaxonSet(const TaxonSet &&other)
    : taxa_set(other.taxa_set),
      taxa(other.taxa),
      index(other.index),
      frozen(other.frozen),
      taxa_bs(other.taxa_bs) {}
TaxonSet &TaxonSet::operator=(const TaxonSet &&other) {
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

void TaxonSet::freeze() { frozen = true; }

int TaxonSet::resize_clades(std::string str) {
  std::stringstream stream(str);
  std::string s;

  while (!stream.eof()) {
    std::getline(stream, s);
    if (s.size() == 0) continue;
    add_clade_taxa(s, taxa_set);
  }
  return taxa_set.size();
}

void TaxonSet::add_clade_taxa(std::string str,
                              std::unordered_set<std::string> &taxa_set) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep("{,}");
  tokenizer tokens(str, sep);

  for (auto &token : tokens) {
    taxa_set.insert(token);
  }
}

Taxon TaxonSet::add(const std::string &str) {
  if (index.count(str)) {
    return index[str];
  }

  if (frozen) {
    std::cerr << "Trying to add " << str << " to frozen taxon set\n";
    for (const std::string &i : taxa_set) {
      std::cerr << (i) << std::endl;
    }
    assert(false);
  }

  int i = taxa.size();
  taxa.push_back(str);
  index[str] = i;
  taxa_bs.set(i);
  return i;
}

size_t TaxonSet::size() const { return taxa.size(); }
