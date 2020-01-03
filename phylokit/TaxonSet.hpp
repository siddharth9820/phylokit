#ifndef TAXONSET_HPP__
#define TAXONSET_HPP__

#include <algorithm>
#include <bitset>
#include <map>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "BitVector.hpp"

typedef int Taxon;
typedef BitVectorFixed clade_bitset;
class Clade;

class TaxonSet {
 private:
  std::unordered_set<std::string> taxa_set;
  std::vector<std::string> taxa;
  std::map<std::string, Taxon> index;
  bool frozen;
  TaxonSet(const TaxonSet &other);
  TaxonSet &operator=(const TaxonSet &other);

 public:
  clade_bitset taxa_bs;

  TaxonSet(int size);
  TaxonSet(std::string str);

  TaxonSet(const TaxonSet &&other);
  TaxonSet &operator=(const TaxonSet &&other);

  BVFIterator begin() const { return taxa_bs.begin(); }

  BVFIterator end() const { return taxa_bs.end(); }

  int resize_clades(std::string str);

  void add_clade_taxa(std::string str,
                      std::unordered_set<std::string> &taxa_set);

  void freeze();

  Taxon operator[](const std::string &str) { return add(str); }

  Taxon operator[](const std::string &str) const { return index.at(str); }
  const std::string &operator[](const Taxon i) const { return taxa.at(i); }

  const std::string &get(const Taxon i) const { return taxa.at(i); }

  bool has(const std::string &str) const { return index.count(str) > 0; }

  size_t size() const;
  Taxon add(const std::string &str);
  std::string str() const {
    std::stringstream ss;
    for (size_t i = 0; i < taxa.size(); i++) {
      ss << i << "\t" << taxa[i] << std::endl;
    }
    return ss.str();
  }
  std::vector<std::string> &sort_taxa() {
    std::sort(taxa.begin(), taxa.end());
    return taxa;
  }
};

#endif  // TAXONSET_HPP__
