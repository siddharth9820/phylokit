#ifndef __NEWICK_HPP__
#define __NEWICK_HPP__

#include <boost/tokenizer.hpp>
#include <boost/multi_array.hpp>
#include <boost/algorithm/string.hpp>

#include <string>
#include <unordered_set>
#include <iostream>

#include "Clade.hpp"
#include "TreeClade.hpp"
#include "TaxonSet.hpp"

typedef boost::multi_array<double, 2> dm_type;

int newick_to_ts(const std::string& s, std::unordered_set<std::string>& taxa);


Clade newick_to_taxa(const std::string& s, TaxonSet& ts);
void newick_to_dm(const std::string& s, TaxonSet& ts, dm_type& dist_mat, dm_type& mask_mat );
void newick_to_clades(const std::string& s, TaxonSet& ts, std::unordered_set<Clade>& clade_set);
Tree newick_to_treeclades(const std::string& s, TaxonSet& ts);
void newick_to_postorder(const std::string& s, TaxonSet& ts, std::vector<Taxon>& order);


std::string map_newick_names(const std::string& s, TaxonSet& ts);
std::string unmap_newick_names(const std::string& s, TaxonSet& ts);
std::string unmap_clade_names(const std::string& s, TaxonSet& ts);
#endif
