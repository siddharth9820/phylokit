#include "newick.hpp"


int newick_to_ts(const string& s, unordered_set<string>& taxa) {
  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";\n", "(),:");
  tokenizer tokens(s, sep);

  int taxon_count = 0;

  string prevtok = "";
  
  for (auto tok : tokens) {
    boost::algorithm::trim(tok);


    if (tok == ":" || tok == "," || tok == "(" || tok == ")") {
    } else {
      if (prevtok == ")" or prevtok == ":") {	
	continue;
      }
      if(tok.find_first_not_of(' ') != string::npos) {
	taxa.insert(tok);
	taxon_count ++;
      }
    }
    prevtok = tok;
  }
  return taxon_count;
}

void newick_to_dm(const string& s, TaxonSet& ts, dm_type& dist_mat, dm_type& mask_mat ) {
  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";\n", "():,");
  
  tokenizer tokens(s, sep);
  
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
	dist_mat[other][id] += dists[other] + 2;
	dist_mat[id][other] += dists[other] + 2;
	mask_mat[other][id] += 1;
	mask_mat[id][other] += 1;	
      }
      seen.push_back(id);
    }
    prevtok = tok;
  }
}

void newick_to_clades(const string& s, TaxonSet& ts, unordered_set<Clade>& clade_set) {
  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";\n", "():,");
  
  tokenizer tokens(s, sep);

  vector<size_t> active;
  vector<Clade> clades;
  
  string prevtok = "";
  
  for (auto tok : tokens) {
    
    if (tok == "(") {
      clades.emplace_back(ts);
      active.push_back(clades.size() - 1);
    }
    else if (tok == ")") {
      active.pop_back();
    }
    else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" or prevtok == ":" or (tok == " " and prevtok == ",")) {
	continue;
      }
      boost::algorithm::trim(tok);
      Taxon id = ts[tok];

      for (size_t a : active) {
	clades.at(a).add(id);
      }

    }
    prevtok = tok;
  }
  for (Clade& c : clades) {
    clade_set.insert(c);
  }
}

string map_newick_names(const string& s, TaxonSet& ts) {
    typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";\n", "():,");
  
  tokenizer tokens(s, sep);

  stringstream output;
  
  string prevtok = "";
  
  for (auto tok : tokens) {
    
    if (tok == "(") {
      output << "(";
    }
    else if (tok == ")") {
      output << ")";
    }
    else if (tok == ":") {
    } else if (tok == ",") {
      output << ",";
    } else {
      if (prevtok == ")" or prevtok == ":" or (tok == " " and prevtok == ",")) {
	continue;
      }
      boost::algorithm::trim(tok);
      Taxon id = ts[tok];
      output << id;
      
    }
    prevtok = tok;
  }
  output << ";";
  return output.str();
}


string unmap_newick_names(const string& s, TaxonSet& ts) {
  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";\n", "():,");
  
  tokenizer tokens(s, sep);

  stringstream output;
  
  string prevtok = "";
  
  for (auto tok : tokens) {
    
    if (tok == "(") {
      output << "(";
    }
    else if (tok == ")") {
      output << ")";
    }
    else if (tok == ":") {
    } else if (tok == ",") {
      output << ",";
    } else {
      if (prevtok == ")" or prevtok == ":" or (tok == " " and prevtok == ",")) {
	continue;
      }
      boost::algorithm::trim(tok);
      stringstream ss(tok);
      int i;
      ss >> i;
      const string& id = ts[i];
      output << id;
      
    }
    prevtok = tok;
  }
  output << ";";
  return output.str();
}


string unmap_clade_names(const string& s, TaxonSet& ts) {
  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";\n", "{},");
  
  tokenizer tokens(s, sep);

  stringstream output;
  
  string prevtok = "";
  
  for (auto tok : tokens) {
    
    if (tok == "{") {
      output << "{";
    }
    else if (tok == "}") {
      output << "}";
    }
    else if (tok == ",") {
      output << ",";
    } else {
      if (prevtok == "}" or prevtok == ":" or (tok == " " and prevtok == ",")) {
	continue;
      }
      boost::algorithm::trim(tok);
      stringstream ss(tok);
      int i;
      ss >> i;
      const string& id = ts[i];
      output << id;
      
    }
    prevtok = tok;
  }
  return output.str();
}
