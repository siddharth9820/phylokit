#include "newick.hpp"
#include "TreeClade.hpp"

int newick_to_ts(const std::string &s, std::unordered_set<std::string> &taxa) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "(),:");
  tokenizer tokens(s, sep);

  int taxon_count = 0;

  std::string prevtok = "";

  for (auto tok : tokens) {
    boost::algorithm::trim(tok);

    if (tok == ":" || tok == "," || tok == "(" || tok == ")") {
    } else {
      if ((prevtok == ")") || (prevtok == ":")) {
        continue;
      }
      if (tok.find_first_not_of(' ') != std::string::npos) {
        taxa.insert(tok);
        taxon_count++;
      }
    }
    prevtok = tok;
  }
  return taxon_count;
}

Clade newick_to_taxa(const std::string &s, TaxonSet &ts) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  Clade clade(ts);

  tokenizer tokens(s, sep);

  std::string prevtok = "";

  for (auto tok : tokens) {
    if (tok == "(") {
    } else if (tok == ")") {
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if ((prevtok == ")") || (prevtok == ":") ||
          (tok == " " && prevtok == ",")) {
        continue;
      }
      boost::algorithm::trim(tok);
      Taxon id = ts[tok];

      clade.add(id);
    }
    prevtok = tok;
  }
  return clade;
}

void newick_to_clades(const std::string &s, TaxonSet &ts,
                      std::unordered_set<Clade> &clade_set) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(s, sep);

  std::vector<size_t> active;
  std::vector<Clade> clades;

  std::string prevtok = "";

  for (auto tok : tokens) {
    if (tok == "(") {
      clades.emplace_back(ts);
      active.push_back(clades.size() - 1);
    } else if (tok == ")") {
      active.pop_back();
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if ((prevtok == ")") || (prevtok == ":") ||
          (tok == " " && prevtok == ",")) {
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
  for (Clade &c : clades) {
    clade_set.insert(c);
  }
}

Tree newick_to_treeclades(const std::string &s, TaxonSet &ts) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(s, sep);

  std::vector<size_t> active;
  std::vector<std::vector<size_t>> children;

  Tree tree(ts);

  std::string prevtok = "";

  for (auto tok : tokens) {
    if (tok == "(") {
      int ind = tree.addNode();

      if (active.size()) {
        tree.node(active.back()).addChild(ind);
      }
      active.push_back(ind);
    } else if (tok == ")") {
      active.pop_back();
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        continue;
      }
      boost::algorithm::trim(tok);
      Taxon id = ts[tok];

      int ind = tree.addNode();

      if (active.size()) {
        tree.node(active.back()).addChild(ind);
      }
      tree.node(ind).add(id);

      for (size_t a : active) {
        tree.node(a).add(id);
      }
    }
    prevtok = tok;
  }
  return tree;
}

void newick_to_postorder(const std::string &s, TaxonSet &ts,
                         std::vector<Taxon> &order) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");
  std::string prevtok = "";

  tokenizer tokens(s, sep);

  std::vector<int> sizes;
  sizes.push_back(0);

  for (auto tok : tokens) {
    if (tok == "(") {
      sizes.back()++;
      sizes.push_back(0);
    } else if (tok == ")") {
      order.push_back(-1 * sizes.back());
      sizes.pop_back();
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        continue;
      }
      sizes.back()++;
      boost::algorithm::trim(tok);
      Taxon id = ts[tok];
      order.push_back(id);
    }
    prevtok = tok;
  }
}

bool is_rooted(const std::string& newick) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  std::string prevtok = "";
  tokenizer tokens(newick, sep);

  int paren_count = 0;
  int root_child_count = 0;
  for (auto tok: tokens) {
    if (tok == "(") {
      if (paren_count == 1) { // at root level
        root_child_count++;
      }
      paren_count++;
    } else if (tok == ")") {
      paren_count--;
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        continue;
      }
      if (paren_count == 1) { //lone taxa attached to root
        root_child_count++;
      }
    }
    prevtok = tok;
  }
  return root_child_count == 2;
}
std::string deroot(const std::string& newick) {
  if (!is_rooted(newick)) {
    return newick;
  }

  std::stringstream outputtree;
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  std::string prevtok = "";
  tokenizer tokens(newick, sep);

  int paren_count = 0;
  int root_child_count = 0;
  int big_root_child_count = 0;

  for (auto tok: tokens) {
    bool copy_token = true;
    if (tok == "(") {
      if (paren_count == 1) { // at root level
        root_child_count++;
        big_root_child_count++;
        if (big_root_child_count == 1) {
          copy_token = false;
        }
      }
      paren_count++;
    } else if (tok == ")") {
      paren_count--;
      if (paren_count == 1 && big_root_child_count == 1) {
        copy_token = false;
      }
    } else if (tok == ":") {
    } else if (tok == ",") {
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        outputtree << tok;
        continue;
      }
      if (paren_count == 1) { //lone taxa attached to root
        root_child_count++;
      }
    }
    if (copy_token) {
      outputtree << tok;
    }
    prevtok = tok;
  }
  return outputtree.str();
}

std::string map_newick_names(const std::string &s, TaxonSet &ts) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(s, sep);

  std::stringstream output;

  std::string prevtok = "";

  for (auto tok : tokens) {
    if (tok == "(") {
      output << "(";
    } else if (tok == ")") {
      output << ")";
    } else if (tok == ":") {
      output << ":";
    } else if (tok == ",") {
      output << ",";
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        output << tok;
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

std::string unmap_newick_names(const std::string &s, TaxonSet &ts) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "():,");

  tokenizer tokens(s, sep);

  std::stringstream output;

  std::string prevtok = "";

  for (auto tok : tokens) {
    if (tok == "(") {
      output << "(";
    } else if (tok == ")") {
      output << ")";
    } else if (tok == ":") {
      output << ":";
    } else if (tok == ",") {
      output << ",";
    } else {
      if (prevtok == ")" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        output << tok;
        continue;
      }
      boost::algorithm::trim(tok);
      std::stringstream ss(tok);
      int i;
      ss >> i;
      const std::string &id = ts[i];
      output << id;
    }
    prevtok = tok;
  }
  output << ";";
  return output.str();
}

std::string unmap_clade_names(const std::string &s, TaxonSet &ts) {
  typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
  boost::char_separator<char> sep(";\n", "{},");

  tokenizer tokens(s, sep);

  std::stringstream output;

  std::string prevtok = "";

  for (auto tok : tokens) {
    if (tok == "{") {
      output << "{";
    } else if (tok == "}") {
      output << "}";
    } else if (tok == ",") {
      output << ",";
    } else {
      if (prevtok == "}" || prevtok == ":" || (tok == " " && prevtok == ",")) {
        continue;
      }
      boost::algorithm::trim(tok);
      std::stringstream ss(tok);
      int i;
      ss >> i;
      const std::string &id = ts[i];
      output << id;
    }
    prevtok = tok;
  }
  return output.str();
}
