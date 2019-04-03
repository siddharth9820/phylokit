#include "TreeClade.hpp"


Clade TreeClade::complement() const {
    return tree.root() - *this;
}

TreeClade& TreeClade::child(int i) {
    return tree.node(children_.at(i));
}
const TreeClade& TreeClade::child(int i) const {
    return tree.node(children_.at(i));
}
bool TreeClade::verify() {
  Clade child_taxa(ts);

    for (int i : children_) {
      if (tree.node(i).parent != index) {
        cout << "Node " << i << " : " << tree.node(i) << " has wrong parent" << endl;
        return false;
      }
      if (!tree.node(i).verify()) {
        return false;
      }

      child_taxa += tree.node(i);
    }

    if (size() > 1 && ! (child_taxa == *this)) {
      cout << "Node " << index << " : " << tree.node(index) << " has taxa " << static_cast<Clade>(*this) << endl;
      cout << "when it should have taxa " << child_taxa << endl;
      cout << "Children : " << endl;
      for (int i : children_) {
        cout << tree.node(i) << " :: ";
        cout << static_cast<Clade>(tree.node(i)) << endl;
      }
      return false;
    }

    return true;
  }

void TreeClade::addChild(int i) {
  children_.push_back(i);
  tree.node(i).parent = index;
}

vector<int>& TreeClade::children() {
  return children_;
}


const vector<int>& TreeClade::children() const {
  return children_;
}

std::ostream& operator<<(std::ostream& os, const TreeClade& tc) {
  if (tc.size() == 1) {
    for (Taxon t : tc.taxa)
      os  << tc.ts[t];
    return os;
  }

  os << "(";
  int first = 1;
  for (int i = 0; i < tc.nchildren(); i++) {
    if (!first)
      os << ",";
    first = 0;
    os << tc.child(i);
  }
  os << ")";
  return os;
}

  ostream& operator<<(ostream& os, const Tree& t){
    os << t.root() << ";";
    return os;
  }


Tree& Tree::binary_root(int a) {
  if (root().nchildren() <= 2) {
    return *this;
  }

  int c1 = root().children().at(a);
  int c2 = addNode();
  node(c2).parent = 0;
  for (int i = 0; i < root().nchildren(); i++) {
    if (i == a) {
      continue;
    }
    node(c2) += root().child(i);
    node(c2).addChild(root().child(i).index);
    root().child(i).parent = c2;
  }
  root().children().clear();
  root().addChild(c1);
  root().addChild(c2);
  return *this;
}



void Tree::swap(int a, int b) {
  int pa = node(a).parent;
  int pb = node(b).parent;


  for (int i = 0; i < node(pa).nchildren(); i++) {
    if (node(pa).children()[i] == a) {
      node(pa).children()[i] = b;
    }
  }


  for (int i = 0; i < node(pb).nchildren(); i++) {
    if (node(pb).children()[i] == b) {
      node(pb).children()[i] = a;
    }
  }

  node(a).parent = pb;
  node(b).parent = pa;


}

Tree& Tree::rotate(int a_i, int b_i) {
    TreeClade& acomp = root().child(1 - a_i);
    TreeClade& a = root().child(a_i);
    TreeClade& b = a.child(b_i);

    swap(acomp.index, b.index);

    a -= b;
    a += acomp;

    return *this;
  }


Tree& Tree::reroot(Taxon x) {
  if (root().nchildren() > 2) {
    binary_root(0);
  }
  binary_root(0);

  while (true) {
    for (int i= 0; i < root().nchildren(); i++) {
        if (root().child(i).contains(x)) {
          if (root().child(i).size() == 1) {
            return *this;
          }
          for (int j = 0; j < root().child(i).nchildren(); j++) {
            if (root().child(i).child(j).contains(x)) {
              rotate(i, j);
              break;
            }
          }
          break;
        }
    }
  }

  return *this;
}


void Tree::LCA(DistanceMatrix& lca) const {


  vector<int> stack;
  stack.push_back(0);

  while (stack.size()) {
    int ind = stack.back();
    const TreeClade& tc = node(ind);
    stack.pop_back();

    for (Taxon i : tc) {
      for (Taxon j : tc) {
        lca(i,j) = ind;
      }
    }

    for (int i = 0; i < tc.nchildren(); i++) {
      stack.push_back(tc.children().at(i));
    }

  }

}

double Tree::RFDist(const Tree& other, bool normalized) const {
  unordered_set<Clade> my_clades;
  for (size_t i = 1; i < clades.size(); i++) {
    Clade ol = clades.at(i).overlap(other.taxa());
    my_clades.emplace(ol);
//    cout << "adding " << ol << endl;
  }

  double matching = 0;
  double count    = 0;

  for (size_t i = 1; i < other.clades.size(); i++) {
    if (other.clades.at(i).overlap_size(taxa()) <= 1) {
      continue;
    }
    count ++;
    if (my_clades.count(other.clades.at(i).overlap(taxa()))) {
      matching ++;
    } else {
  //    cout << "couldn't find " << other.clades.at(i) << endl;
    }
  }
  if (normalized)
    return 1-(matching/count);
  else
    return count - matching;
}
