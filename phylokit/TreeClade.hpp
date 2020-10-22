#ifndef __TREECLADE_HPP__
#define __TREECLADE_HPP__

#include <iostream>
#include "Clade.hpp"
#include "DistanceMatrix.hpp"
class Tree;

class TreeClade : public Clade {
 private:
 public:
  std::vector<int> children_;
  int parent;
  int index;
  Tree &tree;
  using Clade::Clade;
  TreeClade(TaxonSet &ts, Tree &tree, int index)
      : Clade(ts), index(index), tree(tree) {}
  TreeClade(TaxonSet& ts, Tree &tree, const TreeClade& other)
      : Clade(other), 
        children_(other.children_), 
        parent(other.parent),
        index(other.index), 
        tree(tree) {}
  TreeClade(TaxonSet& ts, Tree &tree, const TreeClade&& other)
      : Clade(other), 
        children_(other.children_), 
        parent(other.parent),
        index(other.index), 
        tree(tree) {}
  void addChild(int index);
  std::vector<int> &children();
  const std::vector<int> &children() const;
  int nchildren() const { return children_.size(); }
  TreeClade &child(int i);
  const TreeClade &child(int i) const;
  Clade complement() const;
  bool verify();
};

class Tree {
 private:
 public:
  std::unordered_map<int, TreeClade> clades;
  int next_entry;
  TaxonSet &ts;

  Tree(TaxonSet &ts) : next_entry(0), ts(ts) {}
  Tree(const Tree &other) : 
    next_entry(other.next_entry),
    ts(other.ts) {
      for (auto& i : other.clades) {
        clades.emplace(i.first, TreeClade(ts, *this, i.second));
      }
  }
  Tree(Tree &&other)  : 
    next_entry(other.next_entry),
    ts(other.ts) {
      for (auto& i : other.clades) {
        clades.emplace(i.first, TreeClade(ts, *this, i.second));
      }
  }

  TreeClade &root() { return clades.at(0); }
  const TreeClade &root() const { return clades.at(0); }
  TreeClade &node(int n) { return clades.at(n); }
  const TreeClade &node(int n) const { return clades.at(n); }
  int addNode() {
    Tree &me = *this;
    // clades.emplace(piecewise_construct, make_tuple(next_entry),
    // make_tuple(ts, me, next_entry));
    clades.insert(std::make_pair(next_entry, TreeClade(ts, me, next_entry)));
    clades.at(next_entry).parent = -1;
    next_entry++;
    return next_entry - 1;
  }

  const Clade &taxa() const { return root(); }

  Tree &binary_root(int a);

  void swap(int a, int b);

  Tree &rotate(int a_i, int b_i);

  Tree &reroot(Taxon x);

  void LCA(DistanceMatrix &lca) const;

  double RFDist(const Tree &other, bool normalized = true) const;
};
std::ostream &operator<<(std::ostream &os, const Tree &t);
std::ostream &operator<<(std::ostream &os, const TreeClade &t);

#endif
