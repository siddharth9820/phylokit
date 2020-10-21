/*#ifndef __TREE_HPP_
#define __TREE_HPP_

typedef uid_t int;
typedef explicit edge_uid uid_t;
typedef explicit node_uid uid_t;

class UIDResolver {
  int e_uid;
  int n_uid;
  vector<Edge> edges;
  vector<Node> nodes;
  Node& operator(node_uid uid) {
    return nodes.at(uid);
  }
  Edge& operator(edge_uid uid) {
    return edges.at(uid);
  }
  Node& newNode() {
    nodes.push_back(Node());
    nodes.back().uid = nodes.size() - 1;
    return nodes.back();
  }
  Edge& newEdge(node_uid a, node_uid b) {
    edges.push_back(Edge());
    edges.back().a = a;
    edges.back().b = b;
    edges.back().uid = edges.size() - 1;
    return edges.back();
  }
  Edge& newEdge(Node& a, Node& b) {
    return newEdge(a.uid, b.uid);
  }
};

struct Edge;
struct Node;

struct Edge {
  UIDResolver& r_;
  edge_uid uid;
  node_uid a;
  node_uid b;

  Bipartition split;
  Node& other(Node& x) {
    if (a == x.uid) {
      return r_(b);
    }
    assert(b==x.uid);
    return r_(a);
  }

  node_uid other(node_uid x) {
    if (a == x) {
      return b;
    }
    assert(b==x);
    return a;
  }
};


struct Node {
  UIDResolver& r_;

  node_uid uid;
  Taxon taxon;

  vector<edge_uid> edges;

  int root_dist;

  bool isLeaf() {
    return neighbors.size() == 1;
  }

  Node& addChild() {
    Node& newNode = r_.newNode();
    Edge& newEdge = r_.newEdge(*this, newNode);
    newNode.edges.push_back(newEdge.uid);
    edges.push_back(newEdge.uid);
    newNode.root_dist = root_dist + 1;
  }

  Edge& operator[](Node& other) {
    for (edge_uid u : edges) {
      Edge& e = r_(u);
      if(e.other(uid) == other.uid) {
        return e;
      }
    }
    assert(0);
  }

  struct ChildIterator {
    Node& n;
    int current
    ChildIterator(Node& n) :
    n(n),
    current(0)
    { }

    BVFIterator& operator++() {
      increment();
      return *this;
    }
    void increment();
    int operator*() {
      if n.edges[current]
    }
    bool operator!=(const BVFIterator& other) {
      return (current != other.current);
    }
    bool operator==(const BVFIterator& other) {
      return (current == other.current);
    }
  };
  ChildIterator children() {
    return ChildIterator(*this);
  }



};

template <typename EdgeAttrs, NodeAttrs>
class Tree {
  TaxonSet& ts;
  uid_t root;
  bool splitsGood;
  Clade taxa;
  unordered_map<uid_t, Edge<EdgeAttrs> > edges_;
  unordered_map<uid_t, Node<NodeAttrs> > nodes_;
  UIDResolver r_;
public:
  Tree(TaxonSet& ts) :
  ts(ts),
  root(0),
  splitsGood(0),
  taxa(ts) {}

  Node& operator[](node_uid uid) {
    return r_(uid);
  }



  Node& newNode(Node* other) {
    uid_t uid = gen_uid();
    nodes_[uid] = Node(uid);
    nodes_[uid].neighbors.push_back(active.back());
    active.back().neighbors.push_back(nodes_[uid]);
  }




  void updateSplits() {
    if (splitsGood)
      return;
    splitsGood = true;

    unordered_map<Node*, int> unprocessed_adj;
    unordered_set<uid_t> unprocessed_edges;
    vector<Node*> to_process;

    for (Node& n : nodes_) {
      if (n.adjacent.size() == 1) {
        to_process.push_back(&n);
      } else {
        unprocessed_adj[&n] = n.adjacent.size();
      }

    }

    for (Edge& e : edges) {
      unprocessed_edges.insert(e.uid);
    }


    while (to_process.size()) {
      Node* n = to_process.back();
      to_process.pop_back();
      Clade c1(ts);
      Edge* target;
      for (Edge& e : n->adjacent) {
        if (unprocessed_edges.count(e.uid)) {
          target = &e;
        } else {
          c1 += e.split.a1;
        }
      }

      target->split = Bipartition(c1, taxa - c1);
      unprocessed_edges.remove(target.uid);

      unprocessed_adj[target->other(n)->uid]--;
      if (unprocessed_adj[target->other(n)->uid] == 1) {
        to_process.push_back(target->other(n));
      }
    }
  }

  void rootAt(Edge<EdgeAttrs>& e) {
    root = e.uid;
    vector<pair<Node*, int> > stack;
    e.a->root_dist = 1;
    e.b->root_dist = 1;
    for (auto n : e.a.neighbors) {
      stack.push_back(make_pair(n, 2)));
    }
    for (auto n : e.b.neighbors) {
      stack.push_back(make_pair(n, 2)));
    }

    while(stack.size()) {
      Node* n;
      int dist;
      tie(n, dist) = stack.back();
      stack.pop_back();
      n->root_dist = dist;
      for (auto n1 : n->neighbors) {
        stack.push_back(make_pair(n1, dist + 1));
      }
    }

  }

  stringstream newick(Node* n, Node* from) {
    stringstream ss;
    int count = 0;
    ss << "(";
    for (Node* n1 : n->neighbors) {
      if (n1 == from)
        continue;
      if (count != 0) {
        ss << ",";
      }
      count ++;
      if(n1->isLeaf()) {
        ss << n1.taxon;
      } else {
        ss << "(" << newick(n1, n) << ")";
      }
    }
    return ss;

  }

  string newick() {
    for (Node* n: nodes) {
      if (!n->isLeaf()) {
        return newick(n, 0).str();
      }
    }
  }


  ostream& operator<<(ostream& os, const Clade& c) {
    os << c.str();
    return os;
  }

};

#endif*/
