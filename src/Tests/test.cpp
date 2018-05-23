#include "catch.hpp"
#include  "../TreeClade.hpp"
#include  "../newick.hpp"
#include <iostream>
#include <sstream>


TEST_CASE("Clade", "[clade]") {
  TaxonSet ts(5);
  ts.add("0");
  ts.add("1");
  ts.add("2");
  ts.add("3");
  ts.add("4");

  Clade c1(ts);
  c1.add(1);

  REQUIRE(c1.contains(1));

  Clade c2(ts);
  c2.add(2);

  REQUIRE(c2.contains(2));

  Clade c3(ts);
  c3 += c1;
  c3 += c2;


  REQUIRE(c3.contains(c1));
  REQUIRE(c3.contains(c2));

  REQUIRE(c3.contains(1));
  REQUIRE(c3.contains(2));


}


TEST_CASE("Tree", "[tree]") {
  string treestring = "((a,b),(c,(d,e)),f);";
  unordered_set<string> taxa_set;
  newick_to_ts(treestring, taxa_set);

  TaxonSet ts(taxa_set.size());
  for  (string s : taxa_set ) {
    ts.add(s);
  }

  stringstream ss;

  Tree t = newick_to_treeclades(treestring, ts);


  ss << t;

  REQUIRE(ss.str()  == treestring);
  REQUIRE(t.root().verify());

  t.reroot(ts["b"]);
  cout <<  t << endl;

  REQUIRE(((t.root().child(0).size() == 1 && t.root().child(0).contains(ts["b"])) || (t.root().child(1).size() == 1 && t.root().child(1).contains(ts["b"]))));

  DistanceMatrix lca(ts);
  t.LCA(lca);
  TreeClade& ce_lca = t.node(lca(ts["c"], ts["e"]));
  cout << ce_lca << endl;
  REQUIRE(ce_lca.contains(ts["d"]));


}
