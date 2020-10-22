#include <string>
#include "catch2.hpp"
#include <iostream>
#include "phylokit/Clade.hpp"
#include "phylokit/newick.hpp"

using std::endl;
using std::cerr;

TEST_CASE("rooted tree is_rooted") {
  REQUIRE(is_rooted("((a, b), c)"));
  REQUIRE(is_rooted("(c, (a, b))"));
  REQUIRE(is_rooted("(((a, x), (b, c)), d)"));
  REQUIRE(is_rooted("(d, ((a, x), (b, c)))"));
}

TEST_CASE("unrooted tree is_rooted") {
  REQUIRE(!is_rooted("(a, b, c)"));
  REQUIRE(!is_rooted("((a, x), (b, c), d)"));
  REQUIRE(!is_rooted("((a, x), d, (b, c))"));
  REQUIRE(!is_rooted("(d, (a, x), (b, c))"));
}

TEST_CASE("deroot unrooted tree") {
  REQUIRE(deroot("(a, b, c)") == "(a, b, c)");
  REQUIRE(deroot("((a, x), (b, c), d)") == "((a, x), (b, c), d)");
}

TEST_CASE("deroot rooted tree") {
  REQUIRE(deroot("((a, b), c)") == "(a, b, c)");
  REQUIRE(deroot("(c, (a, b))") == "(c, a, b)");
  REQUIRE(deroot("(((a, x), (b, c)), d)") == "((a, x), (b, c), d)");
  REQUIRE(deroot("(d, ((a, x), (b, c)))") == "(d, (a, x), (b, c))");
  REQUIRE(deroot("((a, x), (b, c))") == "(a, x, (b, c))");
}

TEST_CASE("newick_to_ts") {
  std::unordered_set<std::string> output_set;
  SECTION("Simple names") {
    std::unordered_set<std::string> expected{"a", "b", "c", "d", "e"};
    REQUIRE(newick_to_ts("((a, b), (c, (d, e)))", output_set) == 5);
    REQUIRE(output_set == expected);
  }
  SECTION("names with special characters") {
    std::unordered_set<std::string> expected{"first taxon", "second taxon",
                                             "third taxon", "\"fourth taxon\""};
    REQUIRE(
        newick_to_ts(
            "(first taxon, (second taxon, (third taxon, \"fourth taxon\")))",
            output_set) == 4);
    REQUIRE(output_set == expected);
  }
  SECTION("branch lengths") {
    std::unordered_set<std::string> expected{"a", "b", "c", "d"};
    REQUIRE(newick_to_ts("(a:4.3, (b, (c, d):12):3.23)", output_set) == 4);
    REQUIRE(output_set == expected);
  }
  SECTION("branch labels") {
    std::unordered_set<std::string> expected{"a", "b", "c", "d"};
    REQUIRE(newick_to_ts("(a, (b, (c, d)I1)I2)I3", output_set) == 4);
    REQUIRE(output_set == expected);
  }
}

TEST_CASE("newick_to_taxa") {
  TaxonSet ts(
      "a,b,c,d,e,f,g,first taxon,second taxon,third taxon,\"fourth taxon\"");
  SECTION("Simple names") {
    REQUIRE(newick_to_taxa("((a, b), (c, (d, e)))", ts) ==
            Clade(ts, "a,b,c,d,e"));
  }
  SECTION("names with special characters") {
    REQUIRE(
        newick_to_taxa(
            "(first taxon, (second taxon, (third taxon, \"fourth taxon\")))",
            ts) ==
        Clade(ts, "first taxon,second taxon,third taxon,\"fourth taxon\""));
  }
  SECTION("branch lengths") {
    std::unordered_set<std::string> expected{"a", "b", "c", "d"};
    REQUIRE(newick_to_taxa("(a:4.3, (b, (c, d):12):3.23)", ts) ==
            Clade(ts, "a,b,c,d"));
  }
  SECTION("branch labels") {
    std::unordered_set<std::string> expected{"a", "b", "c", "d"};
    REQUIRE(newick_to_taxa("(a, (b, (c, d)I1)I2)I3", ts) ==
            Clade(ts, "a,b,c,d"));
  }
}

TEST_CASE("newick_to_clades") {
  TaxonSet ts("a,b,c,d,e,f");
  std::unordered_set<Clade> clade_set;
  SECTION("Single tree") {
    newick_to_clades("(a, ((b, c), (d, e)))", ts, clade_set);
    REQUIRE(clade_set == std::unordered_set<Clade>{
                             Clade(ts, "b,c"), Clade(ts, "d,e"),
                             Clade(ts, "b,c,d,e"), Clade(ts, "a,b,c,d,e")});
  }
  SECTION("Polytomy") {
    newick_to_clades("(a, (b, c), (d,e))", ts, clade_set);
    REQUIRE(clade_set == std::unordered_set<Clade>{Clade(ts, "b,c"),
                                                   Clade(ts, "d,e"),
                                                   Clade(ts, "a,b,c,d,e")});
  }
}

TEST_CASE("newick_to_treeclades") {
  TaxonSet ts("a,b,c,d,e,f");

  SECTION("Single tree") {
    Tree tree = newick_to_treeclades("(a, ((b, c), (d, e)))", ts);
    
    CHECK(tree.root() == Clade(ts, "a,b,c,d,e"));
    CHECK(tree.root().nchildren() == 2);
    CHECK(tree.root().child(0) == Clade(ts, "a"));
    CHECK(tree.root().child(1) == Clade(ts, "b,c,d,e"));
    CHECK(tree.root().child(1).child(0) == Clade(ts, "b,c"));
    CHECK(tree.root().child(1).child(1) == Clade(ts, "d,e"));
    CHECK(tree.clades.size() == 9);
    
    for (auto& c_it : tree.clades) {
      for (int i = 0; i < c_it.second.nchildren(); i++) {
        CHECK(c_it.second.child(i).parent == c_it.first);
      }
    }
  }
  SECTION("Polytomy") {
    Tree tree = newick_to_treeclades("(a, (b, c), (d,e))", ts);

    CHECK(tree.root() == Clade(ts, "a,b,c,d,e"));
    CHECK(tree.root().nchildren() == 3);
    CHECK(tree.root().child(0) == Clade(ts, "a"));
    CHECK(tree.root().child(1) == Clade(ts, "b,c"));
    CHECK(tree.root().child(2) == Clade(ts, "d,e"));

    for (auto& c_it : tree.clades) {
      for (int i = 0; i < c_it.second.nchildren(); i++) {
        CHECK(c_it.second.child(i).parent == c_it.first);
      }
    }
  }
}

TEST_CASE("newick_to_postorder") {
  TaxonSet ts("a,b,c,d,e,f");
  std::vector<Taxon> order;
  SECTION("Simple tree") {
    newick_to_postorder("(a, ((b, c), (d, e)))", ts, order);
    REQUIRE(order == std::vector<Taxon>{ts["a"], ts["b"], ts["c"], -2, ts["d"],
                                        ts["e"], -2, -2, -2});
  }
  SECTION("Polytomy") {
    newick_to_postorder("(a, (b, c), (d, e))", ts, order);
    REQUIRE(order == std::vector<Taxon>{ts["a"], ts["b"], ts["c"], -2, ts["d"],
                                        ts["e"], -2, -3});
  }
}

TEST_CASE("map_newick_names") {
  TaxonSet ts(5);
  ts.add("a");
  ts.add("b");
  ts.add("c");
  ts.add("d");
  ts.add("e");
  SECTION("Simple tree") {
    REQUIRE(map_newick_names("(a, (b, c), (d, e))", ts) == "(0, (1,2), (3,4));");
  }
  SECTION("Branch lengths") {
    REQUIRE(map_newick_names("(a, (b, c):1.2, (d, e))", ts) ==
            "(0, (1,2):1.2, (3,4));");
  }
}

TEST_CASE("unmap_newick_names") {
  TaxonSet ts(5);
  ts.add("a");
  ts.add("b");
  ts.add("c");
  ts.add("d");
  ts.add("e");
  SECTION("Simple tree") {
    REQUIRE(unmap_newick_names("(0,(1,2),(3,4));", ts) ==
            "(a,(b,c),(d,e));");
  }
  SECTION("Branch lengths") {
    REQUIRE(unmap_newick_names("(0,(1,2):1.2, (3,4));", ts) ==
            "(a,(b,c):1.2, (d,e));");
  }
}

TEST_CASE("unmap_clade_names") {
  TaxonSet ts(5);
  ts.add("a");
  ts.add("b");
  ts.add("c");
  ts.add("d");
  ts.add("e");
  SECTION("Simple case") {
    REQUIRE(unmap_clade_names("{1,2,3}", ts) == "{b,c,d}");
  }
}