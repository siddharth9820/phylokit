#include "catch2.hpp"
#include "phylokit/newick.hpp"

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