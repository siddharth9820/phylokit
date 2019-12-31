package(default_visibility = ["//:__pkg__", "//test:__subpackages__"])

cc_binary(
    name = "libphylokit.so",
    srcs = [
        ":BitVector",
        ":Clade",
        ":DistanceMatrix",
        ":Quartet",
        ":TaxonSet",
        ":TreeClade",
        ":newick",
        "//phylokit/util:Options",
        "//phylokit/util:Timer",
    ],
    linkshared = 1,
    linkstatic = 1,
)

cc_library(
    name = "BitVector",
    srcs = ["BitVector.cpp"],
    hdrs = ["BitVector.hpp"],
)

cc_library(
    name = "Clade",
    srcs = ["Clade.cpp"],
    hdrs = ["Clade.hpp"],
    deps = [
        ":TaxonSet",
        "@boost//:algorithm",
        "@boost//:tokenizer",
    ],
)

cc_library(
    name = "TaxonSet",
    srcs = ["TaxonSet.cpp"],
    hdrs = ["TaxonSet.hpp"],
    deps = [
        ":BitVector",
        "@boost//:algorithm",
        "@boost//:tokenizer",
    ],
)

cc_library(
    name = "DistanceMatrix",
    srcs = ["DistanceMatrix.cpp"],
    hdrs = ["DistanceMatrix.hpp"],
    deps = [
        ":Clade",
        ":TaxonSet",
        "@boost//:algorithm",
        "@boost//:call_traits",
        "@boost//:multi_array",
        "@boost//:tokenizer",
	"@com_github_google_glog//:glog",
    ],
)

cc_library(
    name = "newick",
    srcs = ["newick.cpp"],
    hdrs = ["newick.hpp"],
    deps = [
        ":Clade",
        ":TaxonSet",
        ":TreeClade",
        "@boost//:algorithm",
        "@boost//:call_traits",
        "@boost//:multi_array",
        "@boost//:tokenizer",
    ],
)

cc_library(
    name = "Quartet",
    srcs = ["Quartet.cpp"],
    hdrs = ["Quartet.hpp"],
    deps = [
        ":TaxonSet",
        "//phylokit/util:Options",
        "@boost//:call_traits",
        "@boost//:multi_array",
        "@com_github_google_glog//:glog",
    ],
)

cc_library(
    name = "TreeClade",
    srcs = ["TreeClade.cpp"],
    hdrs = ["TreeClade.hpp"],
    deps = [
        ":Clade",
        ":DistanceMatrix",
    ],
)
