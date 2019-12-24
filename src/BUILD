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
        "//src/util:Logger",
        "//src/util:Options",
        "//src/util:Timer",
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
        "//src/util:Logger",
    ],
)

cc_library(
    name = "TaxonSet",
    srcs = ["TaxonSet.cpp"],
    hdrs = ["TaxonSet.hpp"],
    deps = [":BitVector"],
)

cc_library(
    name = "DistanceMatrix",
    srcs = ["DistanceMatrix.cpp"],
    hdrs = ["DistanceMatrix.hpp"],
    deps = [
        ":Clade",
        ":TaxonSet",
        "//src/util:Logger",
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
    ],
)

cc_library(
    name = "Quartet",
    srcs = ["Quartet.cpp"],
    hdrs = ["Quartet.hpp"],
    deps = [
        ":TaxonSet",
        "//src/util:Logger",
        "//src/util:Options",
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
