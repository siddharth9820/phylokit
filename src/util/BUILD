package(default_visibility = ["//visibility:public"])

cc_library(
    name = "Options",
    srcs = ["Options.cpp"],
    hdrs = ["Options.hpp"],
)

cc_library(
    name = "Logger",
    srcs = ["Logger.cpp"],
    hdrs = ["Logger.hpp"],
    deps = [":Options"],
)

cc_library(
    name = "Timer",
    srcs = ["Timer.cpp"],
    hdrs = ["Timer.hpp"],
)
