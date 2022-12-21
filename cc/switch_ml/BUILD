load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")


######################
## Switch-ML builds ##
######################

cc_library(
    name = "vi_agent",
    hdrs = ["vi/vi_agent.h"],
    srcs = ["vi/vi_agent.cpp"],
    linkstatic=True,
    visibility = ["//visibility:public"],
)

cc_library(
    name = "vi",
    hdrs = ["vi/vi.h"],
    srcs = ["vi/vi.cpp"],
    linkstatic=True,
    deps = [
        "@gsl//:gsl",
    ],
    visibility = ["//visibility:public"],
)