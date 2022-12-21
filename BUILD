load("@rules_cc//cc:defs.bzl", "cc_binary", "cc_library")

cc_binary (
    name = "main",
    srcs = ["cc/main.cpp"],
    deps = [
        "//cc/switch_ml:vi_agent",
        "//cc/switch_ml:vi",
    ],
    visibility = ["//visibility:public"],
)