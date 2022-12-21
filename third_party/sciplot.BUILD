load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

package(
    default_visibility = ["//visibility:public"],
)

filegroup(
    name = "all_srcs",
    srcs = glob(["**"]),
)

# arrayfire
cmake(
    name = "sciplot",
    cache_entries = {
        "CMAKE_BUILD_TYPE": "Release",
    },
    build_args = [
        "--verbose",
        "--",
        "-j `nproc`",
    ],
    copts = [
        "-fPIC",
    ],
    lib_source = "@sciplot//:all_srcs",
    #out_shared_libs = [
    #  "libaf.so", 
    #  "libaf.so.3", 
    #  "libafcpu.so",
    #  "libafcpu.so.3",
    #],
    visibility = ["//visibility:public"],
)