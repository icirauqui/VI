workspace(name = "main")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")


# Setup rules_foreign_cc (CMake integration)
http_archive(
    name = "rules_foreign_cc",
    strip_prefix = "rules_foreign_cc-8d540605805fb69e24c6bf5dc885b0403d74746a", # 0.9.0
    url = "https://github.com/bazelbuild/rules_foreign_cc/archive/8d540605805fb69e24c6bf5dc885b0403d74746a.tar.gz",
)

load("@rules_foreign_cc//foreign_cc:repositories.bzl", "rules_foreign_cc_dependencies")

rules_foreign_cc_dependencies()





_ALL_CONTENT = """\
filegroup(
    name = "all_srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)
"""



http_archive(
      name = "sciplot",
      build_file = "//:third_party/sciplot.BUILD",
      sha256 = "d5f8fa7783920c0db4f54d8d9cea8a5e50be99153e81dcaccf91532994ebbad0",
      strip_prefix = "sciplot-0.3.1",
      urls = ["https://github.com/sciplot/sciplot/archive/refs/tags/v0.3.1.tar.gz"],
)




http_archive(
      name = "gsl",
      build_file = "//:third_party/gsl.BUILD",
      sha256 = "d719d9af948e6d03424d35f0b272a30e54448e958491ad0415ba4ee69b781ebf",
      strip_prefix = "gsl-20211111",
      urls = ["https://github.com/ampl/gsl/archive/refs/tags/20211111.tar.gz"],
)







#http_archive(
#    name = "gsl",
#    url = "https://mirror.freedif.org/GNU/gsl/gsl-latest.tar.gz",
#    build_file = '@//third_party:gsl.BUILD',
#)
