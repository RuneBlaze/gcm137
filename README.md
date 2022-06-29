gcm137
==============

Experimental re-make of [GCM (graph clustering merger)](https://github.com/vlasmirnov/MAGUS). The short term goal is to design a more self-contained, simply implemented GCM with reasonable (faster than original) performance. The reusable components of this project might be proposed to be reintegrated into the original GCM.

Note that as of right now, the clustering/tracing stage of this GCM implementation is not yet faster than the original implementation.

## Subcommands

### Merging

```
> gcm137 merge --help
gcm137-merge
Run GCM using existing subset and glue alignments

USAGE:
    gcm137 merge [OPTIONS] --output <OUTPUT>

OPTIONS:
    -g, --glues <GLUES>...        Glue alignments
    -h, --help                    Print help information
    -i, --input <INPUT>...        Subset alignments
    -o, --output <OUTPUT>         Output merged alignment path
    -t, --tracer <TRACER>         Tracing strategy [default: auto] [possible values: auto, upgma,
                                  pairwise]
    -w, --weights <WEIGHTS>...    Optional weights to the glues; the order corresponds to the glue alignments
```