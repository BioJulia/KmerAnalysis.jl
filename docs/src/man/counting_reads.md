```@meta
CurrentModule = KmerAnalysis
```

# Counting k-mers in read datasets

MerCounting provides some more dedicated counting algorithms for counting k-mers
in read datasets, particularly ReadDatastores.

This is because as conceptually simple as counting k-mers is, to do it quickly
for large sequencing read datasets output by sequencing machines can be
difficult for some datasets.

There are many ways you can try to optimise the k-mer counting process, and many
k-mer counting tools already exist.

MerCounting provides a `Counters` submodule, which contains nessecery methods
and types required to implement various kinds of k-mer counter, as well as
exporting a selection of "off-the-shelf" methods that use different counting
strategies, one or several of which hopefully suit the user's dataset and
computational resources available. We describe and showcase these below.

## The `SerialMemCounter` counter

The [`SerialMemCounter`](@ref) counter is the simplest k-mer counter for reads that
is provided.

It counts k-mers in each read serially, and stores the k-mers and counts entirely
in RAM.

```@docs
SerialMemCounter
```

It is fairly simple to use. First you need a ReadDatastore, if you need to
recall how to create one then head over
to [here](https://biojulia.dev/ReadDatastores.jl/stable/build-datastores/) in
ReadDatastores.jl's documentation.

The example below opens a datastore before creating a [`SerialMemCounter`](@ref)
with [`serial_mem`](@ref), to count the k-mers in the read datastore.

```@setup serialmem
using ReadDatastores
using FASTX
using BioSequences
fwq = open(FASTQ.Reader, "../../../test/ecoli_tester_R1.fastq")
rvq = open(FASTQ.Reader, "../../../test/ecoli_tester_R2.fastq")
PairedReads{DNAAlphabet{2}}(fwq, rvq, "ecoli-test-paired", "my-ecoli-test", 250, 300, 0, FwRv)
```

```@repl serialmem
using KmerAnalysis, ReadDatastores
ds = @openreads "ecoli-test-paired.prseq"
c = serial_mem(DNAMer{31}, CANONICAL)
kl = c(ds)
```

## The `dist_mem` counter

The [`DistMemCounter`](@ref) counter is a multi-process version of [`SerialMemCounter`](@ref).

```@docs
DistMemCounter
```

It is used in a very similar fashion to [`serial_mem`](@ref), but you need to
have multiple julia worker processes available to gain any benefit.

```@repl distmem
# Spin up some extra julia worker processes, make sure they are using the
# appropriate Project.toml file for your project.
using Distributed
addprocs(8, exeflags="--project=../../../docs/")

# Make sure the packages are loaded on all the workers.
@everywhere using KmerAnalysis, ReadDatastores, BioSequences

ds = @openreads "ecoli-test-paired.prseq"
dmc = dist_mem(DNAMer{31}, CANONICAL)
kl = dmc(ds)
```

Naturally, for such a small and simple example as this, using [`DistMemCounter`](@ref)
for such a small dataset is probably just worse than doing it in memory serially,
because of the overheads. For real datasets you should see a benefit.