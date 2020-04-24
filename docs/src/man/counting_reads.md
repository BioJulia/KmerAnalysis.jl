```@meta
CurrentModule = MerCounting
```

# Counting kmers in read datasets

MerCounting provides some more dedicated counting algorithms for counting kmers
in read datasets, particularly ReadDatastores.

This is because as conceptually simple as counting kmers is, to do it quickly
for large sequencing read datasets output by sequencing machines can be
difficult for some datasets.

There are many ways you can try to optimise the kmer counting process, and many
kmer counting tools already exist.

MerCounting provides a `Counters` submodule, which contains nessecery methods
and types required to implement various kinds of kmer counter, as well as
exporting a selection of "off-the-shelf" methods that use different counting
strategies, one or several of which hopefully suit the user's dataset and
computational resources available. We describe and showcase these below.

## The `serial_mem` counter

The [`serial_mem`](@ref) counter is the simplest kmer counter for reads that
is provided.

It counts kmers in each read serially, and stores the kmers and counts entirely
in RAM.

```@docs
Counters.serial_mem
```

It is fairly simple to use. First you need a ReadDatastore, if you need to
recall how to create one then head over
to [here](https://biojulia.net/ReadDatastores.jl/stable/build-datastores/) in
ReadDatastores.jl's documentation.

The example below opens a datastore before using [`serial_mem`](@ref) to count
the kmers in the read datastore.

```@setup serialmem
using ReadDatastores
using FASTX
using BioSequences
fwq = open(FASTQ.Reader, "../../../test/ecoli_tester_R1.fastq")
rvq = open(FASTQ.Reader, "../../../test/ecoli_tester_R2.fastq")
PairedReads{DNAAlphabet{2}}(fwq, rvq, "ecoli-test-paired", "my-ecoli-test", 250, 300, 0, FwRv)
```

```@repl serialmem
using MerCounting, ReadDatastores
ds = @openreads "ecoli-test-paired.prseq"
kl = Counters.serial_mem(DNAMer{31}, ds, CANONICAL)
```

## The `dist_mem` counter

The [`dist_mem`](@ref) counter is a multi-process version of [`serial_mem`](@ref).

```@docs
Counters.dist_mem
```

It is used in a very similar fashion to [`serial_mem`](@ref), but you need to
have multiple julia worker processes available to gain any benefit.

```@repl distmem
# Spin up some extra julia worker processes, make sure they are using the
# appropriate Project.toml file for your project.
using Distributed
addprocs(8, exeflags="--project=../../../docs/")

# Make sure the packages are loaded on all the workers.
@everywhere using MerCounting, ReadDatastores, BioSequences

ds = @openreads "ecoli-test-paired.prseq"
kl = Counters.dist_mem(DNAMer{31}, ds, CANONICAL)
```

Naturally, for such a small and simple example as this, using [`dist_mem`](@ref)
for such a small dataset is probably just worse than doing it in memory serially,
because of the overheads. For real datasets you should see a benefit.