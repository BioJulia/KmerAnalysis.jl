```@meta
CurrentModule = MerCounting
```

# Basic types and basic kmer counting

MerCounting has a few basic types and methods to allow you do easily do some
basic kmer counting.

The first and perhaps most obvious of these is a type to represent a kmer and
its frequency, for this, MerCounting provides the [`MerCount`](@ref) type.

If you have a [`MerCount`](@ref) variable you can get the kmer value or the
frequency value using the [`mer`](@ref) and [`freq`](@ref) getter methods.

A vector of these [`MerCount`](@ref)s constitute one option for a simple data
structure for representing kmer frequency data. This packages defines other more
dedicated types, but we will get to those in later sections of this manual.

Ok, so let's do some very basic kmer counting for a sequence!

First we need a sequence:

```@repl acount
using MerCounting
using BioSequences
s = randdnaseq(50)
```

Ok, let's say we wanted to count the 7-mers, we can collect all the 7-mers using
BioSequences' kmer iterator:

```@repl acount
sevenmers = collect(each(DNAMer{7}, s))
sevenmers[1]
```

Ok BioSequences' kmer iterator on its own yields a `MerIterResult`, which
contains a tuple of (Base position in sequence, kmer, reverse complement of the
kmer).

Now let's say I was only interested in counting canonical kmers in `s`. For any
kmer and its reverse complement, the canonical version is the one that
is lexicographically less.

So let's revise our `collect` above to give us a vector of the canonical kmers.
Thankfully, BioSequences comes with a method called `canonical`, for exactly
this purpose.

```@repl acount
sevenmers = collect(canonical(x) for x in each(DNAMer{7}, s))
```

If you wanted to collect all the kmers in `s` "as is" as it were, you can use
the `fwmer` method provided with BioSequences for this purpose.

```@repl acount
collect(fwmer(x) for x in each(DNAMer{7}, s))
```

This is quite nice and terse code, but we haven't used MerCounting to help us at
all yet. MerCounting allows us to achieve the above with the use of a single
function called [`collect_mers`](@ref), instead of `collect`ing a generator:

```@repl acount
sevenmers = collect_mers(CANONICAL, each(DNAMer{7}, s))
```

There are a few different methods of [`collect_mers`](@ref) for different types
of input, check out the API reference for more info.

MerCounting exports `CANONICAL` and `NONCANONICAL`, which work with [`collect_mers`](@ref)
and other MerCounting functions to dictate how the kmers are produced.

Ok so we have our canonical kmers, how can we count them? Well one of the ways
to count the number of distinct values in a vector is to sort it and traverse it
in one go. The [`collapse_into_counts!`](@ref) methods do this for you.

!!! note
    [`collapse_into_counts!`](@ref) will sort the input `sevenmers` vector in
    place.

```@repl acount
mercounts = collapse_into_counts!(sevenmers)
```

Looks like every distinct kmer appeared once, with the possible odd exceptions.
Well that's to be expected! We did just do this with a single 50bp long random
sequence after all!

So for any sequence or input we want to count kmers for, we have to [`collect_mers`](@ref)
and then collapse them into counts with [`collapse_into_counts!`](@ref).

For convenience and more terse code, MerCounting provides a constructor for
vectors of [`MerCount`](@ref) that does this for you.

```@repl acount
v = Vector{MerCount{DNAMer{7}}}(CANONICAL, s)
```

That's about all there is for doing basic kmer counting over individual
biological sequences in julia.

Furthur on in this manual we will cover some dedicated kmer counting algorithms
MerCounting provides for larger datasets like those contained in reads
(using [ReadDatastores.jl](https://biojulia.net/ReadDatastores.jl/stable)).
We will also cover dedicated types such as kmer frequency spectra, and types
useful for projecting kmer coverage and other kmer related metrics over a genome
sequence.