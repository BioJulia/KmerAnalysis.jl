# Analysing k-mer count data

Once you have k-mer count data, there are a number of useful analyses that may
be done. Let's go through a few together.

## K-mer frequency spectra

A k-mer frequency spectra is a representation of sequencing data showing how many
k-mers appear a certain number of times.
In other words, how many k-mers in the dataset have a count of 1, how many have a
count of 10, and so on.

When plotted, the frequency of occurence is plotted on the x-axis, and the
number of k-mers on the y-axis.

The k-mer spectra is composed of distributions representing groups of motifs at
different frequencies in the sample, plus biases.
Given not too many biases, the shape of the distributions provides a useful set
of properties describing the biological sample, the sequencing process and the
amount of useful data in the dataset.

Let's take a look at a simple spectra. I'm going to use a small artificial
dataset of a bacterial-like genome I generated using
[Pseudoseq.jl](https://bioinfologics.github.io/Pseudoseq.jl/stable/).

I first created a random single chromosome genome, and then simulated a
paired-end sequencing experiment to generate the kind of raw data we would get
from a sequencing experiment.

!!! note
    I'm doing this with a small artificial genome to be kind to the CI
    infrastructure that builds these docs!!!
    
    But you could follow along and try yourself with paired-end read data from
    a sequencing experiment from e.coli.

Here I'm going to load in these reads, count the k-mers with the [`SerialMemCounter`](@ref)
k-mer counter as the dataset is not that large, and then compute the 1D frequency
spectra.

```@repl onedspectra
using KmerAnalysis, ReadDatastores, BioSequences

reads = @openreads "../../../test/fakemicrobe.prseq"
counter = serial_mem(DNAMer{31}, CANONICAL)
kmer_counts = counter(reads) 

my_spectra = spectra(kmer_counts)
```

I could have also created the counter and passed it to the [`spectra`](@ref)
method along with the reads like so:

```@repl onedspectra
my_spectra = spectra(counter, reads)
```

At this point you can plot the k-mer frequency spectra using Makie.jl, as
KmerAnalysis.jl includes some AbstractPlotting/Makie recipe's, so you can
visualise the k-mer frequency spectra quite easily, with several of Makie.jl's
plotting primitives.

```julia
using Makie

plot(my_spectra)
```

!!! note
    Loading and using Makie without any pre-compilation of packages can be a bit
    slow as the JIT has to compile quite a lot.

## Kmer count projection

KmerAnalysis includes a k-mer count container type called IndexedCounts.
Indexed counts are different to a vector or other container of `MerCount`, in
that IndexedCounts store k-mer counts from many sources (often different read
datasets) indexed against a reference dataset (usually a genome graph or set of
contigs).