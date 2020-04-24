# Analysing kmer count data

Once you have kmer count data, there are a number of useful analyses that may
be done. Let's go through a few together.

## Kmer frequency spectra

A kmer frequency spectra is a (often) graphical representation of a dataset
showing how many kmers in a kmer count dataset appear a certain number of times.
In other words, how many kmers in the dataset have a count of 1, how many have a
count of 10, and so on.

When plotted, the frequency of occurence is plotted on the x-axis, and the
number of kmers on the y-axis.

The kmer spectra is composed of distributions representing groups of motifs at
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

Here I'm going to load in these reads, count the kmers with the [`serial_mem`](@ref)
kmer counter as the dataset is not that large, and then compute the frequency
spectra.

```@repl
using MerCounting, ReadDatastores, BioSequences

reads = @openreads "../../../test/fakemicrobe.prseq"
kmer_counts = Counters.serial_mem(DNAMer{31}, reads, CANONICAL) 

my_spectra = spectra(kmer_counts)
```