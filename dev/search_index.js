var documenterSearchIndex = {"docs":
[{"location":"man/analyses/#Analysing-k-mer-count-data-1","page":"Analysing Kmer count data","title":"Analysing k-mer count data","text":"","category":"section"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"Once you have k-mer count data, there are a number of useful analyses that may be done. Let's go through a few together.","category":"page"},{"location":"man/analyses/#K-mer-frequency-spectra-1","page":"Analysing Kmer count data","title":"K-mer frequency spectra","text":"","category":"section"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"A k-mer frequency spectra is a representation of sequencing data showing how many k-mers appear a certain number of times. In other words, how many k-mers in the dataset have a count of 1, how many have a count of 10, and so on.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"When plotted, the frequency of occurence is plotted on the x-axis, and the number of k-mers on the y-axis.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"The k-mer spectra is composed of distributions representing groups of motifs at different frequencies in the sample, plus biases. Given not too many biases, the shape of the distributions provides a useful set of properties describing the biological sample, the sequencing process and the amount of useful data in the dataset.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"Let's take a look at a simple spectra. I'm going to use a small artificial dataset of a bacterial-like genome I generated using Pseudoseq.jl.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"I first created a random single chromosome genome, and then simulated a paired-end sequencing experiment to generate the kind of raw data we would get from a sequencing experiment.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"note: Note\nI'm doing this with a small artificial genome to be kind to the CI infrastructure that builds these docs!!!But you could follow along and try yourself with paired-end read data from a sequencing experiment from e.coli.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"Here I'm going to load in these reads, count the k-mers with the SerialMemCounter k-mer counter as the dataset is not that large, and then compute the 1D frequency spectra.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"using KmerAnalysis, ReadDatastores, BioSequences\n\nreads = @openreads \"../../../test/fakemicrobe.prseq\"\ncounter = serial_mem(DNAMer{31}, CANONICAL)\nkmer_counts = counter(reads) \n\nmy_spectra = spectra(kmer_counts)","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"I could have also created the counter and passed it to the spectra method along with the reads like so:","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"my_spectra = spectra(counter, reads)","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"At this point you can plot the k-mer frequency spectra using KmerAnalysisMakie.jl. You can visualise the k-mer frequency spectra quite easily, with several of Makie.jl's plotting primitives. By default a barplot will be drawn for a 1D k-mer spectra.","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"using KmerAnalysisMakie\n\nplot(my_spectra)","category":"page"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"note: Note\nLoading and using Makie without any pre-compilation of packages can be a bit slow as the JIT has to compile quite a lot.","category":"page"},{"location":"man/analyses/#Kmer-count-projection-1","page":"Analysing Kmer count data","title":"Kmer count projection","text":"","category":"section"},{"location":"man/analyses/#","page":"Analysing Kmer count data","title":"Analysing Kmer count data","text":"KmerAnalysis includes a k-mer count container type called IndexedCounts. Indexed counts are different to a vector or other container of MerCount, in that IndexedCounts store k-mer counts from many sources (often different read datasets) indexed against a reference dataset (usually a genome graph or set of contigs).","category":"page"},{"location":"api/KmerFrequencySpectra/#","page":"KmerFrequencySpectra","title":"KmerFrequencySpectra","text":"CurrentModule = KmerAnalysis","category":"page"},{"location":"api/KmerFrequencySpectra/#API:-KmerFrequencySpectra-1","page":"KmerFrequencySpectra","title":"API: KmerFrequencySpectra","text":"","category":"section"},{"location":"api/KmerFrequencySpectra/#Types-1","page":"KmerFrequencySpectra","title":"Types","text":"","category":"section"},{"location":"api/KmerFrequencySpectra/#","page":"KmerFrequencySpectra","title":"KmerFrequencySpectra","text":"KmerFrequencySpectra\nSpectraPeak","category":"page"},{"location":"api/KmerFrequencySpectra/#KmerAnalysis.KmerFrequencySpectra","page":"KmerFrequencySpectra","title":"KmerAnalysis.KmerFrequencySpectra","text":"A type for storing a frequency histogram of MerCounts, also sometimes referred to as a Kmer spectra. \n\n\n\n\n\n","category":"type"},{"location":"api/KmerFrequencySpectra/#KmerAnalysis.SpectraPeak","page":"KmerFrequencySpectra","title":"KmerAnalysis.SpectraPeak","text":"A simple struct describing the location and significance (or persistence) of a peak detected in kmer frequency spectra data.\n\nIt has the following fields:\n\npeak: The k-mer frequency at which the peak is at its highest.\nleft: The k-mer frequency constituting the left bound of the peak.\nright: The k-mer frequency constituting the right bound of the peak.\npersistence: The persistence of the peak in persistent homology analysis - essentially a significance value for the peak. Higher value = more significant.\n\n\n\n\n\n","category":"type"},{"location":"api/KmerFrequencySpectra/#Public-/-Safe-methods-1","page":"KmerFrequencySpectra","title":"Public / Safe methods","text":"","category":"section"},{"location":"api/KmerFrequencySpectra/#","page":"KmerFrequencySpectra","title":"KmerFrequencySpectra","text":"spectra\nfind_spectra_peaks","category":"page"},{"location":"api/KmerFrequencySpectra/#KmerAnalysis.spectra","page":"KmerFrequencySpectra","title":"KmerAnalysis.spectra","text":"spectra(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}\n\nBuild a 1 dimensional kmer frequency spectra, from a vector of kmer counts, excluding any kmer counts that don't meet min_count.\n\n\n\n\n\nspectra(xfreqs::Vector{MerCount{M}}, yfreqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}\n\nBuild\n\n\n\n\n\nspectra(x::Pair{I,C}, min_count::Integer = 0) where {I,C<:AbstractKmerCounter}\n\nBuild a 1 dimensional kmer frequency spectra, using a pair of (input data => kmer counter) as the input argument x. The spectra function uses the kmer counter to get the counts from the input data, before computing the spectra.\n\nAny kmer counts that don't meet min_count are excluded.\n\n\n\n\n\n","category":"function"},{"location":"api/KmerFrequencySpectra/#KmerAnalysis.find_spectra_peaks","page":"KmerFrequencySpectra","title":"KmerAnalysis.find_spectra_peaks","text":"find_spectra_peaks(s::KmerFrequencySpectra{1})\n\nAutomatically detect all the peaks present in the signal of a 1D k-mer frequency spectra, using a persistent homology method.\n\nReturns a vector of SpectraPeak.\n\n\n\n\n\n","category":"function"},{"location":"api/KmerFrequencySpectra/#Internal-/-Unsafe-methods-1","page":"KmerFrequencySpectra","title":"Internal / Unsafe methods","text":"","category":"section"},{"location":"api/KmerFrequencySpectra/#","page":"KmerFrequencySpectra","title":"KmerFrequencySpectra","text":"build_persistent_homology","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"CurrentModule = KmerAnalysis","category":"page"},{"location":"man/basic_counting/#Basic-types-and-basic-kmer-counting-1","page":"Basic types and counting","title":"Basic types and basic kmer counting","text":"","category":"section"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"MerCounting has a few basic types and methods to allow you do easily do some basic kmer counting.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"The first and perhaps most obvious of these is a type to represent a kmer and its frequency, for this, MerCounting provides the MerCount type.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"If you have a MerCount variable you can get the kmer value or the frequency value using the mer and freq getter methods.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"A vector of these MerCounts constitute one option for a simple data structure for representing kmer frequency data. This packages defines other more dedicated types, but we will get to those in later sections of this manual.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Ok, so let's do some very basic kmer counting for a sequence!","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"First we need a sequence:","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"using KmerAnalysis\nusing BioSequences\ns = randdnaseq(50)","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Ok, let's say we wanted to count the 7-mers, we can collect all the 7-mers using BioSequences' kmer iterator:","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"sevenmers = collect(each(DNAMer{7}, s))\nsevenmers[1]","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Ok BioSequences' kmer iterator on its own yields a MerIterResult, which contains a tuple of (Base position in sequence, kmer, reverse complement of the kmer).","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Now let's say I was only interested in counting canonical kmers in s. For any kmer and its reverse complement, the canonical version is the one that is lexicographically less.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"So let's revise our collect above to give us a vector of the canonical kmers. Thankfully, BioSequences comes with a method called canonical, for exactly this purpose.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"sevenmers = collect(canonical(x) for x in each(DNAMer{7}, s))","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"If you wanted to collect all the kmers in s \"as is\" as it were, you can use the fwmer method provided with BioSequences for this purpose.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"collect(fwmer(x) for x in each(DNAMer{7}, s))","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"This is quite nice and terse code, but we haven't used KmerAnalysis to help us at all yet. MerCounting allows us to achieve the above with the use of a single function called collect_mers, instead of collecting a generator:","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"sevenmers = collect_mers(CANONICAL, each(DNAMer{7}, s))","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"There are a few different methods of collect_mers for different types of input, check out the API reference for more info.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"MerCounting exports CANONICAL and NONCANONICAL, which work with collect_mers and other MerCounting functions to dictate how the kmers are produced.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Ok so we have our canonical kmers, how can we count them? Well one of the ways to count the number of distinct values in a vector is to sort it and traverse it in one go. The collapse_into_counts! methods do this for you.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"note: Note\ncollapse_into_counts! will sort the input sevenmers vector in place.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"mercounts = collapse_into_counts!(sevenmers)","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Looks like every distinct kmer appeared once, with the possible odd exceptions. Well that's to be expected! We did just do this with a single 50bp long random sequence after all!","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"So for any sequence or input we want to count kmers for, we have to collect_mers and then collapse them into counts with collapse_into_counts!.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"For convenience and more terse code, MerCounting provides a constructor for vectors of MerCount that does this for you.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"v = Vector{MerCount{DNAMer{7}}}(CANONICAL, s)","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"That's about all there is for doing basic kmer counting over individual biological sequences in julia.","category":"page"},{"location":"man/basic_counting/#","page":"Basic types and counting","title":"Basic types and counting","text":"Furthur on in this manual we will cover some dedicated kmer counting algorithms MerCounting provides for larger datasets like those contained in reads (using ReadDatastores.jl). We will also cover dedicated types such as kmer frequency spectra, and types useful for projecting kmer coverage and other kmer related metrics over a genome sequence.","category":"page"},{"location":"api/MerCount/#","page":"MerCount","title":"MerCount","text":"CurrentModule = KmerAnalysis","category":"page"},{"location":"api/MerCount/#API:-MerCount-1","page":"MerCount","title":"API: MerCount","text":"","category":"section"},{"location":"api/MerCount/#Types-1","page":"MerCount","title":"Types","text":"","category":"section"},{"location":"api/MerCount/#","page":"MerCount","title":"MerCount","text":"MerCount\nDNAMerCount\nRNAMerCount","category":"page"},{"location":"api/MerCount/#KmerAnalysis.MerCount","page":"MerCount","title":"KmerAnalysis.MerCount","text":"A simple mer count struct.\n\nMerCount is a simple struct that binds a mer value to a count of the number of times it has been observed.\n\nThis type, (sorted) vectors of them, and some additional utility methods, can form the basic building blocks of the higher-level mer counting functionality. This struct can also be an eltype of kmer hashes or more involved specialized types that store counts of Kmers.\n\nnote: Note\nThe count is stored as an UInt8 because often once the count is more than 255 we hardly care anymore.\n\n\n\n\n\n","category":"type"},{"location":"api/MerCount/#KmerAnalysis.DNAMerCount","page":"MerCount","title":"KmerAnalysis.DNAMerCount","text":"Shorthand for MerCount{DNAMer{K}}\n\n\n\n\n\n","category":"type"},{"location":"api/MerCount/#KmerAnalysis.RNAMerCount","page":"MerCount","title":"KmerAnalysis.RNAMerCount","text":"Shorthand for MerCount{RNAMer{K}}\n\n\n\n\n\n","category":"type"},{"location":"api/MerCount/#Public-/-Safe-methods-1","page":"MerCount","title":"Public / Safe methods","text":"","category":"section"},{"location":"api/MerCount/#","page":"MerCount","title":"MerCount","text":"mer\nfreq\nmerge\ncollect_mers\ncollect_mers!\ncollapse_into_counts\ncollapse_into_counts!\nmerge_into!\ncollapse!","category":"page"},{"location":"api/MerCount/#KmerAnalysis.mer","page":"MerCount","title":"KmerAnalysis.mer","text":"Get the mer from a MerCount.\n\n\n\n\n\n","category":"function"},{"location":"api/MerCount/#KmerAnalysis.freq","page":"MerCount","title":"KmerAnalysis.freq","text":"Get the count from a MerCount.\n\n\n\n\n\nGet the count from a MerCount, and convert it to type R.\n\n\n\n\n\n","category":"function"},{"location":"api/MerCount/#KmerAnalysis.collapse_into_counts!","page":"MerCount","title":"KmerAnalysis.collapse_into_counts!","text":"collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}\n\nBuild a vector of sorted MerCounts from a Vector of a mer type.\n\nThis is a basic kernel function used for any higher level and more complex kmer counting procedures.\n\nThis is like collapse_into_counts, except it's first argument is a result vector that is cleared and filled with the result.\n\nnote: Note\nThe input vector mers will be sorted by this method.\n\n\n\n\n\ncollapse_into_counts(mers::Vector{M}) where {M<:AbstractMer}\n\nBuild a vector of sorted MerCounts from a Vector of a mer type.\n\nThis is a basic kernel function used for any higher level and more complex kmer counting procedures.\n\nnote: Note\nThe input vector mers will be sorted by this method.\n\n\n\n\n\n","category":"function"},{"location":"api/MerCount/#KmerAnalysis.merge_into!","page":"MerCount","title":"KmerAnalysis.merge_into!","text":"merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}\n\nMerge the MerCounts from vector b into the vector a.\n\nnote: Note\nThis will sort the input vectors a and b.\n\n\n\n\n\n","category":"function"},{"location":"api/MerCount/#Internal-/-Unsafe-methods-1","page":"MerCount","title":"Internal / Unsafe methods","text":"","category":"section"},{"location":"api/MerCount/#","page":"MerCount","title":"MerCount","text":"unsafe_collapse_into_counts!\nunsafe_merge_into!","category":"page"},{"location":"api/MerCount/#KmerAnalysis.unsafe_collapse_into_counts!","page":"MerCount","title":"KmerAnalysis.unsafe_collapse_into_counts!","text":"unsafe_collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}\n\nwarning: Warning\nThis method is internal and marked as unsafe because it assumes that the mers input vector is already sorted.\n\n\n\n\n\n","category":"function"},{"location":"api/MerCount/#KmerAnalysis.unsafe_merge_into!","page":"MerCount","title":"KmerAnalysis.unsafe_merge_into!","text":"unsafe_merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}\n\nMerge the MerCounts from vector b into the vector a.\n\nwarning: Warning\nThis method is marked as unsafe as it assumes both of the input vectors a and b are already sorted.\n\n\n\n\n\n","category":"function"},{"location":"#","page":"Home","title":"Home","text":"(Image: Latest release) (Image: MIT license)  (Image: Stable documentation) (Image: Pkg Status) (Image: Chat)","category":"page"},{"location":"#Description-1","page":"Home","title":"Description","text":"","category":"section"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"KmerAnalysis is currently in pre-alpha development. But you can clone KmerAnalysis from the Julia REPL:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"julia> Pkg.clone(\"https://github.com/BioJulia/KmerAnalysis.jl.git\")","category":"page"},{"location":"#Testing-1","page":"Home","title":"Testing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"KmerAnalysis is tested against Julia 1.X on Linux, OS X, and Windows.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Latest build status:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: ) (Image: )","category":"page"},{"location":"#Contributing-1","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"We appreciate contributions from users including reporting bugs, fixing issues, improving performance and adding new features.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Take a look at the contributing files detailed contributor and maintainer guidelines, and code of conduct.","category":"page"},{"location":"#Financial-contributions-1","page":"Home","title":"Financial contributions","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"We also welcome financial contributions in full transparency on our open collective. Anyone can file an expense. If the expense makes sense for the development of the community, it will be \"merged\" in the ledger of our open collective by the core contributors and the person who filed the expense will be reimbursed.","category":"page"},{"location":"#Backers-and-Sponsors-1","page":"Home","title":"Backers & Sponsors","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Thank you to all our backers and sponsors!","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Love our work and community? Become a backer.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: backers)","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Does your company use BioJulia? Help keep BioJulia feature rich and healthy by sponsoring the project Your logo will show up here with a link to your website.","category":"page"},{"location":"#","page":"Home","title":"Home","text":"(Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: ) (Image: )","category":"page"},{"location":"#Questions?-1","page":"Home","title":"Questions?","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"If you have a question about contributing or using BioJulia software, come on over and chat to us on Gitter, or you can try the Bio category of the Julia discourse site.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"CurrentModule = KmerAnalysis","category":"page"},{"location":"man/counting_reads/#Counting-k-mers-in-read-datasets-1","page":"Counters for read datasets","title":"Counting k-mers in read datasets","text":"","category":"section"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"MerCounting provides some more dedicated counting algorithms for counting k-mers in read datasets, particularly ReadDatastores.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"This is because as conceptually simple as counting k-mers is, to do it quickly for large sequencing read datasets output by sequencing machines can be difficult for some datasets.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"There are many ways you can try to optimise the k-mer counting process, and many k-mer counting tools already exist.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"MerCounting provides a Counters submodule, which contains nessecery methods and types required to implement various kinds of k-mer counter, as well as exporting a selection of \"off-the-shelf\" methods that use different counting strategies, one or several of which hopefully suit the user's dataset and computational resources available. We describe and showcase these below.","category":"page"},{"location":"man/counting_reads/#The-SerialMemCounter-counter-1","page":"Counters for read datasets","title":"The SerialMemCounter counter","text":"","category":"section"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"The SerialMemCounter counter is the simplest k-mer counter for reads that is provided.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"It counts k-mers in each read serially, and stores the k-mers and counts entirely in RAM.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"SerialMemCounter","category":"page"},{"location":"man/counting_reads/#KmerAnalysis.SerialMemCounter","page":"Counters for read datasets","title":"KmerAnalysis.SerialMemCounter","text":"SerialMemCounter{M<:AbstractMer,C<:CountMode}\n\nKmerAnalysis' simplest kmer counter.\n\nBuild a sorted list (vector) of kmer counts (MerCount), serially and entirely in memory.\n\nBuild a SerialMemCounter using the serial_mem method. This returns a SerialMemCounter.\n\nA SerialMemCounter can be treated as a function or functor and passed to other functions as an argument, and it can be called on other arguments.\n\nIt is a general purpose kmer counting method:\n\nIt basically collects all the kmers from iterating over the input, before sorting the collected kmers and collapsing them into counts. To reduce allocations and  GC from many repeated calls of a SerialMemCounter, they keep a few buffers that are used often for kmer collection and count collapsing. \n\nWhilst general purpose for a variety of argument types, it is most suitable for:\n\nCounting the kmers in a small number of long-ish BioSequences (e.g. a reference genome or low(ish)-complexity genome graph.).\nCounting the kmers in a ReadDatastore if the dataset size is small enough such that all the kmers collected from every dataset can fit in a single machines memory. If you have a decent MacBook Pro, for example, datasets like E. coli paired end sequencing reads will be no problem.\nOther larger ReadDatastores IF you have the RAM to throw at it, e.g. you're on a HPC machine configured for large memory jobs.\n\n\n\n\n\n","category":"type"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"It is fairly simple to use. First you need a ReadDatastore, if you need to recall how to create one then head over to here in ReadDatastores.jl's documentation.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"The example below opens a datastore before creating a SerialMemCounter with serial_mem, to count the k-mers in the read datastore.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"using ReadDatastores\nusing FASTX\nusing BioSequences\nfwq = open(FASTQ.Reader, \"../../../test/ecoli_tester_R1.fastq\")\nrvq = open(FASTQ.Reader, \"../../../test/ecoli_tester_R2.fastq\")\nPairedReads{DNAAlphabet{2}}(fwq, rvq, \"ecoli-test-paired\", \"my-ecoli-test\", 250, 300, 0, FwRv)","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"using KmerAnalysis, ReadDatastores\nds = @openreads \"ecoli-test-paired.prseq\"\nc = serial_mem(DNAMer{31}, CANONICAL)\nkl = c(ds)","category":"page"},{"location":"man/counting_reads/#The-dist_mem-counter-1","page":"Counters for read datasets","title":"The dist_mem counter","text":"","category":"section"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"The DistMemCounter counter is a multi-process version of SerialMemCounter.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"DistMemCounter","category":"page"},{"location":"man/counting_reads/#KmerAnalysis.DistMemCounter","page":"Counters for read datasets","title":"KmerAnalysis.DistMemCounter","text":"The multi-process parallel version of the SerialMemCounter k-mer counter.\n\nLike the SerialMemCounter this counter does all counting and processing in memory, so it is not designed to limit memory use, but if using a cluster of machines the memory used per worker machine would be less. Of course the data is copied back over to the master process and the results combined and so the master node has to have enough RAM available to hold all the result.\n\n\n\n\n\n","category":"type"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"It is used in a very similar fashion to serial_mem, but you need to have multiple julia worker processes available to gain any benefit.","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"# Spin up some extra julia worker processes, make sure they are using the\n# appropriate Project.toml file for your project.\nusing Distributed\naddprocs(8, exeflags=\"--project=../../../docs/\")\n\n# Make sure the packages are loaded on all the workers.\n@everywhere using KmerAnalysis, ReadDatastores, BioSequences\n\nds = @openreads \"ecoli-test-paired.prseq\"\ndmc = dist_mem(DNAMer{31}, CANONICAL)\nkl = dmc(ds)","category":"page"},{"location":"man/counting_reads/#","page":"Counters for read datasets","title":"Counters for read datasets","text":"Naturally, for such a small and simple example as this, using DistMemCounter for such a small dataset is probably just worse than doing it in memory serially, because of the overheads. For real datasets you should see a benefit.","category":"page"},{"location":"man/installation/#Package-installation-1","page":"Installation","title":"Package installation","text":"","category":"section"},{"location":"man/installation/#","page":"Installation","title":"Installation","text":"KmerAnalysis is a component of the GenomeGraphs framework for graph based genome assembly and analysis. So if you've installed GenomeGraphs, you should already have it.","category":"page"},{"location":"man/installation/#","page":"Installation","title":"Installation","text":"However, KmerAnalysis has been designed to be generally useful in other settings and projects too. So it can be installed on its own as well.","category":"page"},{"location":"man/installation/#","page":"Installation","title":"Installation","text":"You can install KmerAnalysis from the julia REPL. Press ] to enter pkg mode again, and enter the following:","category":"page"},{"location":"man/installation/#","page":"Installation","title":"Installation","text":"pkg> add KmerAnalysis","category":"page"},{"location":"man/installation/#","page":"Installation","title":"Installation","text":"If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.","category":"page"}]
}
