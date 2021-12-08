
# Asynt: R functions for exploring synteny using whole genome alignments

* Make diagonal 'dot' plots
* Plot alignment tracts between a pair of genomes
* Merge adjacent alignments into synteny 'blocks' for cleaner plots

See our [preprint](https://doi.org/10.1101/2021.12.06.471392) for examples of the plots you can make (see Figures 1 and S1)

## How to use this code

If you already have alignment coordinate files from minimap2 (recommended) or mummer (using nucmer and show-coords), you are ready to go.
Open the script `asynt_example_plots.R` in an interactive R session (e.g. Rstudio) and work through it line by line to explore the kinds of plots you can make.

## Where do I get alignments from?

Make alignemnts between two assemblies (or a single assembly) using [minimap2](https://github.com/lh3/minimap2) or [mummer](https://mummer4.github.io/)

Here is an example command for minimap2:

`minimap2 -x asm20 reference.fa query.fasta | gzip > mm2asm20.paf.gz`

`-x asm20` uses presets suited for genomes up to 20% divergent.

## How does asynt infer synteny blocks?

There are some sofisticated tools that use probabilistic approaches for infering synteny blocks. This is not one of those.

The algorithm has three steps:
1. Alignments are split into ‘sub-blocks’ that each correspond to a unique tract of the reference assembly.
2. Sub-blocks below a minimum size are discarded.
3. Adjacent sub-blocks that are in the same orientation and are below some threshold distance apart are merged to yield syntenic blocks.

These three steps can be performed iteratively to first identify regions of fine-scale synteny and build these up into larger syntenic blocks (discarding short overlaps, small inversions etc).
The nature of this approach means that you will get a different result depending on what you use as the reference. If possible, use a reference that represents the ancestral state', such that your query genome is being represented as a new arrangement of ancestral blocks.








