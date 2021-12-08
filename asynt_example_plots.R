###############################################################################
########################  single scaffold alignment plot  #####################
###############################################################################
# We need to import the functions in the asynt.R script.
# You need to have the Intervals package installed on your system for this to work
source("asynt.R")

#if we have a single genomic region we are interested in, we can visualise the alignments diretcly
# First import the alignment data
alignments <- import.paf("examples/DplexMex_Dchry2.2_minimap2.asm20.paf.gz")

# note that you could also have used an alignment generated using mummer
# with the nucmer and show-coords tools), using import.nucmer()

# Now we subset by scaffold to just the  reference and query sequences we are interested in.
alignments <- subset(alignments, Rlen >= 100 & query=="contig30.1" & reference == "mxdp_29")

#and make the plot
par(mar = c(2,0,2,0))
plot.alignments(alignments)

#we can make the plot look a bit more fancy by using sigmoid lines
plot.alignments(alignments, sigmoid=TRUE)

#focus in on a specific region by setting the first and last base in the reference and query
plot.alignments(alignments, sigmoid=TRUE, Rfirst=1500000, Rlast = 1800000, Qfirst=1850000, Qlast=2150000)


###############################################################################
######################  multiple scaffold alignment plot  #####################
###############################################################################
# We need to import the functions in the asynt.R script.
# You need to have the Intervals package installed on your system for this to work
source("asynt.R")

# If we have multiple scaffolds making up a chromosome,
# we can string them together, either automatically or manually

#import alignments
alignments <- import.paf("examples/DplexMex_Dchry2.2_minimap2.asm20.paf.gz")

#Next we import load scaffold length data which is necessary to plot the scaffolds
ref_data <- import.genome(fai_file="examples/dplex_mex.fa.fai")
query_data <- import.genome(fai_file="examples/Dchry2.2.fa.fai")

#now define the scaffolds we're interested in
reference_scafs <- "mxdp_9"
query_scafs <- c("contig4.1", "contig4.2")

#keep only alignments involving these scaffolds
alignments <- subset(alignments, query %in% query_scafs & reference %in% reference_scafs)

#plot alinments
par(mar = c(4,4,4,0), xpd=NA)
plot.alignments.multi(alignments, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T)


###############################################################################
###############  multiple scaffold plot of synetny blocks  ####################
###############################################################################

#Alternatively, we can simplify the alignment by identifying synetny blocks (adjacent alignemnts in the same orientation)
synblocks <- get.synteny.blocks.multi(alignments, min_subblock_size=200)

#PRO TIP:
# You can make your syneeny blocks bigger and cleaner
# by running this command iteratively with increasing minimum subblock size.
# This will discard short alignment overlaps and small inversions and
# thereby increase the size of inferred syntenic blocks.
# Try this by uncommenting the lines below:

# synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
# synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)

plot.alignments.multi(synblocks, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T)

# We can manually change the orientation of the reference scaffolds
# and the plotting function will automatically flip the query scaffolds to optimise the synteny

reference_ori <- c("mxdp_9"="-")

plot.alignments.multi(synblocks, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T,
                      reference_ori=reference_ori)

#or force the orientation by telling it not to reverse query scaffolds

plot.alignments.multi(synblocks, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len, sigmoid=T,
                      reference_ori=reference_ori, no_reverse=TRUE)

###############################################################################
############################  Diagonal dot plot  ##############################
###############################################################################
# We need to import the functions in the asynt.R script.
# You need to have the Intervals package installed on your system for this to work
source("asynt.R")

# We will start with a whole-genome diagnoal 'dot plot'

# First we import the alignment file. In this case we import .paf format, produced my minimap2.
# There are also options for importing blast table format (output format 6) and nucmer coords files
alignments <- import.paf("examples/DplexMex_Dchry2.2_minimap2.asm20.paf.gz")

#we can filter for only long alignments and long scaffolds, which is sometimes necessary when viewing things at a large scale
alignments_5k <- subset(alignments, Rlen>= 5000 & Qlen>= 5000)

# We also need to import information about the assembly contig lengths
# These muct be the reference and query genomes that you aligned
# Note that in blast the reference is called the target
ref_data <- import.genome(fai_file="examples/dplex_mex.fa.fai")
query_data <- import.genome(fai_file="examples/Dchry2.2.fa.fai")

#Now we can plot a diagonal alignment 'dot plot'
par(mar=c(1.5,1.5,0,0), xpd=NA)
plot.alignments.diagonal(alignments_5k, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len)

# The plotting function tries to order and orient the QUERY scaffold to maximise synteny
# But it is clear that some REFERENCE scaffold should also be reversed and reordered
# We can control this manually by providing an additional file that gives
# chromosome and orientation for each reference scaffold.
# We therefore reload the genome information with this data

# Note that we included an additional optional file that
# gives chromosome and orientation information for each scaffold
ref_data <- import.genome(fai_file="examples/dplex_mex.fa.fai", chrom_file = "examples/dplex_mex.scaf_chrom.txt")

par(mar=c(1.5,1.5,0,0), xpd=NA)
plot.alignments.diagonal(alignments_5k, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len,
                         reference_ori=ref_data$seq_ori)

# Finally, if we want to have complete control,
# we can manually set the order and orientation of both the referece and query scaffolds 
# We do this by providing a chromosome file for the quey as well

query_data <- import.genome(fai_file="examples/Dchry2.2.fa.fai", chrom_file = "examples/Dchry2.2.chrom.txt")

#we will plot this without labels, and instead add chromosome information manually
par(mar=c(1.5,1.5,0,0), xpd=NA)
plot.alignments.diagonal(alignments_5k, reference_lens=ref_data$seq_len, query_lens=query_data$seq_len,
                         reference_ori=ref_data$seq_ori, query_ori = query_data$seq_ori,
                         no_reverse=TRUE, no_reorder=TRUE, no_labels=TRUE)

#shading is used to highlight the chromosomes
rect(ref_data$chrom_offset+1, 0, ref_data$chrom_offset+ref_data$chrom_len, query_data$chrom_offset+1, col = c("#00000011", "#00000000"), border=NA)
rect(0, query_data$chrom_offset+1, ref_data$chrom_offset+1, query_data$chrom_offset+query_data$chrom_len, col = c("#00000011", "#00000000"), border=NA)
rect(ref_data$chrom_offset+1, query_data$chrom_offset+1, ref_data$chrom_offset+ref_data$chrom_len, query_data$chrom_offset+query_data$chrom_len, col = c("#00000011", "#00000000"), border=NA)

#chromosome numbers are added mid-way in each chromosome
mtext(1, at=ref_data$chrom_offset+ref_data$chrom_len/2, text=substr(names(ref_data$chrom_len), 4,5), cex=0.5, line=-1.5)
mtext(2, at=query_data$chrom_offset+query_data$chrom_len/2, text=substr(names(query_data$chrom_len), 4,5), cex=0.5, line=-1, las=2)

#and species names are added as axis labels
mtext(1, at = sum(ref_data$chrom_len)/2, text=expression(italic("Danaus plexippus")))
mtext(2, at = sum(query_data$chrom_len)/2, text=expression(italic("Danaus chrysippus")))

