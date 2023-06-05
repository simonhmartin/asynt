### Functions for visualising genome alignemnts and synteny
### simon.martin@ed.ac.uk

###############################################################################
############### This code depends on the Intervals R package ##################
###############################################################################

library(intervals)

###############################################################################
##################### Functions for importing data ############################
###############################################################################

#paf is what minimap2 outputs
import.paf <- function(file, additional_fields=NULL){
    #because paf from minimap2 can have variable numbers of culumns, we use scan and strsplit
    paflines = scan(file, character(), sep="\n")
    
    paflist = strsplit(paflines, split="\t")
    
    output <- data.frame(reference = sapply(paflist, function(line) line[6]),
                         Rstart = as.numeric(sapply(paflist, function(line) line[8]))+1,
                         Rend = as.numeric(sapply(paflist, function(line) line[9])),
                         query = sapply(paflist, function(line) line[1]),
                         Qstart =  sapply(paflist, function(line) ifelse(line[5] == "+", as.numeric(line[3])+1, as.numeric(line[4]))),
                         Qend =  sapply(paflist, function(line) ifelse(line[5] == "+", as.numeric(line[4]), as.numeric(line[3])+1)),
                         strand= sapply(paflist, function(line) line[5]),
                         matches = as.numeric(sapply(paflist, function(line) line[10])),
                         total = as.numeric(sapply(paflist, function(line) line[11])),
                         MQ = as.numeric(sapply(paflist, function(line) line[12])),
                         stringsAsFactors=FALSE)
    
    output$Rlen <- output$Rend - output$Rstart + 1
    output$Qlen <- abs(output$Qend - output$Qstart) + 1
    output$identity <- 100 * output$matches / output$total
    
#     output <- data.frame(reference = paf[,6], Rstart = paf[,8]+1, Rend = paf[,9], Rlen = paf[,9] - paf[,8],
#                          query = paf[,1], Qstart = ifelse(paf[,5] == "+", paf[,3]+1, paf[,4]), Qend = ifelse(paf[,5] == "+", paf[,4], paf[,3]+1),
#                          Qlen = paf[,4] - paf[,3], strand= paf[,5], identity = 100*paf[,10]/paf[,11], MQ=paf[,12], stringsAsFactors=FALSE)
    
    #if any additional values have been requested
    for (f in additional_fields){
        adlist <- lapply(paflist, function(line) parse_key_value(line[-(1:12)]))
        output[,f] <- unlist(lapply(adlist, function(adfields) adfields[f]))
        }
    
    output
    }

#function for parsing the key-value fields in a .paf file
parse_key_value <- function(strings){
    separated <- strsplit(strings, ":")
    output <- lapply(separated, function(x) x[3])
    types <- lapply(separated, function(x) x[2])
    
    names(output) <- lapply(separated, function(x) x[1])
    names(types) <- lapply(separated, function(x) x[1])
    
    for (n in names(output)){
        if (types[[n]] == "i" | types[[n]] == "f") mode(output[[n]]) <- "numeric"
        }
    
    output
    }
    


import.blast <- function(file){
    blast_results <- read.table(file, header = F, as.is=T)
    names(blast_results) <- c("query","reference","identity","len","mismatches","gaps","Qstart","Qend","Rstart","Rend","e","score")
    #ensure that query and reference names are loaded as character strings and not numbers
    for (header in c("query","reference")) blast_results[,header] <- as.character(blast_results[,header])
    blast_results
    }

import.nucmer <- function(file){
    nucmer_results <- read.table(file, as.is=T)
    if (ncol(nucmer_results) == 9){
        names(nucmer_results) <- c("Rstart", "Rend", "Qstart", "Qend", "Rlen", "Qlen", "identity", "reference", "query")
        }
    else names(nucmer_results) <- c("Rstart", "Rend", "Qstart", "Qend", "Rlen", "Qlen", "identity", "Rstrand", "Qstrand", "reference", "query")
    nucmer_results
    }

#genome data is imported using a file that has contig as first column and length as second
#for example the .fai file produced by samtools faidx
get.contig.lengths <- function(fai_file){
    fai <- read.table(fai_file, as.is=T)
    lengths <- fai[,2]
    names(lengths) <- as.character(fai[,1])
    lengths
    }

import.genome <- function(fai_file, chrom_file=NULL){
    seq_len <- get.contig.lengths(fai_file)
    
    seq_ori <- sapply(names(seq_len), function(n) "+")
    
    if (is.null(chrom_file) == FALSE){
        chrom <- read.table(chrom_file, as.is=T, row.names=1)
        seq_by_chrom <- sapply(unique(chrom[is.na(chrom[,1])==FALSE,1]), function(c) rownames(chrom)[which(chrom[,1] == c)])
        if (ncol(chrom) > 1) seq_ori[rownames(chrom)] <- chrom[,2]
        }
    else seq_by_chrom <- sapply(names(seq_len), function(n) n, simplify=F)
    
    chrom_len <- sapply(seq_by_chrom, function(seq_name) sum(seq_len[seq_name]))
    chrom_offset <- cumsum(chrom_len) - chrom_len
    
    seq_names <- unlist(seq_by_chrom, use.names=F)
    seq_names <- c(seq_names, names(seq_len)[!(names(seq_len) %in% seq_names)])
    
    list(seq_names=seq_names, seq_len=seq_len[seq_names], seq_by_chrom=seq_by_chrom, chrom_len=chrom_len, chrom_offset=chrom_offset, seq_ori=seq_ori[seq_names]) 
    }

###############################################################################
############## Functions for analysing and processing alignments ##############
###############################################################################

#total unique length of aligned sequence for each scaffold
get.total.unique.length <- function(alignments, use="reference"){
    if (nrow(alignments) ==0) return(0)
    if (use == "query"){
        starts = apply(alignments[,c("Qstart", "Qend")], 1, min)
        ends = apply(alignments[,c("Qstart", "Qend")], 1, max)
    }
    else{
        starts = apply(alignments[,c("Rstart", "Rend")], 1, min)
        ends = apply(alignments[,c("Rstart", "Rend")], 1, max)
    }
    sum(size(reduce(Intervals(cbind(starts,ends), type="Z"))))
}

get.query.aln.len <- function(alignments){
    sapply(unique(alignments$query), function(q) get.total.unique.length(subset(alignments, query==q), use="query"))
    }

get.reference.aln.len <- function(alignments){
    sapply(unique(alignments$reference), function(t) get.total.unique.length(subset(alignments, reference==t), use="reference"))
    }

get.query.aln.prop <- function(alignments, query_lens){
    query_aln_len <- get.query.aln.len(alignments)
    query_aln_len / query_lens[names(query_aln_len)]
    }

get.reference.aln.prop <- function(alignments, reference_lens){
    ref_aln_len <- get.reference.aln.len(alignments)
    ref_aln_len / reference_lens[names(ref_aln_len)]
    }

#infer the best orientation for a pair of sequences given their alignment start and end positions
infer.orientation <- function(Rstart, Rend, Qstart, Qend){
    Rpos = interleave(Rstart, Rend)
    Qpos = interleave(Qstart, Qend)
    weights = rep(abs(Rend-Rstart)+1 + abs(Qend-Qstart)+1, each=2)
    model = lm(Qpos~Rpos, weights=weights)
    if (model$coefficients["Rpos"] < 0) return("-")
    return("+")
    }

#manually reverse query positions in an alignment table (not often required because plotting functions can reverse on the fly)
reverse.queries <- function(alignments, query_lens, query_names=NULL){
    #this will flip query sequences and thereby reverse their alignments  
    output <- alignments
    
    if (is.null(query_names) == FALSE) rows <- which(alignments[,"query"] %in% query_names)
    else rows <- 1:nrow(alignments)
    
    output[rows,"Qstart"] <- query_lens[alignments[rows,"query"]] - alignments[rows,"Qstart"]
    
    output[rows,"Qend"] <- query_lens[alignments[rows,"query"]] - alignments[rows,"Qend"]
    
    output
    }

#manually reverse reference positions in an alignment table (not often required because plotting functions can reverse on the fly)
reverse.references <- function(alignments, reference_lens, ref_names=NULL){
    #this will flip reference sequences and thereby reverse their alignments  
    output <- alignments
    
    if (is.null(ref_names) == FALSE) rows <- which(alignments[,"reference"] %in% ref_names)
    else rows <- 1:nrow(alignments)
    
    output[rows,"Rstart"] <- reference_lens[alignments[rows,"reference"]] - alignments[rows,"Rstart"]
    
    output[rows,"Rend"] <- reference_lens[alignments[rows,"reference"]] - alignments[rows,"Rend"]
    
    output
    }

#function to make all reference orientations forward and reverse query orientations where necessary and also to reorder by reference position
tidy.alignments <- function(alignments){
    #find those where ref coords are in reverse
    idx <- which(alignments$Rend - alignments$Rstart < 0)
    #flip them
    alignments_tidy <- alignments
    alignments_tidy[idx,c("Rstart","Rend","Qstart","Qend")] <- alignments_tidy[idx,c("Rend","Rstart","Qend","Qstart")]
    
    #reorder by reference start
    for (reference in unique(alignments$reference)){
        idx <- which(alignments$reference == reference)
        new_idx <- idx[order(alignments_tidy[idx,"Rstart"])]
        alignments_tidy[idx,] <- alignments_tidy[new_idx,]
        }
    
    alignments_tidy
    }

get.alignment.depth <- function(alignments, plot=FALSE){
    Rints <- Intervals(t(apply(alignments[,c("Rstart","Rend")], 1, sort)), type="Z", closed=T)
    unique_ints <- get.unique.intervals(Rints)
    unique_ints <- unique_ints[order(unique_ints[,1]),]
    depth <- sapply(interval_overlap(unique_ints, Rints), length)
    if (plot==TRUE) plot.intervals(unique_ints, height=depth,rectangles=TRUE)
    data.frame(start=unique_ints[,1], end=unique_ints[,2], depth=depth)
    }


###############################################################################
################ Functions for infering syntenic blocks #######################
###############################################################################

### function to get synetinic blocks by merging adjacent alignments with the same orientation
get.synteny.blocks <- function(alignments, max_gap=1e5, min_block_size=1e4, min_subblock_size=100){
    #make sure reference alignments have lower number first and are sorted by position
    alignments <- tidy.alignments(alignments)
    #store names of reference and query and make sure there is only one of each
    reference = unique(alignments$reference)
    query = unique(alignments$query)
    if(length(query)!=1 | length(reference)!=1){
        print("Only one reference and query allowed. Try get.synteny.blocks.multi")
        return(NULL)
        }
    
    #record orientation
    ori <- ifelse(sign(alignments$Rend - alignments$Rstart) != sign(alignments$Qend - alignments$Qstart), "-", "+")
    
    #get intervals for reference and query alignments
    Rints <- Intervals(t(apply(alignments[,c("Rstart","Rend")], 1, sort)), type="Z", closed=T)
    Qints <- Intervals(t(apply(alignments[,c("Qstart","Qend")], 1, sort)), type="Z", closed=T)
    
    #get subblocks by chopping up alignment reference intervals
    subblocks <- get.unique.intervals(Rints)
    
    #remove very small subblocks
    subblocks <- subblocks[size(subblocks) >= min_subblock_size,]
    
    #assign subblocks back to their alignments
    aln_subblocks <- interval_overlap(Rints, subblocks)
    
    #get order of query alignments
    Qord <- order(Qints[,1])
    
    #remove from these indices any that lack any subblocks
    Qord <- Qord[sapply(aln_subblocks[Qord], length) >= 1]
    
    #reorder both the query alignments and the alignmnet subblocks
    Qints_ord <- Qints[Qord,]
    Qord_subblocks <- aln_subblocks[Qord]
    ori_ord <- ori[Qord]
    
    #merge consecutive alignments with consecutive subblocks
    merged_alns <- list(1)
    merged_aln_ori <- ori_ord[1]
    i=1 #the current set of merged alignments we are making
    j=2 #the current alignment we are considering
    while (j <= length(Qord_subblocks)){
        merge=FALSE
        #Only consider merging if alignment j it is not too far away from last alignemnt in merged alignments[[i]] (check distance in query alns and in reference subblocks)
        if (Qints_ord[j,1] - Qints_ord[tail(merged_alns[[i]],1),2] <= max_gap &
            min(distance_to_nearest(subblocks[Qord_subblocks[[j]],], subblocks[Qord_subblocks[[tail(merged_alns[[i]],1)]],])) <= max_gap) {
            #if orientation is both forward and this alignment's first subblock is after last of the previous alignment, add it
            if (merged_aln_ori[i] == "+") {
                #if (ori_ord[j] == "+" & min(Qord_subblocks[[j]]) == max(Qord_subblocks[[tail(merged_alns[[i]],1)]]) + 1) {
                if (ori_ord[j] == "+" & min(Qord_subblocks[[j]]) > max(Qord_subblocks[[tail(merged_alns[[i]],1)]])) {
                    merge=TRUE
                    }
                }
            #if orientation is both reverse and this alignment's last subblock is before first of the previous alignment, add it
            else {
                #if (ori_ord[j] == "-" & max(Qord_subblocks[[j]]) == min(Qord_subblocks[[tail(merged_alns[[i]],1)]]) - 1) {
                if (ori_ord[j] == "-" & max(Qord_subblocks[[j]]) < min(Qord_subblocks[[tail(merged_alns[[i]],1)]])) {
                    merge=TRUE
                    }
                }
            }
        #now check whether merge is necessary and otherwise make a separate alignment
        if (merge == TRUE) {
            merged_alns[[i]] <- c(merged_alns[[i]], j)
            j <- j+1
            } else {
            #if not, this will probably become a separate alignment
            #but first check that the last block was not too small
            if (max(Qints_ord[merged_alns[[i]],]) - min(Qints_ord[merged_alns[[i]],]) + 1 >= min_block_size){
                i = i+1
                merged_alns[[i]] <- j
                merged_aln_ori[i] <- ori_ord[j]
                j <- j+1
                } else {
                #if we get here, our alignment cannot be merged to the last one, but the last one is too small
                #so we will discard the last one
                merged_alns <- merged_alns[-i]
                merged_aln_ori <- merged_aln_ori[-i]
                #if it had been the first, then we just use the j'th alignment to start a new first block
                if (i==1){
                    merged_alns[[i]] <- j
                    merged_aln_ori[i] <- ori_ord[j]
                    j <- j+1
                    }
                #otherwise we step back to the previous block 
                else i <- i-1
                }
            }
        }
    
    #check whether the last block was too short and remove if necessary
    if (max(Qints_ord[merged_alns[[i]],]) - min(Qints_ord[merged_alns[[i]],]) + 1 < min_block_size) {
        merged_alns <- merged_alns[-i]
        merged_aln_ori <- merged_aln_ori[-i]
        }
    
    #number of synteny blocks
    n=length(merged_alns)
    
    #if there are no synblocks, return an empty data frame
    if (n == 0) {
        return(data.frame(reference=character(), query=character(),
                          Rstart=numeric(), Rend=numeric(),
                          Qstart=numeric(), Qend=numeric(), stringsAsFactors=FALSE))
        }
    
    #all subblocks for each merged alignment
    merged_aln_subblocks <- lapply(merged_alns, function(alns) subblocks[unlist(Qord_subblocks[alns]),])
    
    #min and max positions for all synblocks
    synblock_Rmin <- sapply(merged_aln_subblocks, min)
    synblock_Rmax <- sapply(merged_aln_subblocks, max)

    synblock_Qmin <- sapply(merged_alns, function(alns) min(Qints_ord[alns,]))
    synblock_Qmax <- sapply(merged_alns, function(alns) max(Qints_ord[alns,]))
    
    #get query starts and ends of merged alignments
    #get ref starts and ends for merged alignments (based on subblocks)
    synblocks <- data.frame(reference=rep(reference,n), query=rep(query,n),
                            Rstart=synblock_Rmin,
                            Rend=synblock_Rmax,
                            Qstart=ifelse(merged_aln_ori=="+",synblock_Qmin,synblock_Qmax),
                            Qend=ifelse(merged_aln_ori=="+",synblock_Qmax,synblock_Qmin),
                            stringsAsFactors=FALSE)
    
    synblocks
    }

#wrapper function for synblocks for when there are multiple referneces and/or queries
get.synteny.blocks.multi <- function(alignments, max_gap=1e5, min_block_size=1e4, min_subblock_size=10){
    #store names of references and queries
    references = unique(alignments$reference)
    queries = unique(alignments$query)
    
    synblocks <- data.frame(reference=character(), query=character(),
                            Rstart=numeric(), Rend=numeric(), Qstart=numeric(), Qend=numeric())
    
    for (reference in references){
        print(reference)
        for (query in queries){
            print(paste("   ", query))
            idx = which(alignments$reference==reference & alignments$query==query)
            if (length(idx) > 0) {
                synblocks <- rbind(synblocks, get.synteny.blocks(alignments[idx,], max_gap=max_gap, min_block_size=min_block_size, min_subblock_size=min_subblock_size))
                }
            }
        }
    
    synblocks
    }

#function that is a bit like reduce, except it chops intervals that overlap rather than merging them
#this is a key function underlying the inference of synteny blocks
get.unique.intervals <- function(ints){
    output = matrix(ncol=2,nrow=0)
    final_end <- max(ints[,2])
    current_start <- min(ints[,1])
    while(current_start <= final_end){
        next_start <- ifelse(any(ints[,1] > current_start), min(ints[ints[,1] > current_start, 1]), Inf)
        next_end <- min(ints[ints[,2] >= current_start, 2])
        if (next_start <= next_end){
            #we have to cut before this next start and start again
            output <- rbind(output, c(current_start, next_start-1))
            current_start <- next_start
            } else {
            #otherwise we use include the end and start at next position
            output <- rbind(output, c(current_start, next_end))
            current_start <- next_end + 1
            }
        }
    output <- Intervals(output[order(output[,1]),], type="Z", closed=T)
    #only return those included in the input
    output[unique(unlist(interval_included(ints, output))),]
    }


###############################################################################
########################### Functions for plotting ############################
###############################################################################

### function for plotting alignments in a parallel arrangement
plot.alignments <- function(alignments, Qfirst=NULL, Qlast=NULL, Rfirst=NULL, Rlast=NULL,
                            cols = c("#0000ff", "#ff0000"), colour_by = "orientation",
                            gap=0, show_outline=TRUE, sigmoid=FALSE, tick_spacing=100000, las=2, cex.axis=0.6,
                            min_colour_value=NA, max_colour_value=NA, lwd=NULL){
    
    if(length(unique(alignments$query)) !=1 | length(unique(alignments$reference)) != 1){
        print("Only one reference and query scaffold allowed. Try plot.alignments.multi() for multiple scaffolds")
        return(NULL)
        }
    
    if(is.null(Qfirst)) Qfirst <- min(c(alignments$Qstart, alignments$Qend))
    if(is.null(Qlast)) Qlast <- max(c(alignments$Qstart, alignments$Qend))
    if(is.null(Rfirst)) Rfirst <- min(c(alignments$Rstart, alignments$Rend))
    if(is.null(Rlast)) Rlast <- max(c(alignments$Rstart, alignments$Rend))
    
    plot_length = max(c(Qlast-Qfirst, Rlast-Rfirst)) + 1
    
    #offset all by first pos
    Qstarts = alignments$Qstart - Qfirst + 1
    Qends = alignments$Qend - Qfirst + 1
    
    Rstarts = alignments$Rstart - Rfirst + 1
    Rends = alignments$Rend - Rfirst + 1
    
    plot(0,cex = 0, xlim = c(1, plot_length), ylim = c(-0.1,1.1), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")
    
    if (colour_by == "identity"){
        border = colorRampPalette(cols)(10)[cut(c(alignments$identity, min_colour_value, max_colour_value), 10)]
        col = paste0(border,"50")
        }
    else{
        border = ifelse(sign(Rends - Rstarts) == sign(Qends - Qstarts), cols[1], cols[2])
        col = paste0(border,"50")
        }
    
    for (i in 1:nrow(alignments)){
        if (sigmoid == FALSE){
            polygon(c(Qstarts[i], Qends[i], Rends[i], Rstarts[i]),
                    c(1-gap,1-gap,0+gap,0+gap), col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
            }
        #curved lines
        else{
            lines.to.poly(sigmoid.connector(Qstarts[i], 1-gap, Rstarts[i], 0+gap, vertical=T),
                        sigmoid.connector(Qends[i], 1-gap, Rends[i], 0+gap, vertical=T),
                        col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
            }
        }
    
    segments(1, 1, Qlast-Qfirst+1, 1, lwd = 5)
    segments(1, 0, Rlast-Rfirst+1, 0, lwd = 5)
    
    mtext(text=alignments$query[1], side=3, at=(Qlast-Qfirst+1)/2,)
    mtext(text=alignments$reference[1], side=1, at=(Rlast-Rfirst+1)/2,)
    
    axis(3, at = (1:plot_length)[which(Qfirst:Qlast %% tick_spacing == 0)],
         labels = (Qfirst:Qlast)[which(Qfirst:Qlast %% tick_spacing == 0)], line = -3, las=las, cex.axis=cex.axis)
    axis(1, at = (1:plot_length)[which(Rfirst:Rlast %% tick_spacing == 0)],
         labels = (Rfirst:Rlast)[which(Rfirst:Rlast %% tick_spacing == 0)], line = -3, las=las, cex.axis=cex.axis)
    }





### Function to plot alignments from multiple scaffolds/chromosomes in a parallel arrangement

plot.alignments.multi <- function(alignments, reference_lens, query_lens, reference_ori=NULL, query_ori=NULL,
                                  only_show_aligned_seqs=TRUE, no_reverse=FALSE, no_reorder=FALSE,
                                  edge_width=0.3, chrom_width=0.1, gap=0, reference_above=FALSE,
                                  cols = c("#0000ff", "#ff0000"), show_connectors=TRUE, show_outline=TRUE, sigmoid=FALSE,
                                  colour_by = "orientation", min_colour_value=NA, max_colour_value=NA,
                                  lwd=NULL, show_contigs=TRUE, show_labels=TRUE, angle_labels=TRUE, labels_cex=0.7, labels_offset=0.02,
                                  centre=TRUE, plot_length = NULL,
                                  show_alignment_tracts=FALSE){
    
    ### sequences to include    
    if (only_show_aligned_seqs == TRUE){
        references <- names(reference_lens)[names(reference_lens) %in% unique(alignments$reference)]
        queries <-  names(query_lens)[names(query_lens) %in% unique(alignments$query)]
        }
    else {
        references <- names(reference_lens)
        queries <-  names(query_lens)
        }
    
    references_total_len <- sum(reference_lens[references])
    queries_total_len <- sum(query_lens[queries])
    
    if (is.null(plot_length) == TRUE) plot_length = max(references_total_len,queries_total_len)
    
    ### New refernece positions given input order and orientation
    reference_offsets <- cumsum(reference_lens[references]) - reference_lens[references]
    if (centre == TRUE) reference_offsets <- reference_offsets + floor((plot_length - references_total_len)/2)

    if (is.null(reference_ori)==TRUE) reference_ori <- sapply(references, function(x) "+")
    
    alignments$Rstart_new <- reference_offsets[alignments$reference] + ifelse(reference_ori[alignments$reference] == "-",
                                                                      reference_lens[alignments$reference] - alignments$Rstart,
                                                                      alignments$Rstart)
    
    alignments$Rend_new <- reference_offsets[alignments$reference] + ifelse(reference_ori[alignments$reference] == "-",
                                                                    reference_lens[alignments$reference] - alignments$Rend,
                                                                    alignments$Rend)
    
    ### Query order, orientation and offset
    alignments <- alignments[order(alignments$Qstart),]
    
    alns_by_query <- sapply(queries, function(query) alignments[which(alignments$query==query),], simplify=F)
    
    if (is.null(query_ori) == TRUE){
        if (no_reverse == FALSE){
            query_ori <- sapply(queries, function(query) infer.orientation(alns_by_query[[query]]$Rstart_new,
                                                                           alns_by_query[[query]]$Rend_new,
                                                                           alns_by_query[[query]]$Qstart,
                                                                           alns_by_query[[query]]$Qend))
            }
        else query_ori <- sapply(queries, function(x) "+")
        }
    
    #get the order for placing queries and then the offset
    if (no_reorder == FALSE){
        midpoints <- sapply(queries, function(query) get.median.from.intervals(alns_by_query[[query]][,c("Rstart_new","Rend_new")]))
        idx <- order(midpoints)
        }
    else idx <- 1:length(queries)
    
    query_offsets <- cumsum(query_lens[queries[idx]]) - query_lens[queries[idx]]
    if (centre == TRUE) query_offsets <- query_offsets + floor((plot_length - queries_total_len)/2)
    
    ### New query positions given order, orientation and offset
    alignments$Qstart_new <- query_offsets[alignments$query] + ifelse(query_ori[alignments$query] == "-",
                                                                      query_lens[alignments$query] - alignments$Qstart,
                                                                      alignments$Qstart)
    
    alignments$Qend_new <- query_offsets[alignments$query] + ifelse(query_ori[alignments$query] == "-",
                                                                    query_lens[alignments$query] - alignments$Qend,
                                                                    alignments$Qend)
    
    ### make plot

    if (reference_above == TRUE) ylim <- c(1+edge_width,-edge_width)
    else ylim <- c(-edge_width,1+edge_width)
    
    plot(0,cex = 0, xlim = c(1, plot_length), ylim = ylim, xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")
    
    if (show_contigs == TRUE){
        rect(reference_offsets[references]+1, 0, reference_offsets[references]+reference_lens[references], -chrom_width, border="gray40", col="gray90")
        rect(query_offsets[queries]+1, 1, query_offsets[queries]+query_lens[queries], 1+chrom_width, border="gray40", col="gray90")
        
        if (show_labels == TRUE){
            if (angle_labels == TRUE){
                angle=45
                adj_ref=ifelse(reference_above==TRUE,0,1)
                adj_qry=ifelse(reference_above==TRUE,1,0)
                }
            else {
                angle=0
                adj_ref=0.5
                adj_qry=0.5
                }
            
            text(reference_offsets[references]+reference_lens[references]/2, -chrom_width-labels_offset,
                 labels = ifelse(reference_ori[references] == "-", paste0(references,"*"), references),
                 cex = labels_cex, srt=angle, adj=adj_ref)
            text(query_offsets[queries]+query_lens[queries]/2, 1+chrom_width+labels_offset,
                 labels = ifelse(query_ori[queries] == "-", paste0(queries,"*"), queries),
                 cex = labels_cex, srt=angle, adj=adj_qry)
            }
        }
    
    if (colour_by == "identity"){
        border = colorRampPalette(cols)(10)[cut(c(alignments$identity, min_colour_value, max_colour_value), 10)]
        col = paste0(border,"50")
        }
    else{
        border = ifelse(sign(alignments$Rend_new-alignments$Rstart_new) == sign(alignments$Qend_new-alignments$Qstart_new), cols[1], cols[2])
        col = paste0(border,"50")
        }
    
    if (show_connectors==TRUE){
        for (i in 1:nrow(alignments)){
            if (sigmoid == FALSE){
                polygon(c(alignments$Qstart_new[i], alignments$Qend_new[i], alignments$Rend_new[i], alignments$Rstart_new[i]),
                        c(1-gap,1-gap,0+gap,0+gap), col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
                }
            #curved lines
            else{
                lines.to.poly(sigmoid.connector(alignments$Qstart_new[i], 1-gap, alignments$Rstart_new[i], 0+gap, vertical=T),
                            sigmoid.connector(alignments$Qend_new[i], 1-gap, alignments$Rend_new[i], 0+gap, vertical=T),
                            col = col[i], border=ifelse(show_outline==FALSE, NA, border[i]), lwd=lwd)
                }
            }
        }
    
    if (show_alignment_tracts == TRUE){
        rect(alignments$Qstart_new, 1, alignments$Qend_new, 1+chrom_width, col = border, border=NA)
        rect(alignments$Rstart_new, 0, alignments$Rend_new, -chrom_width, col = cols[1], border=NA)
        }
    
    }

### Function to plot alignments from multiple scaffolds/chromosomes in a diagonal (dot plot) arrangement

plot.alignments.diagonal <- function(alignments, reference_lens, query_lens, reference_ori=NULL, query_ori=NULL,
                                  only_show_aligned_seqs=TRUE, no_reverse=FALSE, no_reorder=FALSE,
                                  cols = c("#0000ff", "#ff0000"), colour_by = "orientation", min_colour_value=NA, max_colour_value=NA,
                                  lwd=NULL, no_labels=FALSE, angle_labels=TRUE, labels_cex=0.7, labels_offset=0, xmax=NULL, ymax=NULL){
    
    ### sequences to include    
    if (only_show_aligned_seqs == TRUE){
        references <- names(reference_lens)[names(reference_lens) %in% unique(alignments$reference)]
        queries <-  names(query_lens)[names(query_lens) %in% unique(alignments$query)]
        }
    else {
        references <- names(reference_lens)
        queries <-  names(query_lens)
        }
    
    refsum <- sum(reference_lens[references])
    querysum <- sum(query_lens[queries])
    
    ### New refernece positions given input order and orientation
    reference_offsets <- cumsum(reference_lens[references]) - reference_lens[references]

    if (is.null(reference_ori)==TRUE) reference_ori <- sapply(references, function(x) "+")
    
    alignments$Rstart_new <- reference_offsets[alignments$reference] + ifelse(reference_ori[alignments$reference] == "-",
                                                                      reference_lens[alignments$reference] - alignments$Rstart,
                                                                      alignments$Rstart)
    
    alignments$Rend_new <- reference_offsets[alignments$reference] + ifelse(reference_ori[alignments$reference] == "-",
                                                                    reference_lens[alignments$reference] - alignments$Rend,
                                                                    alignments$Rend)
    
    
    ### Query order, orientation and offset (only orient using best reference)
    alignments <- alignments[order(alignments$Qstart),]
    
    alns_by_query <- sapply(queries, function(query) alignments[which(alignments$query==query),], simplify=F)
    
    #get best reference match by query
    best_ref_by_query <- sapply(queries, function(query) names(which.max(get.reference.aln.len(alns_by_query[[query]])))[1])
    
    #for ordering and orienting, make a pruned down version containing best match only
    alns_by_query_BESTREF <- sapply(queries, function(query) alns_by_query[[query]][alns_by_query[[query]]$reference == best_ref_by_query[query],], simplify=F)
    
    if (is.null(query_ori) == TRUE){
        if (no_reverse == FALSE){
            query_ori <- sapply(queries, function(query) infer.orientation(alns_by_query_BESTREF[[query]]$Rstart_new,
                                                                           alns_by_query_BESTREF[[query]]$Rend_new,
                                                                           alns_by_query_BESTREF[[query]]$Qstart,
                                                                           alns_by_query_BESTREF[[query]]$Qend))
            }
        else query_ori <- sapply(queries, function(x) "+")
        }
    
    #get the order for placing queries and then the offset
    if (no_reorder == FALSE){
        midpoints <- sapply(queries, function(query) get.median.from.intervals(alns_by_query[[query]][,c("Rstart_new","Rend_new")]))
        idx <- order(midpoints)
        }
    else idx <- 1:length(queries)
    
    query_offsets <- cumsum(query_lens[queries[idx]]) - query_lens[queries[idx]]
    
    ### New query positions given order, orientation and offset
    alignments$Qstart_new <- query_offsets[alignments$query] + ifelse(query_ori[alignments$query] == "-",
                                                                      query_lens[alignments$query] - alignments$Qstart,
                                                                      alignments$Qstart)
    
    alignments$Qend_new <- query_offsets[alignments$query] + ifelse(query_ori[alignments$query] == "-",
                                                                    query_lens[alignments$query] - alignments$Qend,
                                                                    alignments$Qend)
    
    ### Plot
    
    plot(0,cex = 0, xlim = c(1, ifelse(is.null(xmax)==TRUE, refsum, xmax)), ylim = c(1,ifelse(is.null(ymax)==TRUE, querysum, ymax)), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n", xpd=FALSE)
    
    segments(c(reference_offsets, refsum), 0, c(reference_offsets, refsum), querysum,  col="gray70")
    segments(0, c(query_offsets, querysum), refsum, c(query_offsets, querysum),  col="gray70")
    
    if (no_labels == FALSE){
        if (angle_labels == TRUE){
            angle_ref=45
            angle_qry=45
            adj_ref=c(1,1)
            adj_qry=c(1,0)
            }
        else {
            angle_ref=0
            angle_qry=90
            adj_ref=NULL
            adj_qry=NULL
            }
        
        text(reference_offsets[references]+reference_lens[references]/2, 0-labels_offset,
             labels = ifelse(reference_ori[references]=="-", paste0(references,"*"), references), cex = labels_cex, srt=angle_ref, adj=adj_ref)
        text(0-labels_offset, query_offsets[queries]+query_lens[queries]/2, cex = labels_cex,
             srt=angle_qry, adj=adj_qry, labels = ifelse(query_ori[queries]=="-", paste0(queries,"*"), queries))
        }
    
    if (colour_by == "identity") {
        col = colorRampPalette(cols)(10)[cut(c(alignments$identity, min_colour_value, max_colour_value), 10)]
        }
    else{
        col = ifelse(sign(alignments$Rend_new-alignments$Rstart_new) == sign(alignments$Qend_new-alignments$Qstart_new), cols[1], cols[2])
        }
    
    for (i in 1:nrow(alignments)){
        segments(alignments$Rstart_new[i], alignments$Qstart_new[i], alignments$Rend_new[i], alignments$Qend_new[i], col = col[i], lwd=lwd)
        }
    
    output = list(queries=queries, query_ori=query_ori, query_offsets=query_offsets,
                  referecnes=references, reference_ori=reference_ori, reference_offsets=reference_offsets)
    
    invisible(output)
    }

###############################################################################
####################### Helper functions for plotting #########################
###############################################################################

# calculate the path of a sigmoid connector between two points
sigmoid.connector <- function(x1,y1,x2,y2, curvature=10, steps=50, vertical=FALSE){
    if (vertical==TRUE) {
        vals <- c(x1,y1,x2,y2)
        x1 <- vals[2]
        x2 <- vals[4]
        y1 <- vals[1]
        y2 <- vals[3]
        }
    x <- seq(x1,x2,(x2-x1)/steps)
    x_norm <- (x-(x1+x2)/2)/((x2-x1)/2)
    y_norm <- 1/(1+exp(-x_norm*curvature))
    y <- y1+(y_norm*(y2-y1))
    if (vertical==TRUE) return(cbind(y,x))
    cbind(x,y)
    }

#make a polygon from two lines
lines.to.poly <- function(l1,l2, col=NULL, border=NULL, lwd=NULL){
    polygon(c(l1[,1],rev(l2[,1])), c(l1[,2],rev(l2[,2])), col=col, border=border)
    }

#add scale bar to a plot
draw.scale.bar <- function(x,y, width, height, text, lwd=NULL, offset=NULL, cex=NULL){
    if (is.null(offset)==TRUE) offset<-height
    x_left = x - width/2
    x_right = x+width/2
    y_bot = y-height/2
    y_top = y+height/2
    segments(c(x_left, x_left, x_right), c(y_bot, y, y_bot), c(x_left, x_right, x_right), c(y_top, y, y_top), lwd=lwd)
    text(x,y+offset, labels=text, cex=cex)
    }

#add block arrows to a plot
draw.horizontal.block.arrow <- function(x1,x2,y, width=0, width_scaler=1, length_scaler=0.2, col=NULL, border=NULL){
    x_hinge <- x1 + (x2-x1)*(1-length_scaler)
    y1 <- y-width/2
    y2 <- y+width/2
    y1_hinge <- y1-width_scaler*width/2
    y2_hinge <- y2+width_scaler*width/2
    
    polygon(c(x1,x_hinge,x_hinge,x2,x_hinge,x_hinge,x1), c(y1,y1,y1_hinge,y,y2_hinge,y2,y2), col=col, border=border)
    }


plot.intervals <- function(ints, height=NULL, col="black", rectangles=FALSE, gap=1,
                           xlim=NULL, add=FALSE, return_height=FALSE){
    #sort by start
    ord <- order(ints[,1])
    ints <- ints[ord,]
    n = nrow(ints)
    
    if (is.null(height)==FALSE){
        if (length(height)==1) height <- rep(height, n)
        else height <- height[ord]
        }
    else{
        #if no height value provided, set it so as to distinguish intervals
        offsets <- rep(-gap, n)
        height <- rep(1, n)
        for (i in 1:n){
            while (ints[i,1] < offsets[height[i]] + gap) height[i] <- height[i] + 1
            #now we found an appropriate height value. set the new offset
            offsets[height[i]] <- ints[i,2]
            }
        }
    
    if (add==FALSE){
        if(is.null(xlim) == TRUE) xlim <- c(ints[1,1], max(ints[,2]))
        plot(0,cex=0,xlim=xlim, ylim = c(0, max(height)+1), ylab="depth", xlab="position")
        }
    
    if (rectangles==TRUE) rect(ints[,1], 0, ints[,2], height, lwd=0, border=NA, col=col)
    else segments(ints[,1], height, ints[,2], height, col=col)
    if (return_height) return(height)
    }

sim.intervals <- function(n=1,size_mean=10, size_sd=10, size_min=5, max_start=100){
    starts <- sort(round(runif(n,1,max_start)))
    sizes <- round(rnorm(n, size_mean, size_sd)) -1
    sizes[sizes < size_min] <- size_min
    ends <- starts + sizes -1
    Intervals(cbind(starts,ends), type="Z", closed=T)
    }


#function get get a median from intervals, accounting for overlaps
get.median.from.intervals <- function(intervals){
    #make sure intervals are sorted by
    if (nrow(intervals) == 0) return(NA)
    intervals_sorted <- t(apply(intervals, 1, sort))
    red_intervals <- close_intervals(reduce(Intervals(intervals_sorted, type="Z")))
    sizes <-  size(red_intervals)
    total = sum(sizes)
    size_offsets <- cumsum(sizes) - sizes
    median_index = floor(total/2)
    interval_containing_median <- tail(which(size_offsets < median_index),1)
    
    #if total is odd, median is an integer inside an interval
    if (total %% 2 == 1) return(red_intervals[interval_containing_median,1] + median_index - size_offsets[interval_containing_median])
    
    #if we get here, total is even, first check if median is between two intervals
    if (median_index == size_offsets[interval_containing_median] + sizes[interval_containing_median]){
        return(0.5*(red_intervals[interval_containing_median,2] + red_intervals[interval_containing_median+1,1]))
        }
    
    #otherwise it is inside the interval, but add a half to the final value
    red_intervals[interval_containing_median,1] + median_index - size_offsets[interval_containing_median] - 0.5
    }

#interleave two two equal length vectors into one
interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }

###############################################################################
######################## Additional useful functions ##########################
###############################################################################

#get chromosome positions by stringing scaffolds together
get.chrom.pos <- function(scaffold, position, scaf_len, scaf_ori){
    
    npos <- length(position)
    
    offsets <- cumsum(scaf_len) - scaf_len
    
    chrom_pos <- numeric(length = npos)
    
    for (i in 1:npos){
        if (scaf_ori[scaffold[i]] == "+") chrom_pos[i] <- offsets[scaffold[i]] + position[i]
        else chrom_pos[i] <- offsets[scaffold[i]] + scaf_len[scaffold[i]] - position[i]
        }
    
    chrom_pos
    }


#a function to get the order of a vector given an input vector with the order we want
# so if we start with c("A","A","A","B","B","C") and the order we want is c("C","A","B")
# we will get c("C","B","B","A","A","A")
#to be honest I can't remember why I originally wrote this, but I'm leaving it in becasue it seems useful
order.by.template <- function(values, template){
    template_idx <- 1:length(template)
    names(template_idx) <- template
    order(template_idx[values])
    }

#asembly statistics
get.N50_L50 <- function(contig_lengths){
    lengths_sorted = sort(contig_lengths, decreasing=TRUE)
    cs = cumsum(lengths_sorted)
    L50 = which(cs >= sum(contig_lengths)/2)[1]
    N50 = lengths_sorted[L50]
    return(c(N50=as.numeric(N50), L50=as.numeric(L50)))
    }


