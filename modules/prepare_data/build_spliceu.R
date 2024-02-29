#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# first fasta, second gtf, third output

genome_path = args[1]
gtf_path = args[2]
output_path = args[3]
threads = as.integer(args[4])
min_length = as.integer(args[5])
max_mismatch = as.integer(args[6])

#gtf_path = "/fs/nexus-projects/sc_frag_len/nextflow/data/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#genome_path = "/fs/nexus-projects/sc_frag_len/nextflow/data/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
#output_path = "refdata-gex-GRCh38-2020-A"
# output_path = "refdata-gex-mm10-2020-A"
#output_path = "/fs/nexus-projects/sc_frag_len/nextflow/test/build_spliceu"
#threads = 30
#min_length = 6
#max_mismatch = 0

dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(BSgenome)
    library(GenomicRanges)
    library(doParallel)
    library(eisaR)
    library(Biostrings)
    library(GenomicFeatures)
    library(rtracklayer)
})

# Find the polyA and polyT sites
threads=min(threads,detectCores())
registerDoParallel(cores=threads)

# initialize the grlfull
grlfull <- GenomicRanges::GRangesList()

genome <- Biostrings::readDNAStringSet(file.path(genome_path))
# get the first word as the name
names(genome) <- stringr::word(names(genome), 1)


# chr_polya = IRanges::reduce(
#     as(Biostrings::vmatchPattern(
#         paste0(
#             rep("A", min_length),
#             collapse = ""
#         ),
#         genome,
#         max.mismatch = max_mismatch
#     ),
#     "GRanges"
#     )
# )
# export.bed(chr_polya, ignore.strand = TRUE,con=file.path(output_path, "chr_polya.bed"))


# chr_polyt = IRanges::reduce(
#     as(Biostrings::vmatchPattern(
#         paste0(
#             rep("A", min_length),
#             collapse = ""
#         ),
#         genome,
#         max.mismatch = max_mismatch
#     ),
#     "GRanges"
#     )
# )
# export.bed(chr_polyt, ignore.strand = TRUE,con=file.path(output_path, "chr_polyt.bed"))

# we get exon by transcript GrangesList object
ebt = suppressMessages({
    suppressWarnings({
      grl <- eisaR::getFeatureRanges(
        gtf = file.path(gtf_path),
        featureType = c("spliced"),
        verbose = FALSE
    )})
  })

export.bed(unlist(ebt), con=file.path(output_path, "exon_by_tx.bed"))

t2g = getTx2Gene(ebt)

# single-exon transcripts

batch_size=ceiling(length(ebt)/threads)

num_exon = foreach (batch_id=seq(1,length(ebt),by=batch_size), .combine=c, .inorder=TRUE) %dopar% {
    ebt_batch = ebt[batch_id:min(batch_id+batch_size-1, length(ebt))]
    # check length
    sapply(ebt_batch,length)
}

# Then, use the range as the unspliced transcript
unspliced_gr = unlist(range(ebt[num_exon > 1]))
# We want to remove the single exon transcripts from unspliced list
unspliced_gr$transcript_id = paste0(names(unspliced_gr), "-U")
unspliced_gr$gene_id = t2g[names(unspliced_gr), "gene_id"]
unspliced_gr$exon_id <- unspliced_gr$transcript_id
unspliced_gr$exon_name <- NULL
unspliced_gr$exon_rank <- as.integer(1)
unspliced_gr$type <- "exon"
names(unspliced_gr) <- NULL
mcols(unspliced_gr) <- S4Vectors::mcols(unspliced_gr)[, c("exon_id", "exon_rank", 
                                        "transcript_id", "gene_id", "type")]

unspliced_grl <- split(unspliced_gr, unspliced_gr$transcript_id)

# the spliceu is the combination of the unspliced and spliced
grlfull <- c(grlfull, ebt)
grlfull <- c(grlfull, unspliced_grl)

exportToGtf(grlfull, file.path(output_path, "spliceu.gtf"))

t2g = getTx2Gene(grlfull)
t2g$splicing_status = "S"
t2g$splicing_status[endsWith(t2g$transcript_id, "-U")] = "U"

utils::write.table(t2g,
    file = file.path(output_path, "t2g_3col.tsv"), sep = "\t",
    row.names = FALSE, quote = FALSE,
    col.names = FALSE, append = TRUE)

# Next we find the polyA and polyT sites

# get the sequence of transcripts
seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = genome,
    transcripts = grlfull
)
names(seqs) <- names(grlfull)

Biostrings::writeXStringSet(seqs, file.path(output_path, "spliceu.fa"), format = "fasta")

batch_size=ceiling(length(seqs)/threads)
shuf_ids = sample(1:length(seqs))

# polyat = for(batch_id in seq(1,length(seqs),by=batch_size)) {
polyat = foreach (batch_id=seq(1,length(seqs),by=batch_size), .combine=c, .inorder=TRUE) %dopar% {
    # obtain the batch
    batch = seqs[shuf_ids[batch_id:min(batch_id+batch_size-1, length(seqs))]]

    # get poly-A count
    # This returns a zero based index
    # This is zero based inclusive index
    polya = as(
        Biostrings::vmatchPattern(
        paste0(
            rep("A", min_length),
            collapse = ""
        ),
        batch,
        max.mismatch = max_mismatch
        ),
        "GRanges"
    )

    # filter outboundary records
    seqlevels(polya) = names(batch)
    seqlengths(polya) = width(batch)
    polya = trim(polya)
    polya = polya[width(polya) == min_length]
    # get polya sequence
    polya_seqs = getSeq(batch, polya)
    # keep only the polyA sites that start with A and end with A
    startsWidthA = substring(polya_seqs, 1, 1) == "A"
    endsWithA = substring(polya_seqs, min_length, min_length) == "A"
    polya = polya[startsWidthA & endsWithA]
    polya = IRanges::reduce(polya)
    # polya = shift(polya, 1)
    strand(polya) = "+"

    # get poly-T count
    polyt =as(
            Biostrings::vmatchPattern(
            paste0(
                rep("T", min_length),
                collapse = ""
            ),
            batch,
            max.mismatch = max_mismatch
        ),
        "GRanges"
    )
    
    # filter outboundary records
    seqlevels(polyt) = names(batch)
    seqlengths(polyt) = width(batch)
    polyt = trim(polyt)
    polyt = polyt[width(polyt) == min_length]
    polyt_seqs = getSeq(batch, polyt)
    startsWidthT = substring(polyt_seqs, 1, 1) == "T"
    endsWithT = substring(polyt_seqs, min_length, min_length) == "T"
    polyt = polyt[startsWidthT & endsWithT]
    polyt = IRanges::reduce(polyt)

    strand(polyt) = "-"
    c(polya,polyt)
}

# Then we need to write the polyA and polyT sites to a bed file 
export.bed(polyat[strand(polyat) == "+"], ignore.strand = TRUE, con=file.path(output_path, "tx_polya.bed"))
export.bed(polyat[strand(polyat) == "-"], ignore.strand = TRUE,con=file.path(output_path, "tx_polyt.bed"))
