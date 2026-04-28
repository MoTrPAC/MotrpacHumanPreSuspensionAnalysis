pre_cawg_get_peak_annotations_hs = function(counts_dt, species = "Homo Sapiens", release = 105, txdb = NULL) {
  if (!"feature_id" %in% colnames(counts_dt) & !data.table::is.data.table(counts_dt)) {
    genomic_peaks = data.table::data.table(
      feature_id = rownames(counts_dt),
      chrom = gsub(":.*", "", rownames(counts_dt)),
      start = as.numeric(gsub(".*:|-.*", "", rownames(counts_dt))),
      end = as.numeric(gsub(".*-", "", rownames(counts_dt)))
    )
  } else if (!"feature_id" %in% colnames(counts_dt) & data.table::is.data.table(counts_dt)) {
    counts = counts_dt
    counts[, feature_id := paste0(chrom, ':', start, '-', end)]
    genomic_peaks = counts[, .(chrom, start, end, feature_id)]
  } else if ("feature_id" %in% colnames(counts_dt) & data.table::is.data.table(counts_dt)) {
    counts = counts_dt
    genomic_peaks = counts[, .(chrom, start, end, feature_id)]
  } else {
    counts = data.table::as.data.table(counts_dt)
    genomic_peaks = counts[, .(chrom, start, end, feature_id)]
  }

  if (is.null(txdb)) {
    txdb = txdbmaker::makeTxDbFromEnsembl(organism = species, release = release)
  }

  accepted_chrom = GenomeInfoDb::seqlevels(txdb)
  accepted_chrom = accepted_chrom[!grepl("\\.", accepted_chrom)]

  genomic_peaks = genomic_peaks[!grepl("\\.", chrom)]
  genomic_peaks[, chrom := gsub("^chr", "", as.character(chrom))]

  if (!all(unique(genomic_peaks[, chrom]) %in% accepted_chrom)) {
    stop(sprintf(
      "The following chromosomes are found in the input but not in the txdb object: %s",
      paste0(unique(!genomic_peaks[, chrom] %in% accepted_chrom), collapse = ', ')
    ))
  }

  peak = GenomicRanges::GRanges(
    seqnames = genomic_peaks[, chrom],
    ranges = IRanges::IRanges(as.numeric(genomic_peaks[, start]), as.numeric(genomic_peaks[, end]))
  )
  peakAnno = ChIPseeker::annotatePeak(peak,
                                      level = "gene",
                                      tssRegion = c(-2000, 1000),
                                      TxDb = txdb,
                                      overlap = "all")
  pa = data.table::as.data.table(peakAnno@anno)

  if (nrow(pa) == nrow(genomic_peaks)) {
    pa[, feature_id := genomic_peaks[, feature_id]]
  } else {
    cols = c('seqnames', 'start', 'end')
    pa[, (cols) := lapply(.SD, as.character), .SDcols = cols]
    cols = c('chrom', 'start', 'end')
    genomic_peaks[, (cols) := lapply(.SD, as.character), .SDcols = cols]
    pa = merge(pa, genomic_peaks, by.x = c('seqnames', 'start', 'end'), by.y = c('chrom', 'start', 'end'), all.y = TRUE)
  }

  pa[, short_annotation := annotation]
  pa[grepl('Exon', short_annotation), short_annotation := 'Exon']
  pa[grepl('Intron', short_annotation), short_annotation := 'Intron']

  pa[, c('geneChr', 'strand') := NULL]

  cols = c('start', 'end', 'geneStart', 'geneEnd', 'geneStrand')
  pa[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
  pa[, dist_upstream := ifelse(end - geneStart <= 0, end - geneStart, NA_real_)]
  pa[, dist_downstream := ifelse(start - geneEnd >= 0, start - geneEnd, NA_real_)]
  pa[end >= geneStart & start <= geneEnd, dist_downstream := 0]
  pa[end >= geneStart & start <= geneEnd, dist_upstream := 0]
  pa[, relationship_to_gene := ifelse(is.na(dist_downstream), dist_upstream, dist_downstream)]
  pa[, c('dist_upstream', 'dist_downstream') := NULL]

  pa[relationship_to_gene == 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Overlaps Gene"]
  pa[geneStrand == 1 & relationship_to_gene > 0 & relationship_to_gene < 5000, short_annotation := "Downstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 0 & relationship_to_gene > -5000, short_annotation := "Downstream (<5kb)"]
  pa[geneStrand == 1 & relationship_to_gene > -5000 & relationship_to_gene < 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[geneStrand == 2 & relationship_to_gene < 5000 & relationship_to_gene > 0 & grepl("Downstream|Intergenic", short_annotation), short_annotation := "Upstream (<5kb)"]
  pa[abs(relationship_to_gene) >= 5000, short_annotation := "Distal Intergenic"]

  data.table::setnames(pa,
    c('short_annotation', 'annotation', 'seqnames', 'geneId'),
    c('custom_annotation', 'chipseeker_annotation', 'chrom', 'ensembl_gene')
  )

  return(pa)
}
