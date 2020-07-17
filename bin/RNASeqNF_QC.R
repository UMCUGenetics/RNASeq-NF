# Author:   Job van Riet
# Date:     01-07-2020
# Function: Retrieves output QC of the RNASeq-NF workflow.

# Retrieve input parameters -----------------------------------------------

# Check if all required parameters are given.
if(is.null(params$input) | params$input == '') stop('Please specify the QC input folder (--input)')


# Load libraries ----------------------------------------------------------

suppressPackageStartupMessages(require(futile.logger, quietly = T))
suppressPackageStartupMessages(require(tibble, quietly = T))
suppressPackageStartupMessages(require(tidyr, quietly = T))
suppressPackageStartupMessages(require(ggplot2, quietly = T))
suppressPackageStartupMessages(require(plotly, quietly = T))

# Initialize lists to store all input data and plots.
data.QC <- list()
plots.QC <- list()

futile.logger::flog.info('Importing and plotting QC - RNASeq-NF.')

# Helper functions --------------------------------------------------------

themeQC <- theme(
    text = element_text(family = 'Helvetica'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'top',
    plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.line = element_line(colour="black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = 'grey20', linetype = 'dotted'),
    panel.grid.minor.y = element_line(colour = 'grey50', linetype = 'dotted')
)


# Overview - Input files --------------------------------------------------

files.fq <- base::list.files(params$input, pattern = 'fastqc.html$', full.names = T, recursive = T)

if(length(files.fq) > 0){

    futile.logger::flog.info('\t- FASTQC')

    # Import data.
    data.QC$fastq <- tibble::tibble(file = gsub('_fastqc.*', '', base::basename(files.fq))) %>%
        dplyr::mutate(Direction = ifelse(grepl('_R1_', file), 'R1', 'R2')) %>%
        dplyr::mutate(Sample = gsub('_L.*_.*', '', file))

    # Plot data.
    plots.QC$fastq <- ggplot(data.QC$fastq, aes(x = Sample, fill = Direction)) +
        geom_bar(color = 'black', position = position_dodge2()) +
        labs(x = 'Samples', y = 'Number of files') +
        scale_fill_manual(values = c('#E3072A', '#5FACC8'), guide = guide_legend(title = NULL, title.position = NULL, title.hjust = 0.5, ncol = 1, keywidth = 0.75, keyheight = 0.75)) +
        themeQC
}


# Overview - TrimGalore ---------------------------------------------------

readTrimGaloreFiles <- function(x){

    z <- suppressWarnings(readr::read_delim(x, delim = '\t', col_names = 'Command', trim_ws = T, progress = F, col_types = 'c'))
    z <- z[22:27,] %>%
        tidyr::separate(col = Command, into = c('type', 'value'), sep = ':') %>%
        dplyr::mutate(type = trimws(type), value = as.integer(gsub('[^0-9.-]', '', trimws(gsub(' \\(.*', '', value))))) %>%
        dplyr::mutate(value.rel = 1, type2 = ifelse(grepl('eads', type), 'Read-counts', 'Base-counts')) %>%
        dplyr::mutate(sample = base::gsub('.fastq.gz_trimming_report.txt', '', base::basename(x)))

    z[2:3,]$value.rel <- z[2:3,]$value / z[1,]$value
    z[5:6,]$value.rel <- z[5:6,]$value / z[4,]$value

    return(z)
}

files.TrimGalore <- base::list.files(params$input, pattern = '_trimming_report.txt$', full.names = T, recursive = T)

if(length(files.TrimGalore) > 0){

    futile.logger::flog.info('\t- TrimGalore')

    # Import data.
    data.QC$TrimGalore <- base::do.call(base::rbind, lapply(files.TrimGalore, readTrimGaloreFiles))

    # Sort by decreasing reads with adapters.
    orderSamples <- (data.QC$TrimGalore %>% dplyr::filter(type == 'Reads with adapters') %>% dplyr::arrange(value.rel))$sample
    data.QC$TrimGalore <- data.QC$TrimGalore %>% dplyr::mutate(sample = factor(sample, levels = orderSamples))

    plots.QC$TrimGalore <- ggplot(data.QC$TrimGalore, aes(x = reorder(type, -value), y = sample, fill = value.rel, label = paste0(round(value/1E6, 2), 'M'))) +
        geom_tile(colour = 'grey80', size = 0.5, na.rm = T) +
        facet_wrap(type2 ~ ., ncol = 2, scales = 'free_x')+
        labs(x = NULL, y = NULL) +
        scale_fill_gradient2(low = 'white', mid = 'white', high = 'salmon') +
        geom_text(color = 'grey20', fontface = 'bold', size  = 3) +
        themeQC

}


# Overview - SortMeRNA ----------------------------------------------------

readSortMeRNAFiles <- function(x){

    z <- suppressWarnings(readr::read_delim(x, delim = '\t', col_names = 'Command', trim_ws = T, progress = F, col_types = 'c'))
    z <- z %>% dplyr::filter(grepl('E-value threshold', Command)) %>%
        tidyr::separate(col = Command, into = c('type', 'value'), sep = ' = ') %>%
        dplyr::mutate(type = trimws(type), value = as.integer(trimws(gsub(' \\(.*', '', value)))) %>%
        dplyr::mutate(value.rel = (value / sum(value))) %>%
        dplyr::mutate(sample = base::gsub('_rRNA_report.txt', '', base::basename(x)))

    return(z)
}

files.sortMeRNA <- base::list.files(params$input, pattern = 'rRNA_report.txt$', full.names = T, recursive = T)

if(length(files.sortMeRNA) > 0){
    futile.logger::flog.info('\t- SortMeRNA')

    # Import data.
    data.QC$SortMeRNA <- base::do.call(base::rbind, lapply(files.sortMeRNA, readSortMeRNAFiles))

    # Plot data.
    plots.QC$SortMeRNA <- ggplot(data.QC$SortMeRNA %>% dplyr::filter(grepl('passing', type)), aes(x = reorder(sample, -value.rel), y = value.rel, label = value, fill = type)) +
        geom_bar(color = 'black', stat = 'identity', position = position_dodge2()) +
        scale_y_continuous(labels = scales::percent, breaks = c(0, .01, .025, .05, .1, .25, .5, 1)) +
        labs(x = 'Samples', y = 'Total # reads (Relative)\nRibosomal Databases') +
        scale_fill_manual(values = c('#E3072A', '#5FACC8'), guide = guide_legend(title = NULL, title.position = NULL, title.hjust = 0.5, ncol = 1, keywidth = 0.75, keyheight = 0.75)) +
        themeQC

}


# Overview - Flagstat Metrics ---------------------------------------------

readBAMStatFiles <- function(x){
    z <- readr::read_table(x, skip = 4, col_names = c('type', 'value'), progress = F, col_types = 'cc') %>%
        dplyr::mutate(value = base::as.integer(base::trimws(base::gsub('chrom:', '', value)))) %>%
        dplyr::mutate(type = base::trimws(base::gsub(':', '', type))) %>%
        dplyr::mutate(sample = gsub('_Aligned.*', '', base::basename(x)))
    return(z)
}

files.BAMStat <- base::list.files(params$input, pattern = 'bam_stat.txt$', full.names = T, recursive = T)

if(length(files.BAMStat) > 0){
    futile.logger::flog.info('\t- BAMStat')

    # Import data.
    data.QC$BAMStat <- base::do.call(base::rbind, lapply(files.BAMStat, readBAMStatFiles))

    # Perform centering.
    data.QC$BAMStat <- data.QC$BAMStat %>% dplyr::group_by(type) %>% dplyr::mutate(meanValue = mean(value)) %>% dplyr::ungroup()
    data.QC$BAMStat <- data.QC$BAMStat %>% dplyr::mutate(
        value.centered = value - meanValue,
        value.centered = value.centered / 1E6,
        value.centered = ifelse(is.finite(value.centered), value.centered, 0),
    )

    # Plot data.
    plots.QC$BAMStat <- ggplot(data.QC$BAMStat, aes(x = reorder(type, -value), y = sample, fill = value.centered, label = paste0(round(value/1E6, 2), 'M'))) +
        geom_tile(colour = 'grey80', size = 0.5, na.rm = T) +
        geom_text(color = 'grey20', fontface = 'bold', size  = 3) +
        labs(x = 'BAM Statistics (RSeQC)', y = 'Samples with # reads\n(mean-centered / per million (M))') +
        scale_fill_gradient2(low = '#D23428', mid = 'white', high = '#2C6049', midpoint = 0) +
        themeQC

}


# Overview - TIN ----------------------------------------------------------

files.TIN <- base::list.files(params$input, pattern = '.summary.txt$', full.names = T, recursive = T)

if(length(files.TIN) > 0){
    futile.logger::flog.info('\t- TIN')

    # Import data.
    data.QC$TIN <- base::do.call(base::rbind, base::lapply(files.TIN, function(x) readr::read_delim(x, delim = '\t', col_types = 'cddd')))

    # Plot data.
    plots.QC$TIN <- ggplot(data.QC$TIN, aes(x = reorder(Bam_file, -`TIN(mean)`), y = `TIN(mean)`, label = `TIN(mean)`)) +
        geom_bar(fill = '#2C6049', color = 'black', stat = 'identity', position = position_dodge2()) +
        labs(x = 'Samples', y = 'Transcript Integrity Number\n(Median)') +
        themeQC

}

# Overview - Complexity Curves --------------------------------------------

readComplexityFiles <- function(x){

    z <- readr::read_delim(x, delim = '\t', col_types = 'dddd')
    z$sample <- unique(base::gsub('\\.out.ccurve.txt', '', basename(x)))

    return(z)

}

files.complexity <- base::list.files(params$input, pattern = 'out.ccurve.txt$', full.names = T, recursive = T)

if(length(files.complexity) > 0){
    futile.logger::flog.info('\t- Complexity Curves')

    # Import data.
    data.QC$complexity <- base::do.call(base::rbind, lapply(files.complexity, readComplexityFiles))

    # Plot data.
    plots.QC$complexity <- ggplot(data.QC$complexity, aes(x = TOTAL_READS, y = EXPECTED_DISTINCT, color = sample)) +
        geom_line() +
        scale_color_brewer(palette = 'Dark2', guide = guide_legend(title = NULL, title.position = NULL, title.hjust = 0.5, ncol = 2, keywidth = 0.75, keyheight = 0.75)) +
        labs(x = 'Samples', y = 'Transcript Integrity Number\n(Median)') +
        themeQC
}

