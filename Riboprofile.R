library(optparse)
library(dplyr)
library(GenomicRanges)
library(diffloop)
library(stringr)
library(GenomicAlignments)
library(BuenColors)
library(scales)
"%ni%" <- Negate("%in%")

#---------------------
# Command line options
#---------------------
option_list <- list(
  
  make_option(c("--rg"),  default="hg38", help="Reference genome for the associated file"),
  make_option(c("--outdir"), default="ribo-out",help="Output folder for requisite files"),
  make_option(c("--fastq"),  default="test_1.fastq.gz", help="Reference genome for the associated file"),
  make_option(c("--cores"),  default="4", help="Number of cores for analysis")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Make into more pretty variables
outdir = opt[["outdir"]]
rg = opt[["rg"]]
fastqfile = opt[["fastq"]]
cores = opt[["cores"]]

# Modify
fastq <- basename(fastqfile)

# Set up output files and verify paths
dir.create(file.path(outdir), showWarnings = FALSE)
stopifnot(fastq != "")
stopifnot(rg %in% c("hg38", "mm10"))
print(opt)

# Binary file paths
star <- "/apps/lab/aryee/STAR-2.6.0c/bin/Linux_x86_64_static/STAR"
rob <- "/data/aryee/caleb/ribosome-profiling/riboprofile/ribORF.0.1"
fastqc <- "/source/fastQC/0.11.2/fastqc"
cutadapt <- "/data/aryee/caleb/pythondev/venv3/bin/cutadapt"


# Parse reference genome
if(rg == "hg38"){
  g_base <- "/data/joung/caleb/base_editing/hg38"
  rRNA <- diffloop::bedToGRanges(paste0(g_base, "/hg38_rRNA.bed"))
  tRNA <- diffloop::bedToGRanges(paste0(g_base, "/hg38_tRNA.bed"))
  
  a_base <-  paste0(g_base, "/rds_anno/")
  exons <- readRDS(paste0(a_base, "exons.hg38.gr.rds"))
  introns <- readRDS(paste0(a_base, "introns.hg38.gr.rds"))
  x5utr <- readRDS(paste0(a_base, "x5utr.hg38.gr.rds"))
  x3utr <- readRDS(paste0(a_base, "x3utr.hg38.gr.rds"))
  
} else if(rg == "mm10"){
  g_base <- "/data/joung/caleb/base_editing/GRCm38"
  rRNA <- diffloop::bedToGRanges(paste0(g_base, "/mm10_rRNA.bed")) %>% rmchr()
  tRNA <- diffloop::bedToGRanges(paste0(g_base, "/mm10_tRNA.bed")) %>% rmchr()
  
  a_base <-  paste0(g_base, "/rds_anno/")
  exons <- readRDS(paste0(a_base, "exons.mm10.gr.rds")) %>% rmchr()
  introns <- readRDS(paste0(a_base, "introns.mm10.gr.rds")) %>% rmchr()
  x5utr <- readRDS(paste0(a_base, "x5utr.mm10.gr.rds")) %>% rmchr()
  x3utr <- readRDS(paste0(a_base, "x3utr.mm10.gr.rds")) %>% rmchr()
} else {
  
  stop("Reference genome not found!")
}


# Simple function to parse the STAR output and get the number of uniquely aligned reads
parse_star <- function(log){
  lines <- readLines(log)
  avg_mapped_length <- strsplit(lines[[11]], "\t")[[1]][2]
  reads_out  <- strsplit(lines[[9]], "\t")[[1]][2]
  reads_out
}

#-------------------
# 01 - Run cut adapt
#-------------------
fastq_trimmed <- paste0(outdir, "/", gsub(".fastq.gz", ".trim.fastq.gz", fastq))
ca_call <- paste0(cutadapt, " -u 3 -j ",cores," --report=minimal --discard-untrimmed -a AGATCGGAAGAGCACACGTCTG -q 5 -m 20 ", fastqfile," -o ", fastq_trimmed)
ca_out <- system(ca_call, intern = TRUE) %>% stringr::str_split_fixed("\t", 10)

#--------------
# 02 - Run STAR
#--------------
outname <- gsub(".fastq.gz", "", fastq)
star_call <- paste0(star, " --runMode alignReads --readFilesCommand zcat --outFilterMultimapNmax 1 --outFileNamePrefix ", outdir, "/", outname,
                    " --runThreadN ",cores," --genomeDir ", g_base, " --readFilesIn ", fastq_trimmed, " --outSAMtype BAM Unsorted")
star_out <- system(star_call, intern = TRUE)
star_qc_reads <- parse_star(paste0(outdir, "/", outname, "Log.final.out"))

#-----------------------------------
# 03 - Import uniquely mapping reads / find overlaps
#-----------------------------------
GA <- readGAlignments(paste0(outdir, "/", outname, "Aligned.out.bam"))
tRNA_ov <- findOverlaps(GA, tRNA)
rRNA_ov <- findOverlaps(GA, rRNA)
GA_filt <- GA[1:length(GA) %ni% c(queryHits(rRNA_ov), queryHits(tRNA_ov))]

# Look for overlaps with other annotations
ov_5 <- findOverlaps(GA_filt, x5utr)
ov_3 <- findOverlaps(GA_filt, x3utr)
ov_cds <- findOverlaps(GA_filt, c(exons, introns))

# Classify remaining reads
class <- ifelse(1:length(GA_filt) %in% queryHits(ov_5), "UTR5prime", 
                ifelse(1:length(GA_filt) %in% queryHits(ov_3), "UTR3prime",
                       ifelse(1:length(GA_filt) %in% queryHits(ov_cds), "CDS", "other")))

#-----------------------------------
# 04 - Make qc plots
#-----------------------------------

# Plot 1 -- Read allocation
ws <- c("all read", "w/ adapt", "unique align", "non-t/rRNA") 
qc1df <- data.frame(
  what = factor(ws, levels = rev(ws)),
  n = c(as.numeric(ca_out[2,2]), as.numeric(ca_out[2,7]), as.numeric(star_qc_reads), length(class))
) %>% mutate(prop = n / max(n)) %>%
  mutate(label = paste0("n=", prettyNum(n,big.mark=",",scientific=FALSE), 
                        " (", prettyNum(prop*100,big.mark=",",scientific=FALSE, digits=4),"%)")) %>%
  mutate(text_coord = pmin(prop*100 + 25, 50))

P1 <- ggplot(qc1df, aes(x = what, y = prop*100, label = label)) +
  geom_bar(stat = "identity", fill = "lightgrey", color = "black") + coord_flip() +
  pretty_plot() + L_border() +
  geom_text(size = 5, aes(y = text_coord)) +
  labs(y = "% of reads", x = "Read type") +
  scale_y_continuous(expand = c(0,0))

# Plot 2 -- Fragment plot
qc2df <- data.frame(
  frag_size = qwidth(GA)
) %>% group_by(frag_size) %>% summarize(count = n()) %>%
  mutate(Rbf = frag_size >= 27 & frag_size <= 32) %>%
  mutate(prop = count / sum(count)) 

red_quant <- qc2df %>% filter(Rbf) %>% pull(prop) %>% sum() 
red_quant2 <- paste0(prettyNum(red_quant*100, big.mark=",",scientific=FALSE, digits = 4), "%")
P2 <- ggplot(qc2df, aes(x = frag_size, y = prop*100, fill = Rbf, color = Rbf)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0,0)) +
  pretty_plot() + L_border() + theme(legend.position = "none") +
  scale_fill_manual(values = c("lightgrey", "firebrick")) +
  scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "Fragment size", y = "% of fragments") +
  annotate("text", label = red_quant2, x = 25, y = 10, size = 5, colour = "firebrick")


# Plot 3 -- Pie graph of CDS / UTR / etc
qc3df <- data.frame(
  class = factor(class, levels = rev(c("CDS",  "UTR3prime", "UTR5prime","other"))), stringsAsFactors = FALSE
) %>% group_by(class) %>% summarize(n = n()) %>%
  mutate(value = n / sum(n) * 100) %>%
  mutate(pos = cumsum(value)- value/2)

P3 <- ggplot(qc3df, aes(x="", y=value, fill=class, label = as.character(value)))+
  geom_bar(width = 0.7, stat = "identity", color = "black") + 
  theme(axis.text.x=element_blank()) +
  pretty_plot() + L_border() + scale_y_continuous(expand = c(0,0)) +
  labs (x = "Non t/rRNA fragments", fill = "Class", y = "% of fragments") +
  scale_fill_manual(values = c("grey", "firebrick", "green4", "dodgerblue3"))

cowplot::ggsave(cowplot::plot_grid(P1, P2, P3, nrow = 1,  labels = c("a", "b", "c")),
                file = paste0(outdir, "/", gsub(".fastq.gz", ".qc.pdf", fastq)), 
                width = 12, height = 4)


