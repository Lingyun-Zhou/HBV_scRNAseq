#!/usr/bin/env Rscript
suppressMessages(library(optparse))
version <- '
     VERSION: 1.0
      AUTHOR: Fang Peng (fun_peng@foxmail.com)
Updated date: 2024-06
'

option_list<-list(
	make_option("--input_file", help="cnr file from cnvkit, required", dest='input_file'),
	make_option('--bed', help='Genome BED of HBV CDS region',dest='bed'),
	make_option(c('-o', '--outfile'), default='result.pdf', help="Output file with png format, default: [%default]", dest='outfile'),
	make_option(c("-v", "--version"), action="store_true", default=FALSE, help="Show soft version", dest='version')
)

parse <- OptionParser(option_list = option_list, usage = " %prog --input_file [require] --bed [require] --o [result.pdf]")
opt <- parse_args(parse)

if(opt$version){cat(version); q(status=0)}
if(is.null(opt$input_file) || !file.exists(opt$input_file)) { parse_args(parse, args=c('-h'))}
if(is.null(opt$bed) || !file.exists(opt$bed)) { parse_args(parse, args=c('-h'))}
if(is.null(opt$outfile)) { opt$outfile <- 'result.pdf'}

input_file = opt$input_file
bed = opt$bed
outfile = opt$outfile

# main
suppressWarnings(suppressMessages(library(circlize)))
options(stringsAsFactors = FALSE)

# annotation 
anno = read.table(bed, sep='\t')
anno <- anno[,c('V2', 'V3', 'V6', 'V9')]
colnames(anno) <- c('start', 'end', 'name', 'detail')
anno$color <- NA
anno[anno$name == 'P',]$color <- 'brown'
anno[anno$name == 'X',]$color <- 'deepskyblue'
anno[anno$name == 'C',]$color <- 'purple'
anno[anno$name == 'S',]$color <- 'green'

# input
data <- read.table(input_file, sep ='\t')
colnames(data) <- c('sample', 'readname', 'chr', 'start', 'end', 'seq', 'left_clip', 'right_clip', 'cell_barcode', 'umi', 'adjust_cell_barcode', 'celltype')
# samples = c("hcc_1", "hcc_3", "hcc_4", "hcc_5", "hcc_7", "hcc_9", "hcc_11", "hcc_12", "hcc_14", "hcc_15")

samples = unique(data$sample)

# sample count
plot_data <- data.frame()
for (chr in unique(data$chr)){
    data_sub <- data[data$chr == chr, ]
    plot_data <- rbind(plot_data, data.frame(
                        chr = chr, 
                        start =  min(c(data_sub$start, data_sub$end)):max(c(data_sub$start, data_sub$end)), 
                        end =  min(c(data_sub$start, data_sub$end)):max(c(data_sub$start, data_sub$end))
                    )
                )
}

for(sample in samples){
    data_sub <- data[data$sample == sample, ]
    plot_data[, sample] <- 0
    for(i in 1:nrow(data_sub)){
        plot_data[plot_data$chr == data_sub[i,]$chr & plot_data$start %in% data_sub[i,]$start:data_sub[i,]$end, sample] = plot_data[plot_data$chr == data_sub[i,]$chr & plot_data$start %in% data_sub[i,]$start:data_sub[i,]$end, sample] + 1
    }
}

# annotation plot
cairo_pdf(file = outfile, width = 700, height = 700)
circos.clear()
circos.par("track.height" = 0.1, "gap.degree" = 0, cell.padding = c(0, 0, 0, 0))
circos.initialize(sectors=plot_data$chr, x = plot_data$start)
circos.track(ylim = c(0, 1), track.height = 0.2, panel.fun = function(x, y) {
    y = 0
    for(i in 1:nrow(anno)){
        if(i > 1){
            if(anno[i,]$detail != anno[i-1,]$detail){
                y = y + 0.1
            }
            if(anno[i,]$name != anno[i-1,]$name){
                circos.text(x = anno[i-1,]$start + (anno[i-1,]$end - anno[i-1,]$start) / 5*3, y = y+0.1, labels = anno[i-1,]$name, cex = 1.5, col = anno[i-1,]$color)
            }
        }else{
            y = y + 0.1
        }
        if(i == nrow(anno)){
            circos.text(x = anno[i,]$start + (anno[i,]$end - anno[i,]$start) / 5 *2, y = y + 0.2, labels = anno[i,]$name, cex = 1.5, col = anno[i,]$color)
        }
        circos.segments(anno[i,]$start, y, anno[i,]$end, y, lwd = 5, col = anno[i,]$color)
    }
}, bg.lwd = 0.0001, bg.border = 'white')

for(sample in samples){
    circos.track(ylim = c(0, 1), track.height = 0.05, panel.fun = function(x, y) {
            circos.barplot(value = plot_data[,sample] / max(plot_data[,sample]), pos = plot_data$start, col = 'brown')
    }, bg.lwd = 0.0001, bg.border = 'white')
}

dev.off()
