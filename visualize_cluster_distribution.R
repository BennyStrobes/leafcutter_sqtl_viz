args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(plyr)
library(foreach)
library(gridExtra)
library(Hmisc)
library(reshape2)

#' Make sashimi-esque plot with ggplot2
#'
#' Shows only junction reads and can optionally show splicing variation across groups.
#'
#' @param y [samples] x [introns] matrix of counts
#' @param x [samples] vector of group membership. Can be numeric, factor or character.
#' @param exons_table An optional data frame specifying exons, with columns: chr start end strand gene_name. For hg19 see data/gencode19_exons.txt.gz
#' @param len Number of segments each curve is constructed from, controls smoothness.
#' @param length_transform Function controlling how true genomic space is mapped to the plot for improved visability. Use the identity function (i.e. function(g) g
#' @param main_title Title
#' @param snp_pos An optional list of SNP positions to mark
#' @param summary_func Function to summarize counts within groups: usually colMeans or colSums
#' @param legend_title Match this to summary_func
#' @import ggplot2
#' @import foreach
#' @importFrom gridExtra grid.arrange
#' @importFrom Hmisc bezier
#' @importFrom reshape2 melt
#' @export
make_differential_splicing_plot=function(y, x=numeric(nrow(y))+1, exons_table=NULL, len=500, length_transform=function(g) log(g+1), main_title=NA, snp_pos=NA, summary_func=colMeans, legend_title="Mean count",output_file="temp.pdf") {

  # convert colnames(y) into meta data
  intron_meta=get_intron_meta(colnames(y))

  new_theme_empty <- theme_bw(base_size = 16)
  new_theme_empty$line <- element_blank()
  new_theme_empty$rect <- element_blank()
  new_theme_empty$strip.text <- element_blank()
  new_theme_empty$axis.text <- element_blank()
  
  groups=sort(unique(x))
  
  max_log=.5*ceiling(2*log10( 1+max( unlist( foreach (tis=groups) %do% { intron_meta$counts=summary_func(y[ tis==x,,drop=F]) } ) ) ))

  breaks=if (max_log <= 2.5) seq(0,max_log,by=0.5) else seq(0,ceiling(max_log),by=1)
  limits=c(0.0,max_log)

  intron_meta$id=as.factor(1:nrow(intron_meta))
  temp=intron_meta[,c("id","start","end")]
  m=reshape2::melt(temp, id.vars = "id")

  s=unique(m$value)
  if (!is.na(snp_pos)) s=c(s,snp_pos)
  s=sort(s)
  d=s[2:length(s)]-s[1:length(s)-1]
  trans_d=length_transform(d) # e.g. log(d+1), sqrt(d), d ^ .33
  coords=c(0,cumsum(trans_d)) 
  names(coords)=s

  snp_coord=coords[as.character(snp_pos)]

  total_length=sum(trans_d) # ==max(coords)
  my_xlim=c(-.05*total_length,1.05*total_length)
  first_plot=T
  
  min_height=0
  max_height=0
  
  plots=foreach (tis=groups) %do% {
    intron_meta$counts=summary_func(y[ tis==x,,drop=F])
    group_sample_size=sum(tis==x)
    
    allEdges=do.call(rbind,foreach (i=1:nrow(intron_meta)) %do% {
      if (intron_meta$counts[i]==0) return(NULL)
      start=coords[ as.character(intron_meta$start[i]) ]
      end=coords[ as.character(intron_meta$end[i]) ]
      l=end-start
      #edge = data.frame(bezier(x=c(start, (start + end)/2, end), y=c(0, (1+sqrt(l)) * ( (i %% 2)*2-1 ), 0),evaluation = len))
      h=(1+sqrt(l)) * ( (i %% 2)*2-1 )
      min_height=min(min_height,h)
      max_height=max(max_height,h)
      edge = data.frame(Hmisc::bezier(x=c(start, start-0.1*total_length, (start + end)/2, end+0.1*total_length, end), y=c(0, h*2/3, h, h*2/3, 0),evaluation = len))
      edge$Sequence <- log10(1+intron_meta$counts[i]) * sin( seq(0,pi,length.out=len) ) # For size and colour weighting in plot
      edge$log10counts=log10(intron_meta$counts[i])
      edge$Group <- i
      edge
    })

    if (is.na(main_title) | !first_plot) new_theme_empty$plot.title <- element_blank()
    first_plot=F
    
    g=ggplot(allEdges) + geom_path(aes(x = x, y = y, group = Group, colour=log10counts, size = Sequence, alpha=.9)) + scale_size(breaks=breaks, labels=format(10^breaks,digits=0), limits=limits, range = c(.3, 10), guide = guide_legend(title = legend_title)) + scale_alpha(guide="none",range = c(0.1, 1)) + new_theme_empty + scale_color_gradient(breaks=breaks, limits=limits, labels=format(10^breaks,digits=0),low="blue",high="red", guide = guide_legend(title = legend_title)) + ylab(paste0(tis," (n=",group_sample_size,")")) + xlab("") + xlim(my_xlim) + geom_hline(yintercept=0,alpha=.3) #  + scale_color_discrete(guide="none")
    if (!is.na(snp_coord)) {
        df=data.frame(x=snp_coord,xend=snp_coord,y=0,yend=max_height*1.1)
        g=g + geom_segment(data=df,aes(x=x,y=y,xend=xend,yend=yend)) #+geom_vline(xintercept=snp_coord)
    }
    g
  }
  
  df=data.frame(x=coords, xend=total_length*(s-min(s))/(max(s)-min(s)), y=0, yend=min_height)
  plots[[length(plots)]]=plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend),alpha=.1)
  
  if (!is.null(exons_table)) {
    exons_chr=exons_table[exons_table==intron_meta$chr[1],]
    exons_here=exons_chr[ ( min(s) <= exons_chr$start & exons_chr$start <= max(s) ) | ( min(s) <= exons_chr$end & exons_chr$end <= max(s) ), ]
    
    exons_here$gene_name=factor(exons_here$gene_name)
  
    gene_heights=min_height - ((1:length(levels(exons_here$gene_name)))-1.0) * abs(min_height) *  .15 # 0.05
    heights=gene_heights[ as.numeric(exons_here$gene_name)] # .15
    df=data.frame(x=total_length*(exons_here$start-min(s))/(max(s)-min(s)), xend=total_length*(exons_here$end-min(s))/(max(s)-min(s)), y=heights, yend=heights)
    if (nrow(exons_here)>0)
    #plots[[length(plots)]]=plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend), alpha=.3, size=5)  + geom_hline(yintercept=min_height,alpha=.3) + geom_text(data=data.frame(x=my_xlim[1], y=gene_heights, label=levels(exons_here$gene_name)),aes(x,y,label=label))
    plots[[length(plots)]]=plots[[length(plots)]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend), alpha=.3, size=5)  + geom_hline(yintercept=min_height,alpha=.3) #+ geom_text(data=data.frame(x=my_xlim[1], y=gene_heights, label=levels(exons_here$gene_name)),aes(x,y,label=label))

    invert_mapping=function(pos) if (pos %in% s) coords[as.character(pos)] else 
      if (pos < min(s)) my_xlim[1] else 
      if (pos > max(s)) my_xlim[2] else {
      w=which( pos < s[2:length(s)] & pos > s[1:(length(s)-1)] )
      stopifnot(length(w)==1)
      coords[w] + (coords[w+1]-coords[w])*(pos - s[w])/(s[w+1]-s[w])
      }
    if (nrow(exons_here)>0) {
      df=data.frame(x=sapply(exons_here$start,invert_mapping), xend=sapply(exons_here$end,invert_mapping), y=0, yend=0)
      for (i in 1:length(plots)) plots[[i]]=plots[[i]] + geom_segment(data=df, aes(x=x,y=y,xend=xend,yend=yend), alpha=.3, size=5)
    }
  }

  if (!is.na(main_title)) plots[[1]] = plots[[1]] + ggtitle(main_title)
  pdf(output_file, width=8, height=9)
  do.call( gridExtra::grid.arrange, c(plots, list(ncol=1)))
  dev.off()
  #if (!is.null(exons_table))
    # exons_here
}


#' Make a data.frame of meta data about the introns
#' @param introns Names of the introns
#' @return Data.frame with chr, start, end, cluster id and "middle" 
#' @export
get_intron_meta=function(introns){
  intron_meta=do.call(rbind,strsplit(introns,":"))
  colnames(intron_meta)=c("chr","start","end","clu")
  intron_meta=as.data.frame(intron_meta,stringsAsFactors = F)
  intron_meta$start=as.numeric(intron_meta$start)
  intron_meta$end=as.numeric(intron_meta$end)
  intron_meta$middle=.5*(intron_meta$start+intron_meta$end) # TODO: can remove this now? 
  intron_meta
}

load_data <- function(file_name) {
	data <- read.table(file_name, header=F, stringsAsFactors = F)
	col_string <- as.character(data[1,2:dim(data)[2]])
	data_short <- data[2:dim(data)[1],2:dim(data)[2]]
	# Conver to numeric
	data_short <- sapply(data_short, as.numeric)
	colnames(data_short) <- col_string
	return(data_short)
}





#################################
# Command Line args
#################################
clusters_to_plot_file <- args[1]  # File where each line is a cluster
tissue_name <- args[2]  # Name of tissue
exon_file <- args[3]  # File containing exon positions
visualize_cluster_distribution_dir <- args[4] 


tryCatch(
{
#visualize_cluster_distribution_dir <- "/Users/bennystrobes/Google Drive/Research/gtex_v8/rare_var/visualize/leaf_viz/visualize_cluster_distribution/"

#clusters_to_plot_file <- paste0(visualize_cluster_distribution_dir, clusters_to_plot_file)
#exon_file <- paste0(visualize_cluster_distribution_dir, exon_file)

exons_table <- read.table(exon_file,header=T)

clusters_to_plot <- read.table(clusters_to_plot_file, header=TRUE)



#cluster_assignments <- read.table(paste0(visualize_cluster_distribution_dir,cluster_assignment_file), header=TRUE)

#temper = cluster_assignments[as.character(cluster_assignments$cluster_id) == "cluster10259",]

#print(temper)
#print(temper[1,1])
#print(temper[1,2])
#print(temper[1,3])
#print(as.character(temper[1,4]))
num_clusters <- dim(clusters_to_plot)[1]

if (num_clusters != 0) {
for (cluster_num in 1:num_clusters) {
  cluster_id <- as.character(clusters_to_plot[cluster_num,1])
  ensamble_id <- as.character(clusters_to_plot[cluster_num,2])
  chrom_num <- as.character(clusters_to_plot[cluster_num,3])
	var_pos <- as.numeric(clusters_to_plot[cluster_num,4])
	#dist_to_ss <- as.character(clusters_to_plot[cluster_num,6])
	strand <- as.character(clusters_to_plot[cluster_num,8])
  #annotated_ss <- as.character(clusters_to_plot[cluster_num,8])
  #variant_allele <- as.character(clusters_to_plot[cluster_num,9])
	counts_file <- paste0(visualize_cluster_distribution_dir, tissue_name, "_", cluster_id, "_counts.txt")
	counts <- load_data(counts_file)
  #cluster_assignment_vec <- cluster_assignments[as.character(cluster_assignments$cluster_id) == cluster_id,]
  #es_boolean <- as.character(cluster_assignment_vec[1,2])
  #alt_5_boolean <- as.character(cluster_assignment_vec[1,3])
  #alt_3_boolean <- as.character(cluster_assignment_vec[1,4])
  geno_0_string = as.character(clusters_to_plot[cluster_num,5])
  geno_1_string = as.character(clusters_to_plot[cluster_num,6])
  geno_2_string = as.character(clusters_to_plot[cluster_num,7])

  num_geno_0 <- length(strsplit(geno_0_string,split=",")[[1]])
  num_geno_1 <- length(strsplit(geno_1_string,split=",")[[1]])
  num_geno_2 <- length(strsplit(geno_2_string,split=",")[[1]])


  membership_vector <- factor(c(rep("genotype=0",num_geno_0), rep("genotype=1",num_geno_1), rep("genotype=2",num_geno_2)))
	output_file <- paste0(visualize_cluster_distribution_dir, tissue_name, "_", cluster_id, "_", ensamble_id, "_cluster_viz.pdf")
	make_differential_splicing_plot(counts, x=membership_vector,main_title=paste0(cluster_id, " / ", ensamble_id," / strand=",strand,"\n", chrom_num,':',var_pos, "\n"), snp_pos=var_pos,exons_table=exons_table, output_file=output_file)
}
}

},
error = function(e){
  print("error")
}
)
