library(visutils)
library(Seurat)
library(DropletUtils)
library(plyr)
library(readr)
library(BiocParallel)

theargs <- R.utils::commandArgs(asValues=TRUE)

samples_tsv <- theargs$samples
input_path <- theargs$input
output_path <- theargs$output

##Reading in manifest TSV file
samples = read_tsv(samples_tsv)
samples = samples[samples$directory != '', ]

##Getting list of summary csvs for each cellranger run in this patient
summary_csvs = Sys.glob(file.path(input_path, "*", "*", "summary.csv"))

##Creating dataframe containing all data from each cellranger summary csv
qcs = do.call(rbind,lapply(summary_csvs,function(f)read.csv(f,check.names=F)))
qccols = colnames(qcs)

##Adding joint sanger ID to dataframe
qcs$sanger_ids = samples$`joint sanger id`

##Checking the shared name between atac and gex sample names and adding that to the dataframe
atac <- as.data.frame(gsub("_mA", "", samples$ATAC_SAMPLE_NAME))
colnames(atac) <- "atac_shared_name"
gex <- as.data.frame(gsub("_mG", "", samples$GEX_SAMPLE_NAME))
colnames(gex) <- "gex_shared_name"
if (all(atac$atac_shared_name == gex$gex_shared_name) == TRUE) {
  qcs$sample_name = atac$atac_shared_name
}

##Filter samples with less than 2000 cells estimated
qcs[qcs[,"Estimated number of cells"]<2000,]

##Extracting information from sample name and adding to dataframe
sample.meta = do.call(rbind,strsplit(qcs$sample_name,'-'))
qcs$donor = sample.meta[,1]
qcs$layer = as.numeric(substr(sample.meta[,5],1,1))
qcs$xpos = substr(sample.meta[,5],2,2)
qcs$ypos = as.numeric(substr(sample.meta[,5],3,3))
qcs$replicate = substr(sample.meta[,5],nchar(sample.meta[,5]),nchar(sample.meta[,5]))

qcs = qcs[,c(setdiff(colnames(qcs),qccols),qccols)]
qcs[is.na(qcs$ypos),'sample_name']

##Import per barcode metrics
per_barcode_metrics_csvs = Sys.glob(file.path(input_path, "*", "*", "per_barcode_metrics.csv"))
bc_metrics = lapply(per_barcode_metrics_csvs,function(f){print(f);read.csv(f)})
names(bc_metrics) = qcs$sanger_ids

##Defining cell filtering function
cellFilter = function(bc,gex_umis_count=1000,gex_genes_count=500,atac_fragments=1000){
  bc$is_cell == 1 & 
    bc$gex_umis_count >= gex_umis_count & 
    bc$gex_genes_count >= gex_genes_count &
    bc$atac_fragments >= atac_fragments
}

##Adding filtered cell number to dataframe
qcs$`Estimated number of filtered cells` = sapply(bc_metrics,function(x)sum(cellFilter(x)))

##Generating list of cellranger raw matrices
raw_matrices = Sys.glob(file.path(input_path, "*", "*", "raw_feature_bc_matrix"))

##Calculating empty drops for each cellranger output
emptydrop = llply(1:nrow(as.data.frame(raw_matrices)),function(i){
  print(i)
  d = Read10X(raw_matrices[i])
  BPPARAM = MulticoreParam(workers = 16)
  ed.ge = emptyDrops(d$`Gene Expression`,BPPARAM = BPPARAM)
  ed.at = emptyDrops(d$`Peaks`,BPPARAM = BPPARAM)
  cmn = intersect(rownames(ed.ge),rownames(ed.at))
  
  r=data.frame(ge.fdr = ed.ge[cmn,'FDR'],at.fdr=ed.at[cmn,'FDR'])
  rm(d,ed.ge,ed.at)
  gc()
  r  
}, .parallel = FALSE, .progress = 'text')

dir.create(output_path ,showWarnings = FALSE)
saveRDS(emptydrop, paste0(output_path, '/emptyDrops.rds'))

##Add empty drops information to dataframe
edcell_counts = t(sapply(emptydrop,function(x)apply(x<0.001,2,sum,na.rm=TRUE)))
qcs$`Estimated number of cells emptyDrops_GEX` = edcell_counts[,1]
qcs$`Estimated number of cells emptyDrops_ATAC` = edcell_counts[,2]
qcs$`% cell lost` = (1 - qcs$`Estimated number of filtered cells`/qcs$`Estimated number of cells`)*100

saveRDS(qcs, paste0(output_path, '/qc.rds'))

##Save QC table to path
write_tsv(qcs , paste0(output_path, '/multiome.qc.tsv'))

##Defining statistic functions
getLayerStat = function(qc,donor,layer,metric,na2text=TRUE){
  qc = qc[qc$donor==donor & qc$layer==layer,]
  
  if(na2text){
    qc$xpos[is.na(qc$xpos)] = 'NA'
    qc$ypos[is.na(qc$ypos)] = 'NA'
  }
  qc$ypos = paste0(qc$ypos,qc$replicate)
  xs = sort(unique(qc$xpos))
  ys = sort(unique(qc$ypos))
  mat = array(NA,dim=c(length(xs),length(ys)),dimnames=list(xs,ys))
  for(i in seq_len(nrow(qc))){
    mat[qc$xpos[i],qc$ypos[i]] = qc[i,metric]
  }
  mat
}

plotLayerStat = function(qc,donor,layer,metric,na2text=TRUE,same.lims=TRUE,digits=NULL,ratio2log=40){
  d = getLayerStat(qc,donor,layer,metric,na2text)
  d = d[,rev(colnames(d)),drop=F]
  if(min(dim(d))==0){
    plot.new()
    return(NULL)
  }
  if(same.lims){
    zlim = range(qc[,metric],na.rm=T)
  }else
    zlim = range(d,na.rm=T)
  if(is.null(digits)){
    if(mean(d,na.rm=TRUE)>100)
      digits = 0
    else
      digits = 2
  }else{
    digits = 2
  }
  dt = round(d,digits=digits)
  if(zlim[1] > 0 & zlim[2]/zlim[1] > ratio2log){
    zlim = log1p(zlim)
    d = log1p(d)
  }
  imageWithText(d,dt,las=1,col=num2col(1:100,c('#EEEEEE','orange')),main=paste0('layer ',layer,'\n',metric),zlim=zlim)#,xaxlab=NULL)
  if(nrow(d)>1)
    for(i in 2:nrow(d))
      abline(v=i-0.5,lty=2)
  if(ncol(d)>1)
    for(i in 2:ncol(d)){
      p = colnames(d)[i-1]
      c = colnames(d)[i]
      if(substr(p,1,nchar(p)-1) != substr(c,1,nchar(c)-1) )
        abline(h=i-0.5,lty=2)
    }
  invisible(d)
}

##Specifying stats to look at
stats = c("Estimated number of cells",
          "% cell lost",
          "Estimated number of cells emptyDrops_GEX",
          "ATAC Median high-quality fragments per cell",
          "ATAC Mean raw read pairs per cell",
          "ATAC Percent duplicates",
          "ATAC TSS enrichment score",
          "GEX Median UMI counts per cell",
          "GEX Mean raw reads per cell",
          "GEX Percent duplicates")

##Plotting each defined stat for each run of each block
dir.create(paste0(output_path, '/figures'),showWarnings = FALSE)
pdf(paste0(output_path,'/figures/multiome.qc.pdf'), width=length(stats)*3, height=9)
for(donor in sort(unique(qcs$donor))){
  par(mfcol=c(4,length(stats)),bty='n',oma=c(0,0,1.5,0),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,2.5,0),cex.main=0.9)
  for(stat in stats){
    for(l in 1:4){
      plotLayerStat(qcs,donor,l,stat)
    }
  }
  mtext(donor,side=3,line=0,outer=TRUE)
}
dev.off()

##Defining characteristics for plot
col = char2col(qcs$donor,palette=TRUE)

##Defining plotting function
plotM = function(qc,stat1,stat2,...){
  plot(qc[,stat1],qc[,stat2],xlab=stat1,ylab=stat2,...)
}

##Plotting various scatter comparisons of defined stats
pdf(paste0(output_path,'/figures/multiome.qc.scatters.pdf'), width=9, height=12)
par(mfrow=c(4,3),bty='n',oma=c(0,0,1.5,0),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,2.5,0))
plotM(qcs,"Estimated number of cells","Estimated number of cells emptyDrops_GEX",pch=16,col=col[qcs$donor])
abline(a=0,b=1,lty=2)
plotM(qcs,"Estimated number of cells","Estimated number of cells emptyDrops_ATAC",pch=16,col=col[qcs$donor])
abline(a=0,b=1,lty=2)
plotM(qcs,"Estimated number of cells","% cell lost",pch=16,col=col[qcs$donor])

plotM(qcs,"Estimated number of cells","ATAC Median high-quality fragments per cell",pch=16,col=col[qcs$donor])
plotM(qcs,"Estimated number of cells","ATAC Percent duplicates",pch=16,col=col[qcs$donor])
plotM(qcs,"ATAC Median high-quality fragments per cell","ATAC Percent duplicates",pch=16,col=col[qcs$donor])

plotM(qcs,"Estimated number of cells","GEX Median UMI counts per cell",pch=16,col=col[qcs$donor])
plotM(qcs,"Estimated number of cells","GEX Percent duplicates",pch=16,col=col[qcs$donor])
plotM(qcs,"GEX Median UMI counts per cell","GEX Percent duplicates",pch=16,col=col[qcs$donor])

plotM(qcs,"ATAC Median high-quality fragments per cell","GEX Median UMI counts per cell",pch=16,col=col[qcs$donor])
plotM(qcs,"ATAC Percent duplicates","GEX Percent duplicates",pch=16,col=col[qcs$donor])

legend('bottomright',col=col,pch=16,legend=names(col),bty='n')
dev.off()

##Filter out samples with more than 50 cells lost
qcs[qcs$`% cell lost`>50,]

##Plotting cell counts
pdf(paste0(output_path,'/figures/multiome.cellcounts.scatters.pdf'), width=9, height=6)
par(mfrow=c(2,3),bty='n',oma=c(0,0,1.5,0),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,2.5,0))
plotM(qcs,"Estimated number of cells","Estimated number of filtered cells",pch=16,col=col[qcs$donor])
dev.off()

##Generating list of cellranger filtered matrices
filtered_matrices = Sys.glob(file.path("/nfs/team283/GBM_LEAP/phase_2_multiome_data/AT19", "*", "*", "filtered_feature_bc_matrix"))

##Generating counts for each cellranger output
counts = lapply(filtered_matrices, function(f){print(f);rowSums(Read10X(f)$`Gene Expression`)})
counts = do.call(cbind,counts)
saveRDS(counts, paste0(output_path, '/pbcounts.rds'))

##Plotting GEX vs ATAC
lim=c(0,5)
ybins = xbins = seq(lim[1],lim[2],length.out=50)

zfun = sqrt
pdf(paste0(output_path, '/figures/atac2gex.hm.sqrt.pdf'), width=7*3, height=4*3)
par(mfcol=c(4,7),bty='n',oma=c(0,0,1.5,0),tcl=-0.2,mgp=c(1.3,0.3,0),mar=c(3,3,2.5,0),cex.main=0.9)
for(i in seq_len(nrow(qcs))){
  print(i)
  bc = bc_metrics[[qcs$sanger_ids[i]]]
  f = bc$is_cell==0
  z=hist2D(log10(1+bc$atac_peak_region_fragments[f]),log10(1+bc$gex_umis_count[f]),xbins,ybins,trimZq=0.0,zfun=zfun,cols=c('#0000FF00','#0000FFAA'),legend=FALSE,xlab='ATAC peak region fragments',ylab='GEX umis count',main=c(qcs$sample_name[i],qcs$sanger_ids[i]),xlim=lim,ylim=lim)
  p=hist2D(log10(1+bc$atac_peak_region_fragments[!f]),log10(1+bc$gex_umis_count[!f]),xbins,ybins,trimZq=0.0,zfun=zfun,cols=c('#FF000000','#FF0000AA'),new=FALSE,legend=FALSE)
  legend('topleft',fill=c('red','blue'),legend=paste0(c('Cells','Non-cells'), '(', c(sum(bc$is_cell==1),sum(bc$is_cell==0)), ')'))
}
dev.off()
