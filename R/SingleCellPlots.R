#' Plot cell annotation
#' @param  object The SingCellaR object.
#' @param  type Type of plots including the histogram and boxplot.
#' @export
#' 
plot_cells_annotation <- function(object,type=c("histogram","boxplot")){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	cells.anno<-get_cells_annotation(object)
	##plot number of UMI count per cell##
	type <- match.arg(type)
	
	if(type=="histogram"){
		par(mfrow=c(2,2))
		hist(cells.anno$UMI_count,breaks=100,xlim=c(min(cells.anno$UMI_count),max(cells.anno$UMI_count)),col="cyan",xlab="UMI count",ylab="Cell frequency",main="Number of UMI count per cell")
		abline(v=mean(cells.anno$UMI_count),lwd=2,lty=2,col="red")
		text1<-paste("mean =",format(mean(cells.anno$UMI_count),digit=2),sept="")
		
		legend("topright",text1, col = 2:3, lty = 2:3, lwd = 2)
		
		##plot number of detected genes per cell##
		hist(cells.anno$detectedGenesPerCell,breaks=100,xlim=c(min(cells.anno$detectedGenesPerCell),max(cells.anno$detectedGenesPerCell)),col="blue",xlab="Number of detected genes",ylab="Cell frequency",main="Number of detected genes per cell")
		abline(v=mean(cells.anno$detectedGenesPerCell),lwd=2,lty=2,col="red")
		
		text2<-paste("mean =",format(mean(cells.anno$detectedGenesPerCell),digit=2),sept="")
		legend("topright",text2, col = 2:3, lty = 2:3, lwd = 2)
		
		##plot the percentage of mitocondrial genes##
		hist(cells.anno$percent_mito,breaks=100,xlim=c(min(cells.anno$percent_mito),max(cells.anno$percent_mito)),col="azure",xlab="% of mitocondrial genes",ylab="Cell frequency",main="% of mitocondrial genes per cell")
		abline(v=mean(cells.anno$percent_mito),lwd=2,lty=2,col="red")
		
		text3<-paste("mean =",format(mean(cells.anno$percent_mito),digit=2),sept="")
		legend("topright",text3, col = 2:3, lty = 2:3, lwd = 2)
		
	}else if(type=="boxplot"){
		
		#UMI
		umi.info<-cells.anno["UMI_count"]
		umi.info$Cell<-"cell"
		
		#detectedGenesPerCell
		gene.info<-cells.anno["detectedGenesPerCell"]
		gene.info$Cell<-"cell"
		
		#percent_mito
		mito.info<-cells.anno["percent_mito"]
		mito.info$Cell<-"cell"
		
		my.p<-list()
		my.p[[1]]<-ggplot(umi.info, aes(Cell, UMI_count)) +
				geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
				geom_jitter(aes(Cell, UMI_count),
						position=position_jitter(width=0.3,height=0),
						alpha=0.1,col="cyan",
						size=2)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
						panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x=element_blank())
		
		my.p[[2]]<-ggplot(gene.info, aes(Cell, detectedGenesPerCell)) +
				geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
				geom_jitter(aes(Cell, detectedGenesPerCell),
						position=position_jitter(width=0.3,height=0),
						alpha=0.1,col="orange",
						size=2)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
						panel.background = element_blank(), axis.line = element_line(colour = "black"),
						axis.text.x=element_blank())
		
		my.p[[3]]<-ggplot(mito.info, aes(Cell, percent_mito)) +
				geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
		    scale_y_continuous(breaks=c(0,5,10,15,20,25,30,40,50,60,70,80,90,100))+
				geom_jitter(aes(Cell, percent_mito),
						position=position_jitter(width=0.3,height=0),
						alpha=0.05,col="purple",
						size=2)
		
		grid.arrange(grobs = my.p, ncol = 2, as.table = FALSE)
		
	}else{
		stop("please select the type of plots!")
	}
}

#' Plot cell annotation for Target-Seq data
#' @param  object The SingCellaR object.
#' @param  type Type of plots including the histogram and boxplot.
#' @export

TargetSeq_plot_cells_annotation <- function(object,type=c("histogram","boxplot")){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  cells.anno<-get_cells_annotation(object)
  ##plot number of UMI count per cell##
  type <- match.arg(type)
  
  if(type=="histogram"){
    par(mfrow=c(2,2))
    hist(cells.anno$UMI_count,breaks=100,xlim=c(min(cells.anno$UMI_count),max(cells.anno$UMI_count)),col="cyan",xlab="Read count",ylab="Cell frequency",main="Number of reads per cell")
    abline(v=mean(cells.anno$UMI_count),lwd=2,lty=2,col="red")
    text1<-paste("mean =",format(mean(cells.anno$UMI_count),digit=2),sept="")
    
    legend("topright",text1, col = 2:3, lty = 2:3, lwd = 2)
    
    ##plot number of detected genes per cell##
    hist(cells.anno$detectedGenesPerCell,breaks=100,xlim=c(min(cells.anno$detectedGenesPerCell),max(cells.anno$detectedGenesPerCell)),col="blue",xlab="Number of detected genes",ylab="Cell frequency",main="Number of detected genes per cell")
    abline(v=mean(cells.anno$detectedGenesPerCell),lwd=2,lty=2,col="red")
    
    text2<-paste("mean =",format(mean(cells.anno$detectedGenesPerCell),digit=2),sept="")
    legend("topright",text2, col = 2:3, lty = 2:3, lwd = 2)
    
    ##plot the percentage of mitocondrial genes##
    hist(cells.anno$percent_mito,breaks=100,xlim=c(min(cells.anno$percent_mito,na.rm=T),max(cells.anno$percent_mito,na.rm=T)),col="azure",xlab="% of mitocondrial genes",ylab="Cell frequency",main="% of mitocondrial genes per cell")
    abline(v=mean(cells.anno$percent_mito),lwd=2,lty=2,col="red")
    
    text3<-paste("mean =",format(mean(cells.anno$percent_mito,na.rm=T),digit=2),sept="")
    abline(v=mean(cells.anno$percent_mito,na.rm=T),lwd=2,lty=2,col="red")
    legend("topright",text3, col = 2:3, lty = 2:3, lwd = 2)
    
    ##plot the percentage of ERCC genes##
    hist(cells.anno$percent_ERCC,breaks=100,xlim=c(min(cells.anno$percent_ERCC,na.rm=T),
                                                   max(cells.anno$percent_ERCC,na.rm=T)),
         col="azure",xlab="% of ERCC",ylab="Cell frequency",main="% of ERCC per cell")
    abline(v=mean(cells.anno$percent_ERCC),lwd=2,lty=2,col="red")
    
    text4<-paste("mean =",format(mean(cells.anno$percent_ERCC,na.rm=T),digit=2),sept="")
    abline(v=mean(cells.anno$percent_ERCC,na.rm=T),lwd=2,lty=2,col="red")
    legend("topright",text3, col = 2:3, lty = 2:3, lwd = 2)
    
  }else if(type=="boxplot"){
    
    #UMI
    umi.info<-cells.anno["UMI_count"]
    umi.info$Cell<-"cell"
    colnames(umi.info)<-c("Read_count","Cell")
    
    #detectedGenesPerCell
    gene.info<-cells.anno["detectedGenesPerCell"]
    gene.info$Cell<-"cell"
    
    #percent_mito
    mito.info<-cells.anno["percent_mito"]
    mito.info$Cell<-"cell"
    
    #percent_ERCC
    ERCC.info<-cells.anno["percent_ERCC"]
    ERCC.info$Cell<-"cell"
    
    my.p<-list()
    my.p[[1]]<-ggplot(umi.info, aes(Cell, Read_count)) +
      geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
      geom_jitter(aes(Cell, Read_count),
                  position=position_jitter(width=0.3,height=0),
                  alpha=0.1,col="cyan",
                  size=2)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x=element_blank())
    
    my.p[[2]]<-ggplot(gene.info, aes(Cell, detectedGenesPerCell)) +
      geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
      geom_jitter(aes(Cell, detectedGenesPerCell),
                  position=position_jitter(width=0.3,height=0),
                  alpha=0.1,col="orange",
                  size=2)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                axis.text.x=element_blank())
    
    my.p[[3]]<-ggplot(mito.info, aes(Cell, percent_mito)) +
      geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
      geom_jitter(aes(Cell, percent_mito),
                  position=position_jitter(width=0.3,height=0),
                  alpha=0.05,col="purple",
                  size=2)
    
    my.p[[4]]<-ggplot(ERCC.info, aes(Cell, percent_ERCC)) +
      geom_boxplot(outlier.size=0,outlier.shape = NA,col="black") +
      geom_jitter(aes(Cell, percent_ERCC),
                  position=position_jitter(width=0.3,height=0),
                  alpha=0.05,col="purple",
                  size=2)
    
    grid.arrange(grobs = my.p, ncol = 2, as.table = FALSE)
    
  }else{
    stop("please select the type of plots!")
  }
}

#' Plot UMI count vs detected genes
#' @param  object The SingCellaR object.
#' @export
#'

plot_UMIs_vs_Detected_genes<- function(object){

	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	cells.anno<-get_cells_annotation(object)
	#dev.off()
	plot(cells.anno$UMI_count,cells.anno$detectedGenesPerCell,
			pch=19,
			col=densCols(cells.anno$UMI_count,cells.anno$detectedGenesPerCell,colramp=colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[-(1:4)])),
			cex=0.75,
			main='UMIs Vs number of detected genes per cell',
			xlab='UMI count per cell',
			ylab='Number of detected genes per cell')
	abline(v=mean(cells.anno$UMI_count),lwd=2,lty=2,col="gray")
	abline(h=mean(cells.anno$detectedGenesPerCell),lwd=2,lty=2,col="gray")
	legend("bottomright","Mean", col = "gray", lty = 2:3, lwd = 2)
}

#' Plot UMIs of house keeping genes
#' @param  object The SingCellaR object.
#' @param  selected_top_genes The number of top genes. Default 100
#' @param  selected_top_cells The number of top cells. Default 300
#' @export
#'

plot_UMIs_of_HouseKeepingGenes<- function(object,selected_top_genes=100,selected_top_cells=300){

	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	umi.dat<-log1p(get_umi_count(object))
	##
	genes.info<-get_genes_metadata(object)
	genes.info<-genes.info[order(genes.info$num_detected_cells,decreasing=T),]
	hk.genes<-genes.info[1:selected_top_genes,]
	#hk.genes<-subset(genes.info,num_detected_cells >=(percentage_of_expressing_cells/100)*ncol(umi.dat))
	S.m<-umi.dat[rownames(umi.dat) %in% rownames(hk.genes),]
	S.m2plot<-S.m[,order(colMedians(as.matrix(S.m)),decreasing=T)]

	boxplot(as.matrix(S.m2plot[,1:selected_top_cells]),col="red",ylab="UMIs, Log",las=2,
			main=paste("Top ",nrow(hk.genes)," genes from top ",selected_top_cells," cells", sep=""),cex.axis=0.8,cex=0.4,xaxt='n')
}

#' Plot highly variable genes
#' @param  object The SingCellaR object.
#' @param  quantile_genes_expr_for_fitting The minimum quantile of gene expression used for model fitting. Default 0.2
#' @param  quantile_genes_cv2_for_fitting  The maximum quantile of coefficient of variation used for model fitting. Default 0.8
#' @export


plot_variable_genes<- function(object,quantile_genes_expr_for_fitting=0.20,quantile_genes_cv2_for_fitting=0.80){

	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	cells.info<-get_cells_annotation(object)
	cells.used<-subset(cells.info,IsPassed==TRUE)
	###
	normalized.umi<-get_normalized_umi(object)
	normalized.umi.used<-normalized.umi[,as.character(cells.used$Cell)]
	rm(normalized.umi)
	###
	genes_info<-get_genes_metadata(object)
	selected.genes<-subset(genes_info,IsExpress==TRUE)
	selected.genes$gene<-rownames(selected.genes)
	###
	normalized.umi.used<-normalized.umi.used[rownames(normalized.umi.used) %in% rownames(selected.genes),]
	###
	means <- Matrix::rowMeans(normalized.umi.used)
	
	block.gene.size=2000
	step.size = block.gene.size
	ends <- c(seq(from=block.gene.size,to=nrow(normalized.umi.used),by=step.size),as.integer(nrow(normalized.umi.used)))
	starts <- seq(from=1,to=nrow(normalized.umi.used),by=step.size)
	window.genes<-data.frame(start=starts,end=ends)
	####
	print("Calculate row variance..")
	pb <- txtProgressBar(max = nrow(window.genes), style = 3)
	vars<-c()
	for(i in 1:nrow(window.genes)){
		Sys.sleep(0.5);
		x.start<-window.genes$start[i]
		x.end <-window.genes$end[i]
		vars.w <- matrixStats::rowVars(as.matrix(normalized.umi.used[x.start:x.end,]))
		vars<-append(vars,vars.w)  
		setTxtProgressBar(pb, pb$getVal()+1)
	}
	
	#vars  <-matrixStats::rowVars(as.matrix(normalized.umi.used))
	#vars <- apply(as.matrix(normalized.umi.used),1,var)
	cv2 <- vars/means^2
	#################
	na.index<-which(is.na(cv2)==T)
	if(length(na.index) > 0){
		means<-means[-c(na.index)]
		cv2<-cv2[-c(na.index)]
	}
	minMeanForFit <- unname(quantile( means,quantile_genes_expr_for_fitting))
	maxVarForFit <- unname(quantile(cv2,quantile_genes_cv2_for_fitting))
	###############################################
	genesForFit<-which(means >=minMeanForFit & cv2 <=maxVarForFit)
	##################
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[genesForFit] ),cv2[genesForFit] )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"])
	####
	n.m<-data.frame(gene=names(means),means=means,cv2=cv2)
	var.f<-subset(genes_info,IsVarGenes==TRUE)
	my.var.genes<-as.character(rownames(var.f))
	var.m<-n.m[rownames(n.m) %in% my.var.genes,]

	n.genes<-length(my.var.genes)
	##############################################
	xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
	vfit <- a1/xg + a0
	plot(log(n.m$means),log(n.m$cv2),xlab="Mean expression, Log",ylab="Coefficient of variation, Log",pch=19,col="gray")
	points(log(var.m$mean),log(var.m$cv),pch=19,col="cyan",cex=0.8)
	lines( log(xg), log(vfit), col="red", lwd=3,lty=2 )
	legend("topright",paste("Variable genes : ",n.genes,sep=""), col = "cyan",pch=19)
}

#' Plot PCA's elbow plot
#' @param  object The SingCellaR object.
#' @export

plot_PCA_Elbowplot<-function(object){

	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	res_pca<-get_pca.result(object)
	#compute variance
	pr_var <- res_pca$sdev^2
	#proportion of variance explained
	prop_varex <- pr_var/sum(pr_var)
	plot(prop_varex, xlab = "Principal component",ylab = "Proportion of variance explained",type = "b")
	
	#prop_varex <- res_pca$sdev^2/sum(res_pca$sdev^2)
	#p <- qplot(1:length(prop_varex), prop_varex, alpha = I(0.5)) + theme(legend.position = "top",
	#				legend.key.height = grid::unit(0.35, "in")) +
	#		xlab("components") + ylab("Proportion of Variance Explained")
	#p
}

#' Plot TSNE with QC information
#' @param  object The SingCellaR object.
#' @param  point.size The point size
#' @export

plot_tsne_label_by_qc<-function(object,point.size=0.2){

	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	res.tsne<-get_tsne.result(object)
	p.x <- list()
	p.x[[1]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=UMI_count)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	p.x[[2]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=detectedGenesPerCell)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	p.x[[3]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=percent_mito)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	p.x[[4]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=sampleID)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	do.call(grid.arrange,p.x)
}

#' Plot TSNE with sampleID information
#' @param  object The SingCellaR object.
#' @param  point.size The point size. Default 0.2
#' @export
#' 
plot_tsne_label_by_sampleID<-function(object,point.size=0.2){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	res.tsne<-get_tsne.result(object)
	p.x <- list()
	p.x[[1]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=sampleID)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	do.call(grid.arrange,p.x)
}

#' Plot TSNE and show only selected sampleID
#' @param  object The SingCellaR object.
#' @param  selected.sampleID The vector of selected sampleIDs.
#' @param  title The title of the plot.
#' @param  point.size The point size. Default 0.5
#' @param  point.front.color The point color of selected sampleIDs. Default red
#' @param  point.base.color The point color of unselected samplesIDs.
#' @param  alpha.base The alpha parameter. Default 0.2
#' @param  title.color The title color. Default red
#' @param  title.font.size The font size of the title. Default 16
#' @export

plot_tsne_label_by_selected_sampleID<-function(object,selected.sampleID=c(),title="",point.size=0.5,
                                               point.front.color="red",point.base.color="gray",alpha.base=0.2,
                                               title.color="red",title.font.size=16){

  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(selected.sampleID)==0){
    stop("Please input selected sample IDs!")
  }
  #################################
  cell.info<-get_cells_annotation(object)
  cell.info<-subset(cell.info,IsPassed==T)
  rownames(cell.info)<-cell.info$Cell
  
  res.tsne<-get_tsne.result(object)
  #################################
  res.tsne$color<-point.base.color
  sample.index<-which(res.tsne$sampleID %in% selected.sampleID)
  res.tsne$color[sample.index]<-point.front.color
  
  res.tsne$alpha<-alpha.base
  sample.index<-which(res.tsne$sampleID %in% selected.sampleID)
  res.tsne$alpha[sample.index]<-1
  
  ##################################
  p.x <- ggplot(res.tsne,aes(TSNE1,TSNE2)) + ggtitle(title)+
    geom_point(color=res.tsne$color,size=point.size,alpha=res.tsne$alpha) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  p.x <-p.x+theme(plot.title = element_text(color=title.color, size=title.font.size, face="bold.italic",hjust = 0.5))
  p.x
}

#' Plot TSNE and display only the feature of interest
#' @param  object The SingCellaR object.
#' @param  feature The feature name. The features can be found in the cell metadata.
#' @param  point.size The point size. Default 1
#' @export

plot_tsne_label_by_a_feature_of_interest<-function(object,feature="",point.size=1){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(feature==""){
    stop("Please input a feature to plot!")
  }
  if(feature %in% colnames(object@meta.data)==T){
    ###
    res.tsne<-get_tsne.result(object)
    p.x <- list()
    p.x[[1]] <- qplot(TSNE1,TSNE2, data=res.tsne)+geom_point(aes_string(colour = feature),size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
    do.call(grid.arrange,p.x)
  }else{
    error<-paste("Couldnot find :",feature," in the column names of cell meta data.",sep="")
    stop(error)
  }
}

#' Plot UMAP with QC information
#' @param  object The SingCellaR object.
#' @param  point.size The point size. Default 0.2
#' @export
#' 
plot_umap_label_by_qc<-function(object,point.size=0.2){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	res.umap<-get_umap.result(object)
	p.x <- list()
	p.x[[1]] <- qplot(UMAP1,UMAP2, data=res.umap,colour=UMI_count)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	p.x[[2]] <- qplot(UMAP1,UMAP2, data=res.umap,colour=detectedGenesPerCell)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	p.x[[3]] <- qplot(UMAP1,UMAP2, data=res.umap,colour=percent_mito)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	p.x[[4]] <- qplot(UMAP1,UMAP2, data=res.umap,colour=sampleID)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	do.call(grid.arrange,p.x)
}

#' Plot UMAP with sampleID
#' @param  object The SingCellaR object.
#' @param  point.size The point size. Default 0.2
#' @export

plot_umap_label_by_sampleID<-function(object,point.size=0.2){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	res.umap<-get_umap.result(object)
	p.x <- list()
	p.x[[1]] <- qplot(UMAP1,UMAP2, data=res.umap,colour=sampleID)+geom_point(size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	do.call(grid.arrange,p.x)
}

#' Plot UMAP and display only selected sampleIDs
#' @param  object The SingCellaR object.
#' @param  selected.sampleID The vector of selected sampleIDs. 
#' @param  title The title of the plot.
#' @param  point.size The point size. Default 0.5
#' @param  point.front.color The point color of selected sampleIDs. Default red
#' @param  point.base.color The point color of unselected samplesIDs.
#' @param  alpha.base The alpha parameter. Default 0.2
#' @param  title.color The title color. Default red
#' @param  title.font.size The font size of the title. Default 16
#' @export

plot_umap_label_by_selected_sampleID<-function(object,selected.sampleID=c(),title="",point.size=0.5,point.front.color="red",
                                               point.base.color="gray",alpha.base=0.2,title.color="red",title.font.size=16){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(selected.sampleID)==0){
    stop("Please input selected sample IDs!")
  }
  #################################
  cell.info<-get_cells_annotation(object)
  cell.info<-subset(cell.info,IsPassed==T)
  rownames(cell.info)<-cell.info$Cell
  
  res.umap<-get_umap.result(object)
  #################################
  res.umap$color<-point.base.color
  sample.index<-which(res.umap$sampleID %in% selected.sampleID)
  res.umap$color[sample.index]<-point.front.color
  
  res.umap$alpha<-alpha.base
  sample.index<-which(res.umap$sampleID %in% selected.sampleID)
  res.umap$alpha[sample.index]<-1
  
  ##################################
  p.x <- ggplot(res.umap,aes(UMAP1,UMAP2)) + ggtitle(title)+
    geom_point(color=res.umap$color,size=point.size,alpha=res.umap$alpha) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  p.x <-p.x+theme(plot.title = element_text(color=title.color, size=title.font.size, face="bold.italic",hjust = 0.5))
  p.x
}

#' Plot UMAP and display only the feature of interest
#' @param  object The SingCellaR object.
#' @param  feature The feature name. The features can be found in the cell metadata.
#' @param  point.size The point size. Default 1
#' @export

plot_umap_label_by_a_feature_of_interest<-function(object,feature="",point.size=1){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(feature==""){
    stop("Please input a feature to plot!")
  }
  if(feature %in% colnames(object@meta.data)==T){
    ###
    res.umap<-get_umap.result(object)
    p.x <- list()
    p.x[[1]] <- qplot(UMAP1,UMAP2, data=res.umap)+geom_point(aes_string(colour = feature),size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
    do.call(grid.arrange,p.x)
  }else{
    error<-paste("Couldnot find :",feature," in the column names of cell meta data.",sep="")
    stop(error)
  }
  
}

#' Plot TSNE with gene expression
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of gene names
#' @param  point.size The point size. Default 0.3
#' @param  point.color1 The first color of gene expression gradient
#' @param  point.color2 The second color of gene expression gradient.
#' @export

plot_tsne_label_by_genes<-function(object,gene_list=c(),point.size=0.3,point.color1="blue",point.color2="red"){

	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	#####
	res.tsne<-get_tsne.result(object)
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(res.tsne$Cell)]
	#####
	my.dat<-res.tsne[,(1:3)]
	for(i in 1:length(gene_list)){
	  genes.x<-gene_list[i]
	  g.ind<-which(rownames(umi.used)==genes.x)
	  my.sub.umi<-as.data.frame(log1p(umi.used[g.ind,]))
	  colnames(my.sub.umi)<-genes.x
	  #################################
	  my.dat<-cbind(my.dat,my.sub.umi)
	}
	X.colnames<-colnames(my.dat)
	X.colnames<-gsub("-","_",X.colnames)
	X.colnames<-gsub("\\+","_",X.colnames)
	colnames(my.dat)<-X.colnames

	my.fg.names<-colnames(my.dat)
	my.fg.names<-my.fg.names[4:length(my.fg.names)]

	p.x <- list()
	for(j in 1:length(my.fg.names)){
		genes.x<-my.fg.names[j]
		genes.x<-gsub("-","_",genes.x)
		genes.x<-gsub("\\+","_",genes.x)
		p.x[[j]] <- ggplot(my.dat,aes(TSNE1,TSNE2)) + geom_point(aes_string(colour = genes.x),size=point.size) + scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
	}
	do.call(grid.arrange,p.x)
}

#' Plot UMAP with gene expression
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of gene names
#' @param  point.size The point size. Default 0.3
#' @param  point.color1 The first color of gene expression gradient. Default blue
#' @param  point.color2 The second color of gene expression gradient. Default red
#' @export

plot_umap_label_by_genes<-function(object,gene_list=c(),point.size=0.3,point.color1="blue",point.color2="red"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the scRNAseq object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  #####
  res.umap<-get_umap.result(object)
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(res.umap$Cell)]
  #####
  my.dat<-res.umap[,c("Cell","UMAP1","UMAP2")]
  for(i in 1:length(gene_list)){
    genes.x<-gene_list[i]
    g.ind<-which(rownames(umi.used)==genes.x)
    my.sub.umi<-as.data.frame(log1p(umi.used[g.ind,]))
    colnames(my.sub.umi)<-genes.x
    #################################
    my.dat<-cbind(my.dat,my.sub.umi)
  }
  X.colnames<-colnames(my.dat)
  X.colnames<-gsub("-","_",X.colnames)
  X.colnames<-gsub("\\+","_",X.colnames)
  colnames(my.dat)<-X.colnames
  
  my.fg.names<-colnames(my.dat)
  my.fg.names<-my.fg.names[4:length(my.fg.names)]
  
  p.x <- list()
  for(j in 1:length(my.fg.names)){
    genes.x<-my.fg.names[j]
    genes.x<-gsub("-","_",genes.x)
    genes.x<-gsub("\\+","_",genes.x)
    p.x[[j]] <- ggplot(my.dat,aes(UMAP1,UMAP2)) + geom_point(aes_string(colour = genes.x),size=point.size) + scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  }
  do.call(grid.arrange,p.x)
}

#' Plot diffusion map with gene expression
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of gene names
#' @param  point.size The point size. Default 0.3
#' @param  point.color1 The first color of gene expression gradient. Default blue
#' @param  point.color2 The second color of gene expression gradient. Default red
#' @param  x_eigenVal The diffusion map dimension on x axis. Default DC1
#' @param  y_eigenVal The diffusion map dimension on y axis. Default DC2
#' @export

plot_diffusionmap_label_by_genes<-function(object,gene_list=c(),point.size=0.3,
                                           point.color1="blue",point.color2="red",
                                           x_eigenVal = "DC1",y_eigenVal = "DC2"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  #####
  res.dfm<-get_dfm.result(object)
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(res.dfm$Cell)]
  #####
  my.dat<-res.dfm[,c("Cell",x_eigenVal,y_eigenVal)]
  for(i in 1:length(gene_list)){
    genes.x<-gene_list[i]
    g.ind<-which(rownames(umi.used)==genes.x)
    my.sub.umi<-as.data.frame(log1p(umi.used[g.ind,]))
    colnames(my.sub.umi)<-genes.x
    #################################
    my.dat<-cbind(my.dat,my.sub.umi)
  }
  X.colnames<-colnames(my.dat)
  X.colnames<-gsub("-","_",X.colnames)
  X.colnames<-gsub("\\+","_",X.colnames)
  colnames(my.dat)<-X.colnames
  
  my.fg.names<-colnames(my.dat)
  my.fg.names<-my.fg.names[4:length(my.fg.names)]
  
  p.x <- list()
  for(j in 1:length(my.fg.names)){
    genes.x<-my.fg.names[j]
    genes.x<-gsub("-","_",genes.x)
    genes.x<-gsub("\\+","_",genes.x)
    p.x[[j]] <- ggplot(my.dat,aes_string(x_eigenVal,y_eigenVal)) + 
      geom_point(aes_string(colour = genes.x),size=point.size) + 
      scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
      theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  }
  do.call(grid.arrange,p.x)
}

#' Plot TSNE with the expression of selected gene set
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of gene names.
#' @param  gene_set_name The name of gene set.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  point.size The point size. Default 0.3
#' @param  point.color1 The first color of gene expression gradient. Default blue
#' @param  point.color2 The second color of gene expression gradient. Default red
#' @export

plot_tsne_label_by_a_signature_gene_set<-function(object,gene_list=c(),gene_set_name=c(),
                                                  isNormalizedByHouseKeeping=F,point.size=1,
                                                  point.color1="black",point.color2="red"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	if(length(gene_list) < 2){
		stop("Required more than 2 genes!")
	}
	if(length(gene_set_name)==0){
		stop("Required a name of gene set!")
	}
	#####
	res.tsne<-get_tsne.result(object)
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(res.tsne$Cell)]
	#####
	my.dat<-res.tsne[,(1:3)]
	#####
	exprs<-as.matrix(umi.used[rownames(umi.used) %in% gene_list,])
	########get house keeping genes############
	if(isNormalizedByHouseKeeping==T){
		genes.exp.sum<-Matrix::rowSums(umi.used)
		genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
		hk.genes<-names(genes.exp.sum[1:100])
		hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
		#####
		MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
		MM_score<-MM
	}else{
		MM <- Matrix::colSums(exprs)/(res.tsne$UMI_count)
		MM_score<-MM
	}
	#####
	Genes.score<-data.frame(Genes_score=MM_score)
	updated.tsne<-cbind(res.tsne,Genes.score)
	####################################
	ggplot(updated.tsne,aes(TSNE1,TSNE2, color=Genes_score)) + geom_point(size=point.size)+ggtitle(gene_set_name)+theme(plot.title = element_text(hjust = 0.5))+ scale_colour_gradientn(colours = c("gray85",point.color1,point.color2),values=c(0,0.1,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
}

#' Plot UMAP with the expression of selected gene set
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of gene names.
#' @param  gene_set_name The name of gene set.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  point.size The point size. Default 0.3
#' @param  point.color1 The first color of gene expression gradient. Default blue
#' @param  point.color2 The second color of gene expression gradient. Default red
#' @export
#' 
plot_umap_label_by_a_signature_gene_set<-function(object,gene_list=c(),gene_set_name=c(),
                                                  isNormalizedByHouseKeeping=F,point.size=1,
                                                  point.color1="black",point.color2="red"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	if(length(gene_list) < 2){
		stop("Required more than 2 genes!")
	}
	if(length(gene_set_name)==0){
		stop("Required a name of gene set!")
	}
	#####
	res.umap<-get_umap.result(object)
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(res.umap$Cell)]
	#####
	my.dat<-res.umap[,(1:3)]
	#####
	exprs<-as.matrix(umi.used[rownames(umi.used) %in% gene_list,])
	########get house keeping genes############
	if(isNormalizedByHouseKeeping==T){
		genes.exp.sum<-Matrix::rowSums(umi.used)
		genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
		hk.genes<-names(genes.exp.sum[1:100])
		hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
		#####
		MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
		MM_score<-MM
	}else{
		MM <- Matrix::colSums(exprs)/(res.umap$UMI_count)
		MM_score<-MM
	}
	#####
	Genes.score<-data.frame(Genes_score=MM_score)
	updated.umap<-cbind(res.umap,Genes.score)
	####################################
	ggplot(updated.umap,aes(UMAP1,UMAP2, color=Genes_score)) + geom_point(size=point.size)+ggtitle(gene_set_name)+theme(plot.title = element_text(hjust = 0.5))+ scale_colour_gradientn(colours = c("gray85",point.color1,point.color2),values=c(0,0.1,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
}

#' Plot TSNE with the expression of multiple gene sets
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file for gene sets.
#' @param  show_gene_sets The names of selected gene sets.
#' @param  custom_color The custom color names.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  point.size The point size. Default 2
#' @param  showLegend is logical. If TRUE, the figure legend will be shown.
#' @param  background.color The background color of the plot. Default white
#' @export

plot_tsne_label_by_multiple_gene_sets<-function(object,gmt.file=c(),show_gene_sets=c(),custom_color=c(),
													   isNormalizedByHouseKeeping=T,point.size=2,showLegend=T,background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gmt.file)==0){
		stop("Required the path to GMT file!")
	}
	if(length(show_gene_sets)==0){
		stop("Required gene set names!")
	}
	#if(length(show_gene_sets)==1){
	#	stop("Required at least two gene sets")
	#}
	if(length(show_gene_sets) > 10){
		stop("This is over limitation at 10 gene sets!")
	}
	if(length(custom_color) > 0){
		if(length(custom_color) !=length(show_gene_sets)){
			stop("The number of colors is not equal to the the number of input gene sets.")
		}
	}
	#####
	signature.sets<-get_gmtGeneSets(gmt.file)
	#####
	check.1<-match(show_gene_sets,names(signature.sets))
	check.2<-sum(length(which(is.na(check.1))))
	if(check.2 > 0){
		stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
	}
	res.tsne<-get_tsne.result(object)
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(res.tsne$Cell)]
	#####
	my.dat<-res.tsne[,c("Cell","TSNE1","TSNE2")]
	#####
	genes.exp.sum<-Matrix::rowSums(umi.used)
	genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
	hk.genes<-names(genes.exp.sum[1:100])
	hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
	#####
	Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
	for(i in 1:length(show_gene_sets)){
		set.name<-show_gene_sets[i]
		genes.x<-signature.sets[[set.name]]
		exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
		
		if(isNormalizedByHouseKeeping==T){
			MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
		}else{
			MM<-Matrix::colSums(exprs)/res.tsne$UMI_count
		}
		MM.f<-data.frame(score=MM)
		colnames(MM.f)<-set.name
		Genes.score<-cbind(Genes.score,MM.f)
	}
	Genes.score<-Genes.score[-c(1)]
	####Assign colors#####################
	if(length(custom_color)==0){
		multi.colors<-rainbow(ncol(Genes.score))
	}else{
		multi.colors<-custom_color
	}
	######################################
	my.colors.list<-list()
	for(j in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[j])
		x.color<-multi.colors[j]
		p1<-ggplot(my.dat,aes(TSNE1,TSNE2)) + geom_point(aes(colour = Genes.score[,j])) + scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
		My_color<-ggplot_build(p1)$data[[1]]$colour
		my.colors.list[[set.name]]<-My_color
	}
	
	my.alpha.list<-list()
	for(j in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[j])
		my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
		my.alpha.list[[set.name]]<-my.alpha
	}
	###################################
	p<-ggplot(my.dat,aes(TSNE1,TSNE2))
	for(k in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[k])
		p<-p+geom_point(color=my.colors.list[[set.name]], size = point.size,alpha=my.alpha.list[[set.name]],stroke = 0)
	}
	
	p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
			theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
			theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
			theme(plot.background = element_rect(fill = background.color))
	
	if(showLegend==T){
		df<-data.frame(Name=show_gene_sets)
		df<-as.matrix(df)
	
			tp<-ggtexttable(df,rows = NULL,
			theme = ttheme(
					colnames.style = colnames_style(fill = "white"),
					tbody.style = tbody_style(fill = multi.colors)
			)
		)
		return(ggarrange(p,tp,
					ncol = 2, nrow = 1,
					widths = c(1,0.2)))
	}else{
		return(p)
	}
}

#' Plot UMAP with the expression of multiple gene sets
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file for gene sets.
#' @param  show_gene_sets The names of selected gene sets.
#' @param  custom_color The custom color names.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  point.size The point size. Default 2
#' @param  showLegend is logical. If TRUE, the figure legend will be shown.
#' @param  background.color The background color of the plot. Default white
#' @export

plot_umap_label_by_multiple_gene_sets<-function(object,gmt.file=c(),show_gene_sets=c(),custom_color=c(),
														isNormalizedByHouseKeeping=T,point.size=2,showLegend=T,background.color="white"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gmt.file)==0){
    stop("Required the path to GMT file!")
  }
  if(length(show_gene_sets)==0){
    stop("Required gene set names!")
  }
  #if(length(show_gene_sets)==1){
  #	stop("Required at least two gene sets")
  #}
  if(length(show_gene_sets) > 10){
    stop("This is over limitation at 10 gene sets!")
  }
  if(length(custom_color) > 0){
    if(length(custom_color) !=length(show_gene_sets)){
      stop("The number of colors is not equal to the the number of input gene sets.")
    }
  }
  #####
  signature.sets<-get_gmtGeneSets(gmt.file)
  #####
  check.1<-match(show_gene_sets,names(signature.sets))
  check.2<-sum(length(which(is.na(check.1))))
  if(check.2 > 0){
    stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
  }
  res.umap<-get_umap.result(object)
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(res.umap$Cell)]
  #####
  my.dat<-res.umap[,c("Cell","UMAP1","UMAP2")]
  #####
  genes.exp.sum<-Matrix::rowSums(umi.used)
  genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
  hk.genes<-names(genes.exp.sum[1:100])
  hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
  #####
  Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
  for(i in 1:length(show_gene_sets)){
    set.name<-show_gene_sets[i]
    genes.x<-signature.sets[[set.name]]
    exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
    
    if(isNormalizedByHouseKeeping==T){
      MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
    }else{
      MM<-Matrix::colSums(exprs)/res.umap$UMI_count
    }
    MM.f<-data.frame(score=MM)
    colnames(MM.f)<-set.name
    Genes.score<-cbind(Genes.score,MM.f)
  }
  Genes.score<-Genes.score[-c(1)]
  ####Assign colors#####################
  if(length(custom_color)==0){
    multi.colors<-rainbow(ncol(Genes.score))
  }else{
    multi.colors<-custom_color
  }
  ######################################
  my.colors.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    x.color<-multi.colors[j]
    p1<-ggplot(my.dat,aes(UMAP1,UMAP2)) + geom_point(aes(colour = Genes.score[,j])) + scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
    My_color<-ggplot_build(p1)$data[[1]]$colour
    my.colors.list[[set.name]]<-My_color
  }
  
  my.alpha.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
    my.alpha.list[[set.name]]<-my.alpha
  }
  ###################################
  p<-ggplot(my.dat,aes(UMAP1,UMAP2))
  for(k in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[k])
    p<-p+geom_point(color=my.colors.list[[set.name]], size = point.size,alpha=my.alpha.list[[set.name]],stroke = 0)
  }
  
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(plot.background = element_rect(fill = background.color))
  
  if(showLegend==T){
  		df<-data.frame(Name=show_gene_sets)
  		df<-as.matrix(df)
  
  		tp<-ggtexttable(df,rows = NULL,
                  theme = ttheme(
                    colnames.style = colnames_style(fill = "white"),
                    tbody.style = tbody_style(fill = multi.colors)
                  )
  		)
  		return(ggarrange(p,tp,
                   ncol = 2, nrow = 1,
                   widths = c(1,0.2)))
   }else{
   		return(p)
   }            
}

#' Plot diffusion map with the expression of multiple gene sets
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file for gene sets.
#' @param  show_gene_sets The names of selected gene sets.
#' @param  custom_color The custom color names.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  point.size The point size. Default 1
#' @param  x_eigenVal The diffusion map dimension on x axis. Default DC1
#' @param  y_eigenVal The diffusion map dimension on y axis. Default DC2
#' @param  showLegend is logical. If TRUE, the figure legend will be shown.
#' @param  background.color The background color of the plot. Default white
#' @export

plot_diffusionmap_label_by_multiple_gene_sets<-function(object,gmt.file=c(),show_gene_sets=c(),
                                                custom_color=c(),isNormalizedByHouseKeeping=T,
                                                x_eigenVal = "DC1",y_eigenVal = "DC2",point.size=1,
                                                showLegend=T,background.color="white"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gmt.file)==0){
    stop("Required the path to GMT file!")
  }
  if(length(show_gene_sets)==0){
    stop("Required gene set names!")
  }
  #if(length(show_gene_sets)==1){
  #	stop("Required at least two gene sets")
  #}
  if(length(show_gene_sets) > 10){
    stop("This is over limitation at 10 gene sets!")
  }
  if(length(custom_color) > 0){
    if(length(custom_color) !=length(show_gene_sets)){
      stop("The number of colors is not equal to the the number of input gene sets.")
    }
  }
  #####
  signature.sets<-get_gmtGeneSets(gmt.file)
  #####
  check.1<-match(show_gene_sets,names(signature.sets))
  check.2<-sum(length(which(is.na(check.1))))
  if(check.2 > 0){
    stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
  }
  res.dfm<-get_dfm.result(object)
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(res.dfm$Cell)]
  #####
  my.dat<-res.dfm[,c("Cell",x_eigenVal,y_eigenVal)]
  #####
  genes.exp.sum<-Matrix::rowSums(umi.used)
  genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
  hk.genes<-names(genes.exp.sum[1:100])
  hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
  #####
  Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
  for(i in 1:length(show_gene_sets)){
    set.name<-show_gene_sets[i]
    genes.x<-signature.sets[[set.name]]
    exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
    
    if(isNormalizedByHouseKeeping==T){
      MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
    }else{
      MM<-Matrix::colSums(exprs)/res.dfm$UMI_count
    }
    MM.f<-data.frame(score=MM)
    colnames(MM.f)<-set.name
    Genes.score<-cbind(Genes.score,MM.f)
  }
  Genes.score<-Genes.score[-c(1)]
  ####Assign colors#####################
  if(length(custom_color)==0){
    multi.colors<-rainbow(ncol(Genes.score))
  }else{
    multi.colors<-custom_color
  }
  ######################################
  my.colors.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    x.color<-multi.colors[j]
    p1<-ggplot(my.dat,aes_string(x_eigenVal,y_eigenVal)) + geom_point(aes(colour = Genes.score[,j])) + 
      scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
    My_color<-ggplot_build(p1)$data[[1]]$colour
    my.colors.list[[set.name]]<-My_color
  }
  
  my.alpha.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
    my.alpha.list[[set.name]]<-my.alpha
  }
  ###################################
  p<-ggplot(my.dat,aes_string(x_eigenVal,y_eigenVal))
  for(k in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[k])
    p<-p+geom_point(color=my.colors.list[[set.name]], size = point.size,alpha=my.alpha.list[[set.name]],stroke = 0)
  }
  
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
    theme(plot.background = element_rect(fill = background.color))
  
  if(showLegend==T){
    df<-data.frame(Name=show_gene_sets)
    df<-as.matrix(df)
    
    tp<-ggtexttable(df,rows = NULL,
                    theme = ttheme(
                      colnames.style = colnames_style(fill = "white"),
                      tbody.style = tbody_style(fill = multi.colors)
                    )
    )
    return(ggarrange(p,tp,
                     ncol = 2, nrow = 1,
                     widths = c(1,0.2)))
  }else{
    return(p)
  }
}

#plot_tsne_label_by_density_cluster<-function(object,point.size=2,cluster.font.size=2,add.cluster.id=FALSE,cluster.font.color="yellow"){
#	
#	if(!is(object,"SingCellaR")){
#		stop("Need to initialize the SingCellaR object")
#	}
#	###
#	res.tsne<-get_tsne.result(object)
#	cluster.id<-gsub("cl_","",res.tsne$density_cluster)
#	p.x <- list()
#	if(add.cluster.id==TRUE){
#		p.x[[1]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=density_cluster)+geom_point(aes(fill=density_cluster),colour="black",pch=21,size=point.size)+ geom_text(aes(label = cluster.id),
#				size = cluster.font.size,colour=cluster.font.color)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
#	}else{
#		p.x[[1]] <- qplot(TSNE1,TSNE2, data=res.tsne,colour=density_cluster)+geom_point(aes(fill=density_cluster),colour="black",pch=21,size=point.size)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
#		
#	}
#	do.call(grid.arrange,p.x)
#}

#' Plot TSNE with identified clusters
#' @param  object The SingCellaR object.
#' @param  show_method The clustering method names used to generate clusters.
#' @param  point.size The point size. Default 2
#' @param  mark.clusters is logical. If TRUE, the number of identified clusters will be displayed.
#' @param  mark.font.size The font size of the cluster name.
#' @param  mark.font.color The font color of the cluster name.
#' @param  add.cluster.in.cells is logical. If TRUE, the cluster name will be displayed in each data point.
#' @param  cluster.font.size The font size of cluster name in each data point.
#' @param  cluster.font.color The font color of cluster name in each data point.
#' @export


plot_tsne_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                      point.size=2,mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",add.cluster.in.cells=FALSE,
                                      cluster.font.size=1,cluster.font.color="yellow"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	show_method <- match.arg(show_method)
	###
	res.tsne<-get_tsne.result(object)
	clusters.info<-get_clusters(object)
	
	if(show_method=="walktrap"){
		clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="louvain"){
		clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)
		colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_walktrap"){
	  clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_louvain"){
	  clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="merged_kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else{
		stop("Need community detection or clustering-method names!")
	}
	res.tsne<-merge(res.tsne,clusters.info)
	cluster.id<-gsub("cl","",res.tsne$cluster)
	
	if(add.cluster.in.cells==TRUE){
		p<- qplot(TSNE1,TSNE2, data=res.tsne,colour=cluster)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+ 
				geom_text(aes(label = add.cluster.in.cells),size = cluster.font.size,colour=cluster.font.color)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
		
		if(mark.clusters==TRUE){
			edata <- res.tsne[,c("TSNE1","TSNE2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
		
	}else{
		p<- qplot(TSNE1,TSNE2, data=res.tsne,colour=cluster)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
		
		if(mark.clusters==TRUE){
			edata <- res.tsne[,c("TSNE1","TSNE2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
	}
	p<-p+ggtitle(show_method)
	return(p)
}

#' Plot UMAP with identified clusters
#' @param  object The SingCellaR object.
#' @param  show_method The clustering method names used to generate clusters.
#' @param  point.size The point size. Default 2
#' @param  mark.clusters is logical. If TRUE, the number of identified clusters will be displayed.
#' @param  mark.font.size The font size of the cluster name.
#' @param  mark.font.color The font color of the cluster name.
#' @param  add.cluster.in.cells is logical. If TRUE, the cluster name will be displayed in each data point.
#' @param  cluster.font.size The font size of cluster name in each data point.
#' @param  cluster.font.color The font color of cluster name in each data point.
#' @export

plot_umap_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                      point.size=2,mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",
                                      add.cluster.in.cells=FALSE,cluster.font.size=1,cluster.font.color="yellow"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	show_method <- match.arg(show_method)
	###
	res.umap<-get_umap.result(object)
	clusters.info<-get_clusters(object)
	
	if(show_method=="walktrap"){
	  clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="louvain"){
	  clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_walktrap"){
	  clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_louvain"){
	  clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else{
	  stop("Need community detection or clustering-method names!")
	}
	res.umap<-merge(res.umap,clusters.info)
	cluster.id<-gsub("cl","",res.umap$cluster)
	
	if(add.cluster.in.cells==TRUE){
		p<- qplot(UMAP1,UMAP2, data=res.umap,colour=cluster)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+ 
				geom_text(aes(label = add.cluster.in.cells),size = cluster.font.size,colour=cluster.font.color)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
		
		if(mark.clusters==TRUE){
			edata <- res.umap[,c("UMAP1","UMAP2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
		
	}else{
		p<- qplot(UMAP1,UMAP2, data=res.umap,colour=cluster)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
		
		if(mark.clusters==TRUE){
			edata <- res.umap[,c("UMAP1","UMAP2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
	}
	p<-p+ggtitle(show_method)
	return(p)
}

#' Plot diffusion map with identified clusters
#' @param  object The SingCellaR object.
#' @param  show_method The clustering method names used to generate clusters.
#' @param  point.size The point size. Default 2
#' @param  mark.clusters is logical. If TRUE, the number of identified clusters will be displayed.
#' @param  mark.font.size The font size of the cluster name.
#' @param  mark.font.color The font color of the cluster name.
#' @param  add.cluster.in.cells is logical. If TRUE, the cluster name will be displayed in each data point.
#' @param  cluster.font.size The font size of cluster name in each data point.
#' @param  cluster.font.color The font color of cluster name in each data point.
#' @export

plot_diffusionmap_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                              x_eigenVal="DC1",y_eigenVal="DC2",point.size=2,
                                              mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",
                                              add.cluster.in.cells=FALSE,cluster.font.size=1,cluster.font.color="yellow"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  show_method <- match.arg(show_method)
  ###
  res.dfm<-get_dfm.result(object)
  clusters.info<-get_clusters(object)
  
  if(show_method=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="merged_walktrap"){
    clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="merged_louvain"){
    clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="merged_kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else{
    stop("Need community detection or clustering-method names!")
  }
  res.dfm<-merge(res.dfm,clusters.info)
  cluster.id<-gsub("cl","",res.dfm$cluster)
  
  if(add.cluster.in.cells==TRUE){
    p<-ggplot(aes_string(x=x_eigenVal,y=y_eigenVal),data=res.dfm,colour=cluster)+
       geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+
       geom_text(aes(label = add.cluster.in.cells),size = cluster.font.size,colour=cluster.font.color)+
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
       theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
    
    if(mark.clusters==TRUE){
      edata <- res.dfm[,c(x_eigenVal,y_eigenVal)]
      edata$cluster.id<-cluster.id
      colnames(edata) <- c('x', "y", "z")
      center <- aggregate(cbind(x,y) ~ z, data = edata, median)
      p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
    }
    
  }else{
    p<- ggplot(aes_string(x=x_eigenVal,y=y_eigenVal), data=res.dfm,colour=cluster)+
      geom_point(aes(fill=cluster),colour="black",pch=21,size=point.size)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
    
    if(mark.clusters==TRUE){
      edata <- res.dfm[,c(x_eigenVal,y_eigenVal)]
      edata$cluster.id<-cluster.id
      colnames(edata) <- c('x', "y", "z")
      center <- aggregate(cbind(x,y) ~ z, data = edata, median)
      p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
    }
  }
  p<-p+ggtitle(show_method)
  return(p)
}

#' Plot violin-plot for marker genes
#' @param  object The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  gene_list The vector of gene names.
#' @param  take_log2 is logical. If TRUE, the expression shows in the log2 scale.
#' @param  xlab.text.size The font size of the label in the x-axis. Default 5
#' @param  point.size The point size. Default 0.2
#' @param  point.alpha The alpha parameter. Default 0.1
#' @export
#' 
plot_violin_for_marker_genes<-function(object,cluster.type=c("louvain","walktrap","kmeans"),gene_list=c(),
                                       take_log2=T,xlab.text.size=5,point.size=0.2,point.alpha=0.1){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  
  cluster.type <- match.arg(cluster.type)
  ####################################
  umi.dat<-get_normalized_umi(object)
  clusters.info<-get_clusters(object)
  ####################################
  if(cluster.type=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster")]
  }else if(cluster.type=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster")]
  }else if(cluster.type=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)
    clusters.info<-clusters.info[,c("Cell","kmeans_cluster")]
  }else{
    stop("Need a clustering-method name!")
  }
  
  clusters.ids<-1:length(unique(clusters.info[,2]))
  my.p<-list()
  
  for(i in 1:length(gene_list)){
    my.gene<-as.character(gene_list)[i]
    my.sub.dat<-umi.dat[rownames(umi.dat)==my.gene,]
    my.new.dat<-data.frame()
    for(k in 1:length(clusters.ids)){
      cluster.id<-paste("cl",k,sep="")
      x.cell.names<-subset(clusters.info,clusters.info[,2]==cluster.id)
      x.expr<-my.sub.dat[names(my.sub.dat) %in% as.character(x.cell.names$Cell)]
      ##################
      exp.val<-as.numeric(x.expr)
      n.exp<-length(exp.val[exp.val>0])
      n.total<-length(exp.val)
      my.label<-paste(cluster.id,"\n",n.exp,"/",n.total,sep="")
      ##################
      y.legend<-""
      ##################
      if(take_log2==T){
        my.f<-data.frame(cluster=my.label,normUMI=log2(as.numeric(x.expr)+1))
        y.legend<-"log2(normalized UMI)"
      }else{
        my.f<-data.frame(cluster=my.label,normUMI=as.numeric(x.expr))
        y.legend<-"normalized UMI"
      }
      my.new.dat<-rbind(my.new.dat,my.f)
    }
    my.p[[i]]<-ggplot(
      data = my.new.dat,aes(x=cluster, y=normUMI,fill=cluster))+labs(y=y.legend)+
      geom_violin(trim=TRUE,scale="width")+
      geom_jitter(height = 0,size=point.size,alpha = point.alpha)+
      ggtitle(my.gene)+
      guides(fill=FALSE)+
      theme(axis.text=element_text(size=xlab.text.size))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } 
  do.call(grid.arrange,my.p)
}

#' Plot violin-plot for differential genes
#' @param  objectA The SingCellaR objectA.
#' @param  objectB The SingCellaR objectB.
#' @param  cellsA The vector of cell names in group A.
#' @param  cellsB The vector of cell names in group B.
#' @param  gene_list The vector of gene names.
#' @param  take_log2 is logical. If TRUE, the expression shows in the log2 scale.
#' @param  groupA.name The group A name.
#' @param  groupB.name The group B name.
#' @param  xlab.text.size The font size of the label in the x-axis. Default 5
#' @param  point.size The point size. Default 0.2
#' @param  point.alpha The alpha parameter. Default 0.1
#' @export

plot_violin_for_differential_genes<-function(objectA=objectA,objectB=objectB,cellsA=c(),cellsB=c(),gene_list=c(),
                                             take_log2=T,groupA.name="",groupB.name="",
                                             xlab.text.size=5,point.size=0.5,point.alpha=0.1){
  
  if(!is(objectA,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(!is(objectB,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  if(length(cellsA)==0 | length(cellsB)==0){
    stop("Required input cells!")
  }
  ####################################
  umiA.dat<-get_normalized_umi(objectA)
  umiB.dat<-get_normalized_umi(objectB)
  
  cellsA.m<-umiA.dat[,colnames(umiA.dat) %in% cellsA]
  cellsB.m<-umiB.dat[,colnames(umiB.dat) %in% cellsB]
  
  my.p<-list()
  
  for(i in 1:length(gene_list)){
    
    my.gene<-as.character(gene_list)[i]
    ##################
    A.expr<-cellsA.m[rownames(cellsA.m)==my.gene,]
    A.exp.val<-as.numeric(A.expr)
    A.n<-length(A.exp.val[A.exp.val>0])
    A.total<-length(A.exp.val)
    if(groupA.name==""){
      A.label<-paste("groupA","\n",A.n,"/",A.total,sep="")
    }else{
      A.label<-paste(groupA.name,"\n",A.n,"/",A.total,sep="")
    }
    ##################
    B.expr<-cellsB.m[rownames(cellsB.m)==my.gene,]
    B.exp.val<-as.numeric(B.expr)
    B.n<-length(B.exp.val[B.exp.val>0])
    B.total<-length(B.exp.val)
    
    if(groupB.name==""){
      B.label<-paste("groupB","\n",B.n,"/",B.total,sep="")
    }else{
      B.label<-paste(groupB.name,"\n",B.n,"/",B.total,sep="")
    }
    ##################
    y.legend<-""
    ##################
    if(take_log2==T){
      A.f<-data.frame(cell=A.label,normUMI=log2(as.numeric(A.expr)+1))
      B.f<-data.frame(cell=B.label,normUMI=log2(as.numeric(B.expr)+1))
      y.legend<-"log2(normalized UMI)"
    }else{
      A.f<-data.frame(cell=A.label,normUMI=as.numeric(A.expr))
      B.f<-data.frame(cell=B.label,normUMI=as.numeric(B.expr))
      y.legend<-"normalized UMI"
    }
    
    my.new.dat<-rbind(A.f,B.f)
    
    my.p[[i]]<-ggplot(
      data = my.new.dat,aes(x=cell, y=normUMI,fill=cell))+labs(y=y.legend)+
      geom_violin(trim=TRUE,scale="width")+
      geom_jitter(height = 0,size=point.size,alpha = point.alpha)+
      ggtitle(my.gene)+
      guides(fill=FALSE)+
      theme(axis.text=element_text(size=xlab.text.size))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } 
  do.call(grid.arrange,my.p)
}

#' Plot heatmap for differential genes
#' @param  objectA The SingCellaR objectA.
#' @param  objectB The SingCellaR objectB.
#' @param  cellsA The vector of cell names in group A.
#' @param  cellsB The vector of cell names in group B.
#' @param  gene_list The vector of gene names.
#' @param  groupA.name The group A cell names.
#' @param  groupB.name The group B cell names.
#' @param  show.cell_annotation The vector of cell annotation names
#' @param  isClusterByRow is logical. If TRUE, the clustering by row will be performed.
#' @param  isClusterByCol is logical. If TRUE, the clustering by column will be performed.
#' @param  distance_row The distance method for row. Default euclidean
#' @param  distance_column The distance method for column. Default euclidean
#' @param  show_row_dend is logical. If TRUE, the row dendrogram will be shown.
#' @param  show_column_dend is logical. If TRUE, the column dendrogram will be shown.
#' @param  col.scaled.min The minimum scaled value. Default -2.5
#' @param  col.scaled.max The maximum scaled value. Default 2.5
#' @param  col.low The low color gradient. Default blue
#' @param  col.mid The mid color gradient. Default black
#' @param  col.high The high color gradient. Default red
#' @param  rowFont.size The row font size. Default 6
#' @export

plot_heatmap_for_differential_genes<-function(objectA=objectA,objectB=objectB,cellsA=c(),cellsB=c(),gene_list=c(),
                                              groupA.name="",groupB.name="",show.cell_annotation=c(),
                                              isClusterByRow=T,isClusterByCol=T,distance_row="euclidean",
                                              distance_column="euclidean",show_row_dend = T,show_column_dend = T,
                                              col.scaled.min = -2.5,col.scaled.max = 2.5,col.low="blue",col.mid = "black",
                                              col.high = "red",rowFont.size=6){
  
  if(!is(objectA,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(!is(objectB,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  if(length(cellsA)==0 | length(cellsB)==0){
    stop("Required input cells!")
  }
  ####################################
  umiA.dat<-get_normalized_umi(objectA)
  umiB.dat<-get_normalized_umi(objectB)
  
  cellsA.m<-umiA.dat[as.character(gene_list),colnames(umiA.dat) %in% cellsA]
  cellsB.m<-umiB.dat[as.character(gene_list),colnames(umiB.dat) %in% cellsB]
  
  cbind.umi<-cbind(as.matrix(cellsA.m),as.matrix(cellsB.m))
  cbind.umi.log<-log1p(cbind.umi)
  mat_scaled = t(apply(cbind.umi.log, 1, scale))
  ####################################
  cellA.anno<-get_cells_annotation(objectA)
  rownames(cellA.anno)<-cellA.anno$Cell
  cellB.anno<-get_cells_annotation(objectB)
  rownames(cellB.anno)<-cellB.anno$Cell
  
  cellA.anno<-cellA.anno[cellsA,]
  if(groupA.name==""){
    cellA.anno$cell_group<-"A"
  }else{
    cellA.anno$cell_group<-groupA.name
  }
  cellB.anno<-cellB.anno[cellsB,]
  if(groupB.name==""){
    cellB.anno$cell_group<-"B"
  }else{
    cellB.anno$cell_group<-groupB.name
  }
  
  df<-rbind(cellA.anno,cellB.anno)
  if(length(show.cell_annotation)==0){
    df<-df["cell_group"]
    col.xx <- c("blue","brown")
    names(col.xx) <- names(table(df$cell_group))
    cluster.annotation = HeatmapAnnotation(df = df, col = list(cell_group = col.xx))
  }else{
    show.cell_annotation<-c(show.cell_annotation,"cell_group")
    df<-df[show.cell_annotation]
    if(ncol(df) > 0){
      col.list<-list()
      for(i in 1:ncol(df)){
        col.name<-colnames(df)[i]
        col.xx<-colorRampPalette(brewer.pal(10,"Paired"))(length(names(table(df[i]))))
        names(col.xx)<-names(table(df[i]))
        col.list[[col.name]]<-col.xx
      }
      col.yy <- c("blue","brown")
      names(col.yy) <- names(table(df$cell_group))
      col.list[["cell_group"]]<-col.yy
      cluster.annotation = HeatmapAnnotation(df = df,col=col.list)
    }
  
  }
  ###############
  Heatmap(mat_scaled, name ="expression",col = colorRamp2(c(col.scaled.min, 0, col.scaled.max), c(col.low,col.mid,col.high)), 
          cluster_rows = isClusterByRow, 
          cluster_columns = isClusterByCol,
          show_row_dend = show_row_dend,
          show_column_dend = show_column_dend,
          clustering_distance_rows = distance_row,
          clustering_distance_columns= distance_column,
          show_row_names = T,show_column_names = F,
          top_annotation = cluster.annotation,
          row_names_gp = gpar(fontsize = rowFont.size))
}

#' Plot heatmap for identified marker genes
#' @param  object The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  n.TopGenes The number of top differential genes. Default 5
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.5
#' @param  min.expFraction The minimum cutoff for the fraction of expressing cell frequency. Default 0.3
#' @param  isClusterByRow is logical. If TRUE, the clustering by row will be performed.
#' @param  col.scaled.min The minimum scaled value. Default -2.5
#' @param  col.scaled.max The maximum scaled value. Default 2.5
#' @param  col.low The low color gradient. Default blue
#' @param  col.mid The mid color gradient. Default black
#' @param  col.high The high color gradient. Default red
#' @param  rowFont.size The row font size. Default 6
#' @param  split.line.col The split line color. Default white
#' @param  split.line.type The split line type. Default 1
#' @param  split.line.lwd The split line lwd. Default 1
#' @export

plot_heatmap_for_marker_genes<-function(object,cluster.type=c("walktrap","louvain","kmeans","merged_louvain","merged_walktrap","merged_kmeans"),
                                        n.TopGenes=5,min.log2FC=0.5,min.expFraction=0.3,isClusterByRow=F,col.scaled.min = -2.5,col.scaled.max = 2.5,
		col.low="blue",col.mid = "black",col.high = "red",rowFont.size=6,split.line.col="white", split.line.type=1,split.line.lwd=1){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
  
	clusters.info<-get_clusters(object)
	cluster.type <- match.arg(cluster.type)
	####################################
	if(cluster.type=="walktrap"){
		clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
	}else if(cluster.type=="louvain"){
		clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
	}else if(cluster.type=="kmeans"){
		clusters.info<-get_knn_graph.kmeans.cluster(object)
	}else if(cluster.type=="merged_walktrap"){
	  clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
	}else if(cluster.type=="merged_louvain"){
	  clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
	}else if(cluster.type=="merged_kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
	}else{
		stop("Need a clustering-method name!")
	}
	clusters.info$id<-0
	clusters.info$id<-as.numeric(gsub("cl","",clusters.info[,2]))
	
	clusters.ids<-table(clusters.info[,4])
	TopGenes.info<-data.frame()
	for(k in 1:length(clusters.ids)){
		my.cluster<-paste("cl",names(clusters.ids)[k],sep="")
		##################
		sig.genes<-getMarkerGenesInfo(object,cluster.type = cluster.type,cluster_id = my.cluster)
		if(nrow(sig.genes) >0){
		sig.genes<-subset(sig.genes,fishers.pval< 0.05 
						& wilcoxon.pval < 0.05 
						& ExpFractionA > ExpFractionB 
						& ExpFractionA > min.expFraction 
						& log2FC >= min.log2FC)
		}
		if(nrow(sig.genes) > 0){
			sig.genes<-sig.genes[order(sig.genes$FoldChange,decreasing = T),]	
			top.genes<-head(sig.genes,n=n.TopGenes)
			top.genes$cluster<-my.cluster
			TopGenes.info<-rbind(TopGenes.info,top.genes)
		}else{
			message<-paste("There is no gene passed the cut-off for:",my.cluster,"!.",sept="")
			print(message)
		}
	}
	topGenes<-unique(as.character(TopGenes.info$Gene))
	umi.dat<-get_normalized_umi(object)
	###############
	clusters.info.ranked<-clusters.info[order(as.numeric(clusters.info$id)),]
	umi.dat<-umi.dat[as.character(topGenes),as.character(clusters.info.ranked$Cell)]
	umi.dat<-log1p(umi.dat)
	###############
	mat_scaled = t(apply(umi.dat, 1, scale))
	df<-data.frame(cluster = as.character(clusters.info.ranked[,2]))
	###get color###
	clusters.info.part<-clusters.info.ranked[,c(2:3)]
	clusters.info.part<-unique(clusters.info.part)
	my.cols<-clusters.info.part[,2]
	names(my.cols)<-clusters.info.part[,1]
	cluster.cols<-list(cluster=my.cols)
	###############
	cluster.annotation = HeatmapAnnotation(df = df,col=cluster.cols)
	###############
	ht1<-Heatmap(mat_scaled, name ="expression",col = colorRamp2(c(col.scaled.min, 0, col.scaled.max), c(col.low,col.mid,col.high)), 
			cluster_rows = isClusterByRow, cluster_columns = F,show_row_dend = F,clustering_distance_rows = "euclidean",show_row_names = T,show_column_names = F,
			top_annotation = cluster.annotation,row_names_gp = gpar(fontsize = rowFont.size))
	
	draw(ht1)
	loc<-cumsum(as.numeric(clusters.ids))
	decorate_heatmap_body("expression", {
				for(l in 1:length(loc)){
					i = loc[l]
					x = i/ncol(mat_scaled)
					grid.lines(c(x, x), c(0, 1), gp = gpar(col = split.line.col, lty = split.line.type,lwd=split.line.lwd))
				}
			})
}

#' Plot KNN-Graph trajectory in 3D with the gene expression of a signature gene set
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of gene names.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  vertext.size The node size. Default 0.4
#' @param  point.color1 The low color gradient of gene expression. Default black
#' @param  point.color2 The high color gradient of gene expression. Default red
#' @export

plot_3D_knn_graph_label_by_a_signature_gene_set<-function(object,gene_list=c(),isNormalizedByHouseKeeping=F,
                                                          vertext.size=0.4,point.color1="black",point.color2="red"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  #####
  s.layout<-get_knn_graph.layout(object)
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(rownames(s.layout))]
  #####
  my.dat<-as.data.frame(s.layout)
  colnames(my.dat)<-c("Dim1","Dim2","Dim3")
  #####
  exprs<-as.matrix(umi.used[rownames(umi.used) %in% gene_list,])
  cells.info<-get_cells_annotation(object)
  rownames(cells.info)<-cells.info$Cell
  cells.info<-cells.info[as.character(rownames(s.layout)),]
  #####
  if(isNormalizedByHouseKeeping==T){
    genes.exp.sum<-Matrix::rowSums(umi.used)
    genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
    hk.genes<-names(genes.exp.sum[1:100])
    hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
    #####
    MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
    MM_score<-MM
  }else{
    MM <- Matrix::colSums(exprs)/(cells.info$UMI_count)
    MM_score<-MM
  }
  #####
  Genes.score<-data.frame(Genes_Score=MM_score)
  my.dat<-cbind(my.dat,Genes.score)
  ####################################
  p<-ggplot(my.dat,aes(Dim1,Dim2, color=Genes_Score)) + geom_point() + scale_colour_gradientn(colours = c("gray85",point.color1,point.color2),
                                                                                              values=c(0,0.01,1))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  gradient_score_color<-ggplot_build(p)$data[[1]]$colour
  #######################################
  g <- set_vertex_attr(get_knn_graph.graph(object), "color", value=gradient_score_color)
  graphjs(g,vertex.size=vertext.size,edge.color="gray45",
          layout=s.layout,
          fpl=300)
}

#' Plot KNN-Graph trajectory in 3D with the gene expression of selected genes
#' @param  object The SingCellaR object.
#' @param  show.gene The vector of gene names.
#' @param  vertext.size The node size. Default 0.4
#' @param  point.color1 The low color gradient of gene expression. Default black
#' @param  point.color2 The high color gradient of gene expression. Default red
#' @param  edge.color The edge color. Default gray45
#' @export

plot_3D_knn_graph_label_by_gene<-function(object,show.gene=c(),vertext.size=0.4,point.color1="black",point.color2="red",edge.color="gray45"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(show.gene)==0){
    stop("Required a gene!")
  }
  if(length(show.gene) > 1){
    stop("Please input only one gene name!")
  }
  #####
  s.layout<-get_knn_graph.layout(object)
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(rownames(s.layout))]
  #####
  my.dat<-as.data.frame(s.layout)
  colnames(my.dat)<-c("Dim1","Dim2","Dim3")
  #####
  exprs<-as.matrix(umi.used[rownames(umi.used) %in% show.gene,])
  #####
  my.dat<-cbind(my.dat,log1p(exprs))
  colnames(my.dat)<-c("Dim1","Dim2","Dim3",show.gene)
  ####################################
  p<-ggplot(my.dat,aes(Dim1,Dim2)) + geom_point(aes_string(colour = show.gene)) + scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))
  gradient_score_color<-ggplot_build(p)$data[[1]]$colour
  #######################################
  g <- set_vertex_attr(get_knn_graph.graph(object), "color", value=gradient_score_color)
  graphjs(g,vertex.size=vertext.size,edge.color=edge.color,
          layout=s.layout,fpl=300)
}

#' Plot KNN-Graph trajectory in 3D with the cluster name
#' @param  object The SingCellaR object.
#' @param  show_method The vector of clustering method names.
#' @param  vertext.size The node size. Default 0.4
#' @param  edge.color The edge color. Default gray45
#' @export

plot_3D_knn_graph_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                              vertext.size=0.4,edge.color="gray45"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	clusters.info<-get_clusters(object)
	show_method <- match.arg(show_method)
	####################################
	if(show_method=="walktrap"){
	  clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="louvain"){
	  clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_walktrap"){
	  clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_louvain"){
	  clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else{
	  stop("Need community detection or clustering-method names!")
	}
	rownames(clusters.info)<-clusters.info$Cell
	
	s.layout<-get_knn_graph.layout(object)
	clusters.info.ordered<-clusters.info[as.character(rownames(s.layout)),]
	
	my.cols<-clusters.info.ordered[,3]
	my.cluster<-clusters.info.ordered[,2]
	
	g <- set_vertex_attr(get_knn_graph.graph(object), "color", value=my.cols)
	
	graphjs(g,vertex.size=vertext.size,edge.color=edge.color,vertex.label=my.cluster,
			layout=s.layout,fpl=300)
}


#' Plot force-directed graph with the QC information
#' @param  object The SingCellaR object.
#' @param  vertext.size The node size. Default 0.5
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray45
#' @param  vertex.base.color The node base color. Default gray35
#' @export

plot_forceDirectedGraph_label_by_qc<-function(object,vertex.size=0.5,edge.size=0.2,edge.color="gray",vertex.base.color="gray35"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  #####
  cell.info<-get_cells_annotation(object)
  cell.info<-subset(cell.info,IsPassed==T)
  rownames(cell.info)<-cell.info$Cell
  
  fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
  colnames(fa2.layout) <- c("X1","X2")
  #########
  cell.info<-cell.info[rownames(fa2.layout),c("sampleID","UMI_count","detectedGenesPerCell","percent_mito")]
  fa2.layout<-cbind(fa2.layout,cell.info)
  #########
  rownames(fa2.layout) <- 1:nrow(fa2.layout)
  fa2.layout.x<-fa2.layout[,c(1:2)]
  edgelist <- get.edgelist(get_igraph.graph(object))
  edges <- data.frame(fa2.layout.x[edgelist[,1],], fa2.layout.x[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  
  p.x <- list()
  p.x[[1]] <- ggplot(fa2.layout,aes(X1,X2))+
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
      geom_point(aes(colour = UMI_count),size=vertex.size)+scale_colour_gradient()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
 
  p.x[[2]] <- ggplot(fa2.layout,aes(X1,X2)) + 
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
    geom_point(aes(colour = detectedGenesPerCell),size=vertex.size) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  
  p.x[[3]] <- ggplot(fa2.layout,aes(X1,X2)) + 
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
    geom_point(aes(colour = percent_mito),size=vertex.size) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  
  p.x[[4]] <- ggplot(fa2.layout,aes(X1,X2)) + 
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
    geom_point(aes(colour = sampleID),size=vertex.size) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  
  do.call(grid.arrange,p.x)
}

#' Plot force-directed graph with the sampleID
#' @param  object The SingCellaR object.
#' @param  vertext.size The node size. Default 0.5
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray45
#' @param  vertex.base.color The node base color. Default gray35
#' @export

plot_forceDirectedGraph_label_by_sampleID<-function(object,vertex.size=0.5,edge.size=0.2,edge.color="gray",vertex.base.color="gray35"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  #####
  cell.info<-get_cells_annotation(object)
  cell.info<-subset(cell.info,IsPassed==T)
  rownames(cell.info)<-cell.info$Cell
  
  fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
  colnames(fa2.layout) <- c("X1","X2")
  #########
  cell.info<-cell.info[rownames(fa2.layout),c("sampleID","UMI_count","detectedGenesPerCell","percent_mito")]
  fa2.layout<-cbind(fa2.layout,cell.info)
  #########
  rownames(fa2.layout) <- 1:nrow(fa2.layout)
  fa2.layout.x<-fa2.layout[,c(1:2)]
  edgelist <- get.edgelist(get_igraph.graph(object))
  edges <- data.frame(fa2.layout.x[edgelist[,1],], fa2.layout.x[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  
  p.x <- ggplot(fa2.layout,aes(X1,X2)) + 
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
    geom_point(aes(colour = sampleID),size=vertex.size) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  p.x
}

#' Plot force-directed graph with selected sampleIDs
#' @param  object The SingCellaR object.
#' @param  selected.sampleID The vector of selected sampleIDs.
#' @param  title The title of the plot.
#' @param  vertext.size The node size. Default 0.5
#' @param  vertex.front.color The color of selected samples. Default red
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default white
#' @param  vertex.base.color The node base color. Default gray
#' @param  alpha.base The alpha parameter of node. Default 0.2
#' @export

plot_forceDirectedGraph_label_by_selected_sampleID<-function(object,selected.sampleID=c(),title="",vertex.size=0.5,
                                                               vertex.front.color="red",edge.size=0.2,edge.color="white",
                                                               vertex.base.color="gray",alpha.base=0.2){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(selected.sampleID)==0){
    stop("Please input selected sample IDs!")
  }
  #####
  cell.info<-get_cells_annotation(object)
  cell.info<-subset(cell.info,IsPassed==T)
  rownames(cell.info)<-cell.info$Cell
  
  fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
  colnames(fa2.layout) <- c("X1","X2")
  #########
  cell.info<-cell.info[rownames(fa2.layout),c("sampleID","UMI_count","detectedGenesPerCell","percent_mito")]
  fa2.layout<-cbind(fa2.layout,cell.info)
  fa2.layout$color<-vertex.base.color
  sample.index<-which(fa2.layout$sampleID %in% selected.sampleID)
  fa2.layout$color[sample.index]<-vertex.front.color
  
  fa2.layout$alpha<-alpha.base
  sample.index<-which(fa2.layout$sampleID %in% selected.sampleID)
  fa2.layout$alpha[sample.index]<-1
  
  #########
  rownames(fa2.layout) <- 1:nrow(fa2.layout)
  fa2.layout.x<-fa2.layout[,c(1:2)]
  edgelist <- get.edgelist(get_igraph.graph(object))
  edges <- data.frame(fa2.layout.x[edgelist[,1],], fa2.layout.x[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  
  p.x <- ggplot(fa2.layout,aes(X1,X2)) + ggtitle(title)+
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
    geom_point(color=fa2.layout$color,size=vertex.size,alpha=fa2.layout$alpha) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
  p.x <-p.x+theme(plot.title = element_text(color="red", size=16, face="bold.italic",hjust = 0.5))
  p.x
}

#' Plot force-directed graph with a feature of interest
#' @param  object The SingCellaR object.
#' @param  feature The selected feature of interest.
#' @param  title The title of the plot.
#' @param  vertext.size The node size. Default 0.5
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray85
#' @param  background.color The background color. Default white
#' @export

plot_forceDirectedGraph_label_by_a_feature_of_interest<-function(object,feature="",title="",vertex.size=0.5,
                                                                 edge.size=0.2,edge.color="gray85",background.color="white"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(feature==""){
    stop("Please input a feature to plot!")
  }
  if(feature %in% colnames(object@meta.data)==T){
    cell.info<-get_cells_annotation(object)
    cell.info<-subset(cell.info,IsPassed==T)
    rownames(cell.info)<-cell.info$Cell
    
    fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
    colnames(fa2.layout) <- c("X1","X2")
    #########
    cell.info<-cell.info[rownames(fa2.layout), colnames(object@meta.data)]
    fa2.layout<-cbind(fa2.layout,cell.info)
    #########
    rownames(fa2.layout) <- 1:nrow(fa2.layout)
    fa2.layout.x<-fa2.layout[,c(1:2)]
    edgelist <- get.edgelist(get_igraph.graph(object))
    edges <- data.frame(fa2.layout.x[edgelist[,1],], fa2.layout.x[edgelist[,2],])
    colnames(edges) <- c("X1","Y1","X2","Y2")
    
    p.x <- ggplot(fa2.layout,aes(X1,X2)) + ggtitle(title)+
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
      geom_point(aes_string(colour = feature),size=vertex.size) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())
    p.x <-p.x+theme(plot.title = element_text(color="red", size=16, face="bold.italic",hjust = 0.5))
    p.x <-p.x+theme(plot.background = element_rect(fill = background.color))
    p.x
    
  }else{
    error<-paste("Couldnot find :",feature," in the column names of cell meta data.",sep="")
    stop(error)
  }
}

#' Plot force-directed graph with gene expression
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of genes.
#' @param  vertext.size The node size. Default 0.5
#' @param  vertex.color1 The first color of gene expression gradient. Default blue
#' @param  vertext.color2 The second color of gene expression gradient. Default red
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray
#' @param  vertex.base.color The node base color. Default gray35
#' @param  background.color The background color. Default white
#' @export

plot_forceDirectedGraph_label_by_genes<-function(object,gene_list=c(),vertex.size=0.5,vertex.color1="blue",
		vertex.color2="red",edge.size=0.2,edge.color="gray",vertex.base.color="gray35",background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	#####
	fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
	colnames(fa2.layout) <- c("X1","X2")
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(rownames(fa2.layout))]
	#####
	rownames(fa2.layout) <- 1:nrow(fa2.layout)
	#####
	my.dat<-fa2.layout[,(1:2)]
	
	for(i in 1:length(gene_list)){
	  genes.x<-gene_list[i]
	  g.ind<-which(rownames(umi.used)==genes.x)
	  my.sub.umi<-as.data.frame(log1p(umi.used[g.ind,]))
	  colnames(my.sub.umi)<-genes.x
	  #################################
	  my.dat<-cbind(my.dat,my.sub.umi)
	}
	
	X.colnames<-colnames(my.dat)
	X.colnames<-gsub("-","_",X.colnames)
	X.colnames<-gsub("\\+","_",X.colnames)
	colnames(my.dat)<-X.colnames
	
	my.fg.names<-colnames(my.dat)
	my.fg.names<-my.fg.names[3:length(my.fg.names)]
	#########
	edgelist <- get.edgelist(get_igraph.graph(object))
	edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	
	p.x <- list()
	for(j in 1:length(my.fg.names)){
		genes.x<-my.fg.names[j]
		genes.x<-gsub("-","_",genes.x)
		genes.x<-gsub("\\+","_",genes.x)
		
		p.x[[j]] <- ggplot(my.dat,aes(X1,X2)) + 
				geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
				geom_point(aes_string(colour = genes.x),size=vertex.size) + 
				scale_colour_gradientn(colours = c(vertex.base.color,vertex.color1,vertex.color2),values=c(0,0.01,1))+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+
				theme(plot.background = element_rect(fill = background.color))
				
	}
	do.call(grid.arrange,p.x)
}

#' Plot force-directed graph with the expression of a gene set
#' @param  object The SingCellaR object.
#' @param  gene_list The vector of genes.
#' @param  gene_set_name The name of a gene set.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  vertext.size The node size. Default 0.5
#' @param  vertex.color1 The first color of gene expression gradient. Default blue
#' @param  vertext.color2 The second color of gene expression gradient. Default red
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray
#' @param  vertex.base.color The node base color. Default gray35
#' @param  background.color The background color. Default white
#' @export

plot_forceDirectedGraph_label_by_a_signature_gene_set<-function(object,gene_list=c(),gene_set_name=c(),
                                                                isNormalizedByHouseKeeping=F,vertex.size=0.5,
                                                                vertex.color1="blue",vertex.color2="red",
                                                                edge.size=0.2,edge.color="gray",
                                                                vertex.base.color="gray35",background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gene_list)==0){
		stop("Required list of genes!")
	}
	if(length(gene_list) < 2){
		stop("Required more than 2 genes!")
	}
	if(length(gene_set_name)==0){
		stop("Required a name of gene set!")
	}
	#####
	fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
	colnames(fa2.layout) <- c("X1","X2")
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(rownames(fa2.layout))]
	#####
	UMI_count<-Matrix::colSums(umi.used)
	rownames(fa2.layout) <- 1:nrow(fa2.layout)
	#####
	my.dat<-fa2.layout[,(1:2)]
	#####
	exprs<-as.matrix(umi.used[rownames(umi.used) %in% gene_list,])
	########get house keeping genes#####
	if(isNormalizedByHouseKeeping==T){
		genes.exp.sum<-Matrix::rowSums(umi.used)
		genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
		hk.genes<-names(genes.exp.sum[1:100])
		hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
		#####
		MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
		MM_score<-MM
	}else{
		MM <- Matrix::colSums(exprs)/(UMI_count)
		MM_score<-MM
	}
	####################################
	####################################
	edgelist <- get.edgelist(get_igraph.graph(object))
	edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	####################################
	Genes.score<-data.frame(Genes_score=MM_score)
	fa2.layout<-cbind(fa2.layout,Genes.score)
	####################################
	ggplot(fa2.layout,aes(X1,X2, color=Genes_score)) + 
			geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
			geom_point(size=vertex.size) + 
			scale_colour_gradientn(colours = c(vertex.base.color,vertex.color1,vertex.color2),values=c(0,0.01,1))+
			ggtitle(gene_set_name)+
			theme(plot.title = element_text(hjust = 0.5))+ 
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
			theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank(),axis.title.x=element_blank())+
			theme(plot.background = element_rect(fill = background.color))
	
}

#' Plot force-directed graph with the expression of selected multiple gene sets
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file for gene sets.
#' @param  show_gene_sets The names of selected gene sets.
#' @param  custom_color The custom color names.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  vertext.size The node size. Default 1.5
#' @param  showEdge is logical. If TRUE, the graph's edges will be displayed.
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray
#' @param  showLegend is logical. If TRUE, the figure legend will be displayed.
#' @param  background.color The background color. Default white
#' @export

plot_forceDirectedGraph_label_by_multiple_gene_sets<-function(object,gmt.file=c(),show_gene_sets=c(),custom_color=c(),
															 isNormalizedByHouseKeeping=T,vertex.size=1.5,showEdge=T,
															 edge.size=0.2,edge.color="gray",showLegend = T,background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gmt.file)==0){
		stop("Required the path to GMT file!")
	}
	if(length(show_gene_sets)==0){
		stop("Required gene set names!")
	}
	if(length(show_gene_sets)==1){
		stop("Required at least two gene sets")
	}
	if(length(show_gene_sets) > 10){
		stop("This is over limitation at 10 gene sets!")
	}
	if(length(custom_color) > 0){
		if(length(custom_color) !=length(show_gene_sets)){
			stop("The number of colors is not equal to the the number of input gene sets.")
		}
	}
	#####
	signature.sets<-get_gmtGeneSets(gmt.file)
	#####
	check.1<-match(show_gene_sets,names(signature.sets))
	check.2<-sum(length(which(is.na(check.1))))
	if(check.2 > 0){
		stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
	}
	fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
	colnames(fa2.layout) <- c("X1","X2")
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(rownames(fa2.layout))]
	#####
	UMI_count<-Matrix::colSums(umi.used)
	my.dat<-fa2.layout[,(1:2)]
	#####
	genes.exp.sum<-Matrix::rowSums(umi.used)
	genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
	hk.genes<-names(genes.exp.sum[1:100])
	hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
	#####
	Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
	for(i in 1:length(show_gene_sets)){
		set.name<-show_gene_sets[i]
		genes.x<-signature.sets[[set.name]]
		exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
		
		if(isNormalizedByHouseKeeping==T){
			MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
		}else{
			MM<-Matrix::colSums(exprs)/UMI_count
		}
		MM.f<-data.frame(score=MM)
		colnames(MM.f)<-set.name
		Genes.score<-cbind(Genes.score,MM.f)
	}
	Genes.score<-Genes.score[-c(1)]
	####Assign colors#####################
	if(length(custom_color)==0){
		multi.colors<-rainbow(ncol(Genes.score))
	}else{
		multi.colors<-custom_color
	}
	######################################
	my.colors.list<-list()
	for(j in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[j])
		x.color<-multi.colors[j]
		p1<-ggplot(my.dat,aes(X1,X2)) + geom_point(aes(colour = Genes.score[,j])) + scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
		My_color<-ggplot_build(p1)$data[[1]]$colour
		my.colors.list[[set.name]]<-My_color
	}
	my.alpha.list<-list()
	for(j in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[j])
		my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
		my.alpha.list[[set.name]]<-my.alpha
	}
	###################################
	rownames(fa2.layout) <- 1:nrow(fa2.layout)
	####################################
	edgelist <- get.edgelist(get_igraph.graph(object))
	edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	###################################
	if(showEdge==T){
		p<-ggplot(my.dat,aes(X1,X2))+
			geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)
	}else{
		p<-ggplot(my.dat,aes(X1,X2))
	}
	for(k in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[k])
		p<-p+geom_point(color=my.colors.list[[set.name]], size = vertex.size,alpha=my.alpha.list[[set.name]],stroke = 0)
	}
	
	p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
			theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
			theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
			theme(plot.background = element_rect(fill = background.color))
			
	
	if(showLegend==T){
		df<-data.frame(Name=show_gene_sets)
		df<-as.matrix(df)
	
		tp<-ggtexttable(df,rows = NULL,
			theme = ttheme(
					colnames.style = colnames_style(fill = "white"),
					tbody.style = tbody_style(fill = multi.colors)
			)
		)
		return(ggarrange(p,tp,ncol = 2, nrow = 1,widths = c(1,0.2)))
	}else{
		return(p)
	}
}

#' Plot force-directed graph with identified clusters
#' @param  object The SingCellaR object.
#' @param  show_method The clustering method names used to generate clusters.
#' @param  vertext.size The node size. Default 2
#' @param  mark.clusters is logical. If TRUE, the number of identified clusters will be displayed.
#' @param  mark.font.size The font size of the cluster name. Default 10
#' @param  mark.font.color The font color of the cluster name. Default yellow
#' @param  add.cluster.in.cells is logical. If TRUE, the cluster name will be displayed in each data point.
#' @param  cluster.font.size The font size of cluster name in each data point. Default 1
#' @param  cluster.font.color The font color of cluster name in each data point. Default yellow
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray
#' @param  background.color The background color. Default white
#' @export

plot_forceDirectedGraph_label_by_clusters<-function(object,show_method=c("walktrap","louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                                    vertex.size=2,mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",
                                                    add.cluster.in.cells=FALSE,cluster.font.size=1,cluster.font.color="yellow",
                                                    edge.size=0.2,edge.color="gray",background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	show_method <- match.arg(show_method)
	###
	fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
	clusters.info<-get_clusters(object)
	
	if(show_method=="walktrap"){
	  clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="louvain"){
	  clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_walktrap"){
	  clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_louvain"){
	  clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else if(show_method=="merged_kmeans"){
	  clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
	  colnames(clusters.info)<-c("Cell","cluster","cluster_color")
	}else{
	  stop("Need community detection or clustering-method names!")
	}
	colnames(fa2.layout)<-c("X1","X2")
	fa2.layout$Cell<-rownames(fa2.layout)
	
	res.fa2<-merge(fa2.layout,clusters.info)
	cluster.id<-gsub("cl","",res.fa2$cluster)
	####################################
	rownames(fa2.layout) <- 1:nrow(fa2.layout)
	fa2.layout<-fa2.layout[,c(1:2)]
	####################################
	edgelist <- get.edgelist(get_igraph.graph(object))
	edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
	colnames(edges) <- c("X1","Y1","X2","Y2")
	####################################
	if(add.cluster.in.cells==TRUE){
		p<- qplot(X1,X2, data=res.fa2,colour=cluster)+
				geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=vertex.size)+ 
				geom_text(aes(label = add.cluster.in.cells),size = cluster.font.size,colour=cluster.font.color)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
				theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
				theme(plot.background = element_rect(fill = background.color))
		
		if(mark.clusters==TRUE){
			edata <- res.fa2[,c("X1","X2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
		
	}else{
		p<- qplot(X1,X2, data=res.fa2,colour=cluster)+
				geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
				geom_point(aes(fill=cluster),colour="black",pch=21,size=vertex.size)+
				theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
				theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
				theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
				theme(plot.background = element_rect(fill = background.color))
		
		if(mark.clusters==TRUE){
			edata <- res.fa2[,c("X1","X2")]
			edata$cluster.id<-cluster.id
			colnames(edata) <- c('x', "y", "z")
			center <- aggregate(cbind(x,y) ~ z, data = edata, median)
			p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
		}
	}
	p<-p+ggtitle(show_method)
	return(p)
}

#' Plot jaccard similar across indetified clusters
#' @param  object The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.3
#' @param  min.expFraction The minimum fraction of expressing cell quency. Default 0.3
#' @export

plot_jaccard_similarity_among_clusters <- function(object,cluster.type=c("louvain","walktrap","kmeans"),
                                                   min.log2FC=0.30,min.expFraction=0.30){
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  if(cluster.type=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster")]
  }else if(cluster.type=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster")]
  }else if(cluster.type=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)
  }else{
    stop("Need a clustering-method name!")
  }
  clusters.ids<-paste("cl",1:length(unique(clusters.info[,2])),sep="")
  #######################################
  clusters.combn<-t(combn(clusters.ids,2))
  clusters.combn.ids<-t(combn(1:length(unique(clusters.info[,2])),2))
  #######################################
  jaccard.mat <- matrix(nrow=length(clusters.ids),ncol=length(clusters.ids)) 
  for(k in 1:nrow(clusters.combn)){
    cl.x1<-clusters.combn[k,1]
    cl.x2<-clusters.combn[k,2]
    id.x1<-clusters.combn.ids[k,1]
    id.x2<-clusters.combn.ids[k,2]
    markers.x1<-getMarkerGenesInfo(object,cluster.type=cluster.type,cluster_id=cl.x1)
    markers.x2<-getMarkerGenesInfo(object,cluster.type=cluster.type,cluster_id=cl.x2)
    
    if(nrow(markers.x1) > 0 & nrow(markers.x2) > 0){
      markers.x1<-subset(markers.x1,fishers.pval< 0.05 
                       & wilcoxon.pval < 0.05 
                       & ExpFractionA > ExpFractionB
                       & ExpFractionA > min.expFraction 
                       & log2FC >= min.log2FC)
    
    
    markers.x2<-subset(markers.x2,fishers.pval< 0.05 
                       & wilcoxon.pval < 0.05 
                       & ExpFractionA > ExpFractionB
                       & ExpFractionA > min.expFraction 
                       & log2FC >= min.log2FC)
    
    a<-intersect(markers.x1$Gene,markers.x2$Gene)
    b<-union(markers.x1$Gene,markers.x2$Gene)
    my.jaccard.index<-(length(a)/length(b))
    jaccard.mat[id.x1,id.x2]<-my.jaccard.index
    }else{
      jaccard.mat[id.x1,id.x2]<-0
    }
  }
  rownames(jaccard.mat)<-clusters.ids
  colnames(jaccard.mat)<-clusters.ids
  diag(jaccard.mat)<-1
  jaccard.mat[lower.tri(jaccard.mat)] <- t(jaccard.mat)[lower.tri(t(jaccard.mat))]
  ##########pheatmap#######################
  pheatmap(jaccard.mat,display_numbers = T)
  ##########################################
}

#' Plot the shortest-path on the force-directed graph with identified clusters
#' @param  object The SingCellaR object.
#' @param  start_cluster The starting cluster.
#' @param  end_cluster  The ending cluster.
#' @param  show_method The clustering method names used to generate clusters.
#' @param  vertex.size The node size. Default 2
#' @param  path.size The path size. Default 2
#' @param  path.type The path type. Default 1
#' @param  path.color The path color. Default red
#' @param  mark.clusters is logical. If TRUE, the number of identified clusters will be displayed.
#' @param  mark.font.size The font size of the cluster name. Default 10
#' @param  mark.font.color The font color of the cluster name. Default yellow
#' @param  add.cluster.in.cells is logical. If TRUE, the cluster name will be displayed in each data point.
#' @param  cluster.font.size The font size of cluster name in each data point. Default 1
#' @param  cluster.font.color The font color of cluster name in each data point. Default yellow
#' @param  edge.size  The size of the edge. Default 0.2
#' @param  edge.color The edge color. Default gray
#' @export

plot_ShortestPath_on_forceDirectedGraph_label_by_clusters<-function(object,start_cluster = "",end_cluster = "",
                                                                    show_method=c("walktrap","louvain","kmeans"),
                                                                    vertex.size=2,path.size=2,path.type=1,path.color="red",
                                                                    mark.clusters=TRUE,mark.font.size=10,mark.font.color="yellow",
                                                                    add.cluster.in.cells=FALSE,cluster.font.size=1,
                                                                    cluster.font.color="yellow",edge.size=0.2,edge.color="gray"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  show_method <- match.arg(show_method)
  
  if(start_cluster==""){
    stop("Please input a cluster!")
  }
  if(end_cluster==""){
    stop("Please input a cluster!")
  }
  ####################################
  fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
  clusters.info<-get_clusters(object)
  
  if(show_method=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else if(show_method=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)
    colnames(clusters.info)<-c("Cell","cluster","cluster_color")
  }else{
    stop("Need community detection or clustering-method names!")
  }
  colnames(fa2.layout)<-c("X1","X2")
  fa2.layout$Cell<-rownames(fa2.layout)
  
  res.fa2<-merge(fa2.layout,clusters.info)
  cluster.id<-gsub("cl","",res.fa2$cluster)
  ####################################
  rownames(fa2.layout) <- 1:nrow(fa2.layout)
  fa2.layout<-fa2.layout[,c(1:2)]
  ####################################
  edgelist <- get.edgelist(get_igraph.graph(object))
  edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  ####################################
  st_path<-findShortestPath(object,start_cluster = start_cluster,end_cluster = end_cluster,cluster.type = show_method)
  st_res<-res.fa2[st_path$members,]
  ####################################
  p<- qplot(X1,X2, data=res.fa2,colour=cluster)+
    geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)+
    geom_point(aes(fill=cluster),colour="black",pch=21,size=vertex.size)+
    geom_line(aes(x=X1,y=X2),data=st_res,size=path.size,color=path.color,linetype = path.type)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())
  if(mark.clusters==TRUE){
    edata <- res.fa2[,c("X1","X2")]
    edata$cluster.id<-cluster.id
    colnames(edata) <- c('x', "y", "z")
    center <- aggregate(cbind(x,y) ~ z, data = edata, median)
    p<-p+geom_text(data=center, aes_string(x = "x", y = "y", label = "z"), size = mark.font.size, colour = mark.font.color)
  }
  p<-p+ggtitle(show_method)
  return(p)
}

#' Plot the shortest-path on the 3D KNN-Graph
#' @param  object The SingCellaR object.
#' @param  start_cluster The starting cluster.
#' @param  end_cluster  The ending cluster.
#' @param  cluster.type The clustering method name.
#' @param  path.color The path color. Default red
#' @param  vertext.size The node size. Default 0.4
#' @param  edge.color The edge color. Default gray
#' @export

plot_ShortestPath_on_3D_knn_graph<-function(object,start_cluster="",end_cluster="",cluster.type=c("louvain","walktrap","kmeans"),
                                            path.color="red",vertext.size=0.4,edge.color="gray45"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  cluster.type <- match.arg(cluster.type)
  if(start_cluster==""){
    stop("Please input a cluster!")
  }
  if(end_cluster==""){
    stop("Please input a cluster!")
  }
  clusters.info<-get_clusters(object)
  ####################################
  if(cluster.type=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
  }else if(cluster.type=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
  }else if(cluster.type=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)
  }else{
    stop("Need a clustering-method name!")
  }
  rownames(clusters.info)<-clusters.info$Cell
  
  s.layout<-get_knn_graph.layout(object)
  clusters.info.ordered<-clusters.info[as.character(rownames(s.layout)),]
  ##################################
  st_path<-findShortestPath(object,start_cluster = start_cluster,end_cluster = end_cluster,cluster.type = cluster.type)
  
  clusters.info.ordered$path.cols<-"gray85"
  clusters.info.ordered$path.cols[st_path$members]<-path.color
  g <- set_vertex_attr(get_knn_graph.graph(object), "color", value=clusters.info.ordered$path.cols)
  graphjs(g,vertex.size=vertext.size,edge.color=edge.color,layout=s.layout,fpl=300)
}

#' Plot TSNE with AUCell score
#' @param  object The SingCellaR object.
#' @param  geneSets.AUC The dataframe of AUCell score per gene set.
#' @param  show_gene_sets The vector of gene set names
#' @param  isUseCutOffScore is logical. If TRUE, the cutoff score will be applied.
#' @param  cutOffScore The cutoff score. Default 0
#' @param  isSubtract.random.gene is logical. If TRUE, the AUC score will be subtracted by the AUCell score from the random gene set.
#' @param  point.size The point size. Default 1
#' @param  point.color1 The color for a low AUCell score. Default gray
#' @param  point.color2 The color for a high AUCell score . Default red
#' @export

plot_tsne_label_by_AUCell_score<-function(object,geneSets.AUC,show_gene_sets=c(),isUseCutOffScore=F,
                                          cutOffScore=0,isSubtract.random.gene=T,point.size=1,
                                          point.color1="gray",
                                          point.color2="red"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  if(nrow(geneSets.AUC)==0){
    stop("Please input the data frame object of the AUC gene set score!")
  }
  res.tsne<-get_tsne.result(object)
  my.dat<-res.tsne[,c("Cell","TSNE1","TSNE2")]
  #####
  geneSets.AUC.updated<-geneSets.AUC[as.character(my.dat$Cell),]
  Genes.score<-geneSets.AUC.updated
  
  if(isSubtract.random.gene==T){
    Genes.score<-Genes.score-Genes.score$random_gene
    Genes.score<-Genes.score[show_gene_sets]
  }else{
    Genes.score<-Genes.score[show_gene_sets]
  }
  if(isUseCutOffScore==T){
    Genes.score[Genes.score < cutOffScore]<-0
    gg.title<-paste("Cells with AUCell >",cutOffScore)
  }else{
    gg.title<-""
  }
  my.dat<-cbind(my.dat,Genes.score)
  p.x <- list()
  for(j in 1:length(show_gene_sets)){
    genes.x<-show_gene_sets[j]
    p.x[[j]] <- ggplot(my.dat,aes(TSNE1,TSNE2)) + geom_point(aes_string(colour = genes.x),size=point.size) + 
      scale_colour_gradientn(colours = c("gray85",point.color1,point.color2),values=c(0,0.01,1))+
      ggtitle(gg.title)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
  }
  do.call(grid.arrange,p.x)
}

#' Plot heatmap for GSEA for all comparisons
#' @param  object The SingCellaR object.
#' @param  fGSEA_results.data.frame The input dataframe for fGSEA analysis results.
#' @param  isApplyCutoff is logical. If TRUE, the cutoff will be applied.
#' @param  adjusted_pval The cutoff adjusted p-value. Default 0.25
#' @param  use_pvalues_for_clustering is logical. If TRUE, the p-values will be used for clustering analysis.
#' @param  show_NES_score is logical. If TRUE, the enrichment score will be displayed.
#' @param  show_only_NES_positive_score is logical. If TRUE, the only NES positive scores will be displayed.
#' @param  show_only_NES_negative_score is logical. If TRUE, the only NES negative scores will be displayed.
#' @param  custom_order_samples The vector of input custom order of GSEA comparison on the heatmap.
#' @param  fontsize_row The row font size. Default 5
#' @param  format.digits The format of digit number shows on the heatmap. Default 2
#' @param  cluster_rows is logical. If TRUE, the clustering analysis will be performed by row.
#' @param  cluster_cols is logical. If TRUE, the clustering analysis will be performed by column.
#' @param  clustering_method The HCL clustering method. Default complete
#' @param  clustering_distance_rows The distance of rows for the clustering. Default euclidean
#' @param  clustering_distance_cols The distance of columns for the clustering. Default euclidean.
#' @param  gaps_row The vector of gaps for row.
#' @param  gaps_col The vector of gaps for column.
#' @param  show_text_for_ns is logical. If TRUE. The ns or 'non-significant' will be displayed. 
#' @export

plot_heatmap_for_fGSEA_all_clusters <- function(fGSEA_results.data.frame,isApplyCutoff=F,adjusted_pval=0.25,
                                                add_small_value=0.0001,use_pvalues_for_clustering=T,
                                                show_NES_score=T,show_only_NES_positive_score=F,
                                                show_only_NES_negative_score=F,custom_order_samples=c(),
                                                fontsize_row=5,format.digits=2,cluster_rows = TRUE,
                                                cluster_cols = TRUE,clustering_method="complete",
                                                clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",
                                                gaps_row = NULL, gaps_col = NULL,show_text_for_ns=T){
  
  if(nrow(fGSEA_results.data.frame)==0){
    stop("Requires fGSEA results in the data frame object!")
  }
  
  if(isApplyCutoff==T){
    fGSEA_results.data.frame<-subset(fGSEA_results.data.frame,padj < adjusted_pval)
  }
  
  if(show_only_NES_positive_score==T & show_only_NES_negative_score==T){
    stop("Please set show_only_NES_positive_score or show_only_NES_negative_score to false.")
  }
  
  if(nrow(fGSEA_results.data.frame)==0){
    stop("There is no pathway that passed this adjust p-value.")
  }
  
  FM<- as.data.frame(dcast(fGSEA_results.data.frame, pathway~cluster,value.var="padj",na.rm=T))
  rownames(FM)<-FM$pathway
  FM<-FM[-c(1)]
  FM[is.na(FM)==T]<-1
  FM <- -log10(FM+add_small_value)
  
  if(show_NES_score==T){
    
    FM2<- as.data.frame(dcast(fGSEA_results.data.frame, pathway~cluster,value.var="NES",na.rm=T))
    rownames(FM2)<-FM2$pathway
    FM2<-FM2[-c(1)]
    FM2[is.na(FM2)==T]<-0
    
    if(show_only_NES_positive_score==T){
      FM2[FM2 < 0]<-0
      FM2<-FM2[rowSums(FM2) > 0,]
      FM<-FM[rownames(FM2),]
      FM[FM2 == 0]<-0
      
      info.y<-as.data.frame(FM2)
      info.y<-format(info.y,digits = format.digits,na.encode = TRUE)
      if(show_text_for_ns==T){
     	 info.y[info.y=="0.0" | info.y=="0.00" | info.y=="0.000"]<-"ns"
      }else{
      	 info.y[info.y=="0.0" | info.y=="0.00" | info.y=="0.000"]<-""
      }

    }else if(show_only_NES_negative_score==T){
      FM2[FM2 > 0]<-0
      FM2<-FM2[rowSums(FM2) < 0,]
      FM<-FM[rownames(FM2),]
      FM[FM2 == 0]<-0
      
      info.y<-as.data.frame(FM2)
      info.y<-format(info.y,digits = format.digits)
      if(show_text_for_ns==T){
      	info.y[info.y==" 0.0" | info.y==" 0.00" | info.y==" 0.000"]<-"ns"
       }else{
       	info.y[info.y==" 0.0" | info.y==" 0.00" | info.y==" 0.000"]<-""
       }
    }else{
      
      info.y<-as.matrix(FM2)
      info.y<-format(info.y,digits = format.digits)
      if(show_text_for_ns==T){
       info.y[info.y==" 0.0" | info.y==" 0.00" | info.y==" 0.000"]<-"ns"
      }else{
       info.y[info.y==" 0.0" | info.y==" 0.00" | info.y==" 0.000"]<-""
      }
    }
    if(length(custom_order_samples) > 0){
        FM<-FM[,custom_order_samples]
        FM2<-FM2[,custom_order_samples]
        info.y<-info.y[,custom_order_samples]
    }
    if(use_pvalues_for_clustering==TRUE){
             pheatmap(FM,clustering_method=clustering_method,clustering_distance_rows=clustering_distance_rows,
             clustering_distance_cols=clustering_distance_cols,cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,color=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
             display_numbers=info.y,fontsize_row=fontsize_row,gaps_row = gaps_row, gaps_col = gaps_col)
    }else{
            pheatmap(FM2,clustering_method=clustering_method,clustering_distance_rows=clustering_distance_rows,
                 clustering_distance_cols=clustering_distance_cols,cluster_rows = cluster_rows,
                 cluster_cols = cluster_cols,color=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
                 display_numbers=info.y,fontsize_row=fontsize_row,gaps_row = gaps_row, gaps_col = gaps_col)
    }
    
  }else{
    
    FM2<- as.data.frame(dcast(fGSEA_results.data.frame, pathway~cluster,value.var="NES",na.rm=T))
    rownames(FM2)<-FM2$pathway
    FM2<-FM2[-c(1)]
    FM2[is.na(FM2)==T]<-0
    
    if(length(custom_order_samples) > 0){
      FM<-FM[,custom_order_samples]
      FM2<-FM2[,custom_order_samples]
    }
    if(use_pvalues_for_clustering==TRUE){
      pheatmap(FM,clustering_method=clustering_method,clustering_distance_rows=clustering_distance_rows,
             clustering_distance_cols=clustering_distance_cols,cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,color=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
             fontsize_row=fontsize_row,gaps_row = gaps_row, gaps_col = gaps_col)
    }else{
      pheatmap(FM2,clustering_method=clustering_method,clustering_distance_rows=clustering_distance_rows,
               clustering_distance_cols=clustering_distance_cols,cluster_rows = cluster_rows,
               cluster_cols = cluster_cols,color=colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
               fontsize_row=fontsize_row,gaps_row = gaps_row, gaps_col = gaps_col)
    }
  }
}

#' Plot UMAP with AUCell score
#' @param  object The SingCellaR object.
#' @param  AUCell_gene_set_name The vector of gene set names.
#' @param  AUCell_score The dataframe of AUCell score per gene set.
#' @param  AUCell_cutoff The cutoff score. Default 0
#' @param  IsShowOnlySampleIDs is logical. If TRUE, only selected sample IDs will be shown.
#' @param  selected.sampleID The selected sample IDs.
#' @param  IsLimitedAUCscoreByClusters is logical. If TRUE, AUCell score will be limited by selected clusters.
#' @param  selected.limited.clusters The selected clusters of interest.
#' @param  clustering_method The clustering method.
#' @param  point.size The point size. Default 1
#' @param  point.color1 The color for a low AUCell score. Default gray
#' @param  point.color2 The color for a high AUCell score . Default red
#' @export

plot_umap_label_by_AUCell_score<-function(object,AUCell_gene_set_name=c(),AUCell_score,AUCell_cutoff=0,
                                          IsShowOnlySampleIDs=FALSE,selected.sampleID=c(),IsLimitedAUCscoreByClusters=FALSE,
                                          selected.limited.clusters=c(),clustering_method="louvain",
                                          point.size=1,point.color1="black",point.color2="red"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(AUCell_gene_set_name)==0){
    stop("Required a name of gene set!")
  }
  #####
  res.umap<-get_umap.result(object)
  #####
  my.AUCell<-AUCell_score[as.character(res.umap$Cell),]
  my.AUCell.score<-my.AUCell[,which(colnames(my.AUCell)==AUCell_gene_set_name)]
  #####
  Genes.score<-data.frame(AUCell_Score=my.AUCell.score)
  f.umap<-cbind(res.umap,Genes.score)
  plot.umap<-f.umap
  plot.umap$AUCell_Score[plot.umap$AUCell_Score < AUCell_cutoff]<-0
  AUCell_gene_set_name<-paste(AUCell_gene_set_name," (Cells with AUC >",AUCell_cutoff, ")",sep="")
  
  if(IsShowOnlySampleIDs==T & IsLimitedAUCscoreByClusters==F){
    sample.index<-which(plot.umap$sampleID %in% selected.sampleID)
    plot.umap<-plot.umap[sample.index,]
    names.sampleIDs<-paste(names(table(plot.umap$sampleID)), collapse=",") 
    AUCell_gene_set_name<-paste(AUCell_gene_set_name," [show only: ",names.sampleIDs,"]",sep="")
    
  }else if(IsShowOnlySampleIDs==F & IsLimitedAUCscoreByClusters==T){
    
    clusters.info<-get_clusters(object)
    
    if(clustering_method=="walktrap"){
      clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="louvain"){
      clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_walktrap"){
      clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_louvain"){
      clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else{
      stop("Need community detection or clustering-method names!")
    }
    
    plot.umap<-merge(plot.umap,clusters.info)
    sample.index<-which(plot.umap$cluster %in% selected.limited.clusters)
    edited_score<-rep(0,nrow(plot.umap))
    edited_score[sample.index]<-plot.umap$AUCell_Score[sample.index]
    plot.umap$AUCell_Score<-edited_score
    total.n.cells<-nrow(plot.umap)
    #################
    z<-plot.umap[sample.index,]
    z<-subset(z,AUCell_Score > AUCell_cutoff)
    AUC.n.cells<-nrow(z)
    #################
    txt<-paste("total cells=",total.n.cells," ; cells with AUCell_score=",AUC.n.cells,sep="")
    AUCell_gene_set_name<-paste(AUCell_gene_set_name,txt,sep="\n")
    
  }else if(IsShowOnlySampleIDs==T & IsLimitedAUCscoreByClusters==T){
    sample.index<-which(plot.umap$sampleID %in% selected.sampleID)
    plot.umap<-plot.umap[sample.index,]
    ##########
    names.sampleIDs<-paste(names(table(plot.umap$sampleID)), collapse=",") 
    AUCell_gene_set_name<-paste(AUCell_gene_set_name," [show only: ",names.sampleIDs,"]",sep="")
    ##########
    clusters.info<-get_clusters(object)
    
    if(clustering_method=="walktrap"){
      clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="louvain"){
      clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_walktrap"){
      clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_louvain"){
      clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else{
      stop("Need community detection or clustering-method names!")
    }
    plot.umap<-merge(plot.umap,clusters.info)
    sample.index<-which(plot.umap$cluster %in% selected.limited.clusters)
    edited_score<-rep(0,nrow(plot.umap))
    edited_score[sample.index]<-plot.umap$AUCell_Score[sample.index]
    plot.umap$AUCell_Score<-edited_score
    total.n.cells<-nrow(plot.umap)
    #################
    z<-plot.umap[sample.index,]
    z<-subset(z,AUCell_Score > AUCell_cutoff)
    AUC.n.cells<-nrow(z)
    #################
    txt<-paste("total cells=",total.n.cells," ; cells with AUCell_score=",AUC.n.cells,sep="")
    AUCell_gene_set_name<-paste(AUCell_gene_set_name,txt,sep="\n")
  }
  ####################################
  ggplot(plot.umap,aes(UMAP1,UMAP2, color=AUCell_Score)) + geom_point(size=point.size)+ggtitle(AUCell_gene_set_name)+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_colour_gradientn(colours = c("gray85",point.color1,point.color2),values=c(0,0.1,1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
}

#' Plot TSNE with the expression of multiple gene sets on selected sample IDs
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file.
#' @param  selected.sampleID Selected sample IDs
#' @param  show_gene_sets The name of gene sets to be shown.
#' @param  custom_color The vector of cutom color to specify gene sets.
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  IsDownsample is logical. If TRUE, the cells will be downsampled.
#' @param  downsample.size The number of cells for the downsample. Default 0
#' @param  title.color The color of the title. Default red
#' @param  title.font.size The title font size. Default 16
#' @param  point.size The point size. Default 2
#' @param  showLegend is logical. If TRUE, the figure legend will be displayed.
#' @param  background.color The background color. Default white
#' @export

plot_tsne_label_by_multiple_gene_sets_with_limited_sampleID<-function(object,gmt.file=c(),
                                selected.sampleID=c(),title="",show_gene_sets=c(),custom_color=c(),
																isNormalizedByHouseKeeping=T,IsDownsample=F,downsample.size=0,
																title.color="red",title.font.size=16,point.size=2,showLegend=T,
																background.color="white"){
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gmt.file)==0){
		stop("Required the path to GMT file!")
	}
	if(length(show_gene_sets)==0){
		stop("Required gene set names!")
	}
	#if(length(show_gene_sets)==1){
	#	stop("Required at least two gene sets")
	#}
	if(length(show_gene_sets) > 10){
		stop("This is over limitation at 10 gene sets!")
	}
	if(length(custom_color) > 0){
		if(length(custom_color) !=length(show_gene_sets)){
			stop("The number of colors is not equal to the the number of input gene sets.")
		}
	}
	#####
	signature.sets<-get_gmtGeneSets(gmt.file)
	#####
	check.1<-match(show_gene_sets,names(signature.sets))
	check.2<-sum(length(which(is.na(check.1))))
	if(check.2 > 0){
		stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
	}
	res.tsne<-get_tsne.result(object)
  	n.cells<-0
  	if(IsDownsample==FALSE){
    	sample.index<-which(res.tsne$sampleID %in% selected.sampleID)
    	res.tsne<-res.tsne[sample.index,]
    	n.cells<-nrow(res.tsne)
    	
  	}else if(IsDownsample==TRUE & downsample.size>0){
    	sample.index<-which(res.tsne$sampleID %in% selected.sampleID)
    	sample.index.s<-sample(sample.index,size=downsample.size)
    	res.tsne<-res.tsne[sample.index.s,]
    	n.cells<-nrow(res.tsne)
  	}
  	title<-paste(title," [",n.cells," cells]",sep="")
	
	#####
	umi<-get_umi_count(object)
	umi.used<-umi[,as.character(res.tsne$Cell)]
	#####
	my.dat<-res.tsne[,c("Cell","TSNE1","TSNE2")]
	#####
	genes.exp.sum<-Matrix::rowSums(umi.used)
	genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
	hk.genes<-names(genes.exp.sum[1:100])
	hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
	#####
	Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
	for(i in 1:length(show_gene_sets)){
		set.name<-show_gene_sets[i]
		genes.x<-signature.sets[[set.name]]
		exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
		
		if(isNormalizedByHouseKeeping==T){
			MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
		}else{
			MM<-Matrix::colSums(exprs)/res.tsne$UMI_count
		}
		MM.f<-data.frame(score=MM)
		colnames(MM.f)<-set.name
		Genes.score<-cbind(Genes.score,MM.f)
	}
	Genes.score<-Genes.score[-c(1)]
	####Assign colors#####################
	if(length(custom_color)==0){
		multi.colors<-rainbow(ncol(Genes.score))
	}else{
		multi.colors<-custom_color
	}
	######################################
	my.colors.list<-list()
	for(j in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[j])
		x.color<-multi.colors[j]
		p1<-ggplot(my.dat,aes(TSNE1,TSNE2)) + geom_point(aes(colour = Genes.score[,j])) + scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
		My_color<-ggplot_build(p1)$data[[1]]$colour
		my.colors.list[[set.name]]<-My_color
	}
	
	my.alpha.list<-list()
	for(j in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[j])
		my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
		my.alpha.list[[set.name]]<-my.alpha
	}
	###################################
	p<-ggplot(my.dat,aes(TSNE1,TSNE2))
	for(k in 1:ncol(Genes.score)){
		set.name<-colnames(Genes.score[k])
		p<-p+geom_point(color=my.colors.list[[set.name]], size = point.size,alpha=my.alpha.list[[set.name]],stroke = 0)
	}
	
	p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
			theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
			theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+ggtitle(title)+
    		theme(plot.title = element_text(color=title.color, size=title.font.size, face="bold.italic",hjust = 0.5))+
    		theme(plot.background = element_rect(fill = background.color))
	
	if(showLegend==T){
		df<-data.frame(Name=show_gene_sets)
		df<-as.matrix(df)
	
			tp<-ggtexttable(df,rows = NULL,
			theme = ttheme(
					colnames.style = colnames_style(fill = "white"),
					tbody.style = tbody_style(fill = multi.colors)
			)
		)
		return(ggarrange(p,tp,
					ncol = 2, nrow = 1,
					widths = c(1,0.2)))
	}else{
		return(p)
	}
}

#' Plot UMAP with the expression of multiple gene sets on selected sample IDs
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file.
#' @param  selected.sampleID Selected sample IDs
#' @param  title The title of the plot.
#' @param  show_gene_sets The name of gene sets to be shown.
#' @param  custom_color The vector of cutom color to specify gene sets.
#' @param  IsDownsample is logical. If TRUE, the cells will be downsampled.
#' @param  downsample.size The number of cells for the downsample. Default 0
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  point.size The point size. Default 2
#' @param  title.color The color of the title. Default red
#' @param  title.font.size The title font size. Default 16
#' @param  showLegend is logical. If TRUE, the figure legend will be displayed.
#' @param  background.color The background color. Default white
#' @export


plot_umap_label_by_multiple_gene_sets_with_limited_sampleID<-function(object,gmt.file=c(),selected.sampleID=c(),title="",show_gene_sets=c(),
                                                                      custom_color=c(),IsDownsample=F,downsample.size=0,
                                                                      isNormalizedByHouseKeeping=T,point.size=2,title.color="red",
                                                                      title.font.size=16,showLegend=T,background.color="white"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(selected.sampleID)==0){
    stop("Please input selected sample IDs!")
  }
  
  if(length(gmt.file)==0){
    stop("Required the path to GMT file!")
  }
  if(length(show_gene_sets)==0){
    stop("Required gene set names!")
  }
  
  #if(length(show_gene_sets)==1){
  #	stop("Required at least two gene sets")
  #}
  if(length(show_gene_sets) > 10){
    stop("This is over limitation at 10 gene sets!")
  }
  if(length(custom_color) > 0){
    if(length(custom_color) !=length(show_gene_sets)){
      stop("The number of colors is not equal to the the number of input gene sets.")
    }
  }
  #####
  signature.sets<-get_gmtGeneSets(gmt.file)
  #####
  check.1<-match(show_gene_sets,names(signature.sets))
  check.2<-sum(length(which(is.na(check.1))))
  if(check.2 > 0){
    stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
  }
  
  res.umap<-get_umap.result(object)
  n.cells<-0
  if(IsDownsample==FALSE){
    sample.index<-which(res.umap$sampleID %in% selected.sampleID)
    res.umap<-res.umap[sample.index,]
    n.cells<-nrow(res.umap)
  }else if(IsDownsample==TRUE & downsample.size>0){
    sample.index<-which(res.umap$sampleID %in% selected.sampleID)
    sample.index.s<-sample(sample.index,size=downsample.size)
    res.umap<-res.umap[sample.index.s,]
    n.cells<-nrow(res.umap)
  }
  title<-paste(title," [",n.cells," cells]",sep="")
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(res.umap$Cell)]
  #####
  my.dat<-res.umap[,c("Cell","UMAP1","UMAP2")]
  #####
  genes.exp.sum<-Matrix::rowSums(umi.used)
  genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
  hk.genes<-names(genes.exp.sum[1:100])
  hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
  #####
  Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
  for(i in 1:length(show_gene_sets)){
    set.name<-show_gene_sets[i]
    genes.x<-signature.sets[[set.name]]
    exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
    
    if(isNormalizedByHouseKeeping==T){
      MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
    }else{
      MM<-Matrix::colSums(exprs)/res.umap$UMI_count
    }
    MM.f<-data.frame(score=MM)
    colnames(MM.f)<-set.name
    Genes.score<-cbind(Genes.score,MM.f)
  }
  Genes.score<-Genes.score[-c(1)]
  ####Assign colors#####################
  if(length(custom_color)==0){
    multi.colors<-rainbow(ncol(Genes.score))
  }else{
    multi.colors<-custom_color
  }
  ######################################
  my.colors.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    x.color<-multi.colors[j]
    p1<-ggplot(my.dat,aes(UMAP1,UMAP2)) + geom_point(aes(colour = Genes.score[,j])) + scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
    My_color<-ggplot_build(p1)$data[[1]]$colour
    my.colors.list[[set.name]]<-My_color
  }
  
  my.alpha.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
    my.alpha.list[[set.name]]<-my.alpha
  }
  ###################################
  p<-ggplot(my.dat,aes(UMAP1,UMAP2))
  for(k in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[k])
    p<-p+geom_point(color=my.colors.list[[set.name]], size = point.size,alpha=my.alpha.list[[set.name]],stroke = 0)
  }
  
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+ggtitle(title)+
    theme(plot.title = element_text(color=title.color, size=title.font.size, face="bold.italic",hjust = 0.5))+
    theme(plot.background = element_rect(fill = background.color))
  
  if(showLegend==T){
  		df<-data.frame(Name=show_gene_sets)
  		df<-as.matrix(df)
  
  		tp<-ggtexttable(df,rows = NULL,
                  theme = ttheme(
                    colnames.style = colnames_style(fill = "white"),
                    tbody.style = tbody_style(fill = multi.colors)
                  )
  		)
 
  		return(ggarrange(p,tp,
                   ncol = 2, nrow = 1,
                   widths = c(1,0.2)))
  }else{
  		return(p)
  }
}

#' Plot force-directed graph with the expression of multiple gene sets on selected sample IDs
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file.
#' @param  selected.sampleID Selected sample IDs
#' @param  title The title of the plot.
#' @param  show_gene_sets The name of gene sets to be shown.
#' @param  custom_color The vector of cutom color to specify gene sets.
#' @param  IsDownsample is logical. If TRUE, the cells will be downsampled.
#' @param  downsample.size The number of cells for the downsample. Default 0
#' @param  isNormalizedByHouseKeeping is logical. If TRUE, the gene expression will be normalized by the expression of house keeping genes.
#' @param  vertex.size The point size. Default 1.5
#' @param  edge.size The edge size. Default 0.2
#' @param  edge.color The edge color. Default gray
#' @param  title.color The color of the title. Default red
#' @param  title.font.size The title font size. Default 16
#' @param  showLegend is logical. If TRUE, the figure legend will be displayed.
#' @param  background.color The background color. Default white
#' @export

plot_forceDirectedGraph_label_by_multiple_gene_sets_with_limited_sampleID<-function(object,gmt.file=c(),selected.sampleID=c(),title="",
                                                                                    show_gene_sets=c(),custom_color=c(),IsDownsample=F,
                                                                                    downsample.size=0,isNormalizedByHouseKeeping=T,
                                                                                    vertex.size=1.5,edge.size=0.2,edge.color="gray",
                                                                                    title.color="red",title.font.size=16,
                                                                                    showLegend=T,background.color="white"){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gmt.file)==0){
    stop("Required the path to GMT file!")
  }
  if(length(show_gene_sets)==0){
    stop("Required gene set names!")
  }
  if(length(show_gene_sets)==1){
    stop("Required at least two gene sets")
  }
  if(length(show_gene_sets) > 10){
    stop("This is over limitation at 10 gene sets!")
  }
  if(length(custom_color) > 0){
    if(length(custom_color) !=length(show_gene_sets)){
      stop("The number of colors is not equal to the the number of input gene sets.")
    }
  }
  if(length(selected.sampleID)==0){
    stop("Please input selected sample IDs!")
  }
  #####
  signature.sets<-get_gmtGeneSets(gmt.file)
  #####
  check.1<-match(show_gene_sets,names(signature.sets))
  check.2<-sum(length(which(is.na(check.1))))
  if(check.2 > 0){
    stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
  }
  fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
  colnames(fa2.layout) <- c("X1","X2")
  #########
  cell.info<-get_cells_annotation(object)
  cell.info<-subset(cell.info,IsPassed==T)
  rownames(cell.info)<-cell.info$Cell
  #########
  cell.info<-cell.info[rownames(fa2.layout),c("sampleID","UMI_count","detectedGenesPerCell","percent_mito")]
  fa2.layout<-cbind(fa2.layout,cell.info)
  
  n.cells<-0
  if(IsDownsample==FALSE){
    sample.index<-which(fa2.layout$sampleID %in% selected.sampleID)
    fa2.layout<-fa2.layout[sample.index,]
    n.cells<-nrow(fa2.layout)
  }else if(IsDownsample==TRUE & downsample.size>0){
    sample.index<-which(fa2.layout$sampleID %in% selected.sampleID)
    sample.index.s<-sample(sample.index,size=downsample.size)
    fa2.layout<-fa2.layout[sample.index.s,]
    n.cells<-nrow(fa2.layout)
  }
  title<-paste(title," [",n.cells," cells]",sep="")
  
  #####
  umi<-get_umi_count(object)
  umi.used<-umi[,as.character(rownames(fa2.layout))]
  #####
  UMI_count<-Matrix::colSums(umi.used)
  my.dat<-fa2.layout[,(1:2)]
  #####
  genes.exp.sum<-Matrix::rowSums(umi.used)
  genes.exp.sum<-genes.exp.sum[order(genes.exp.sum,decreasing=T)]
  hk.genes<-names(genes.exp.sum[1:100])
  hk_exprs<-as.matrix(umi.used[rownames(umi.used) %in% hk.genes,])
  #####
  Genes.score<-data.frame(score=rep(0,nrow(my.dat)))
  for(i in 1:length(show_gene_sets)){
    set.name<-show_gene_sets[i]
    genes.x<-signature.sets[[set.name]]
    exprs<-as.matrix(umi.used[rownames(umi.used) %in% genes.x,])
    
    if(isNormalizedByHouseKeeping==T){
      MM <- Matrix::colSums(exprs)/Matrix::colSums(hk_exprs)
    }else{
      MM<-Matrix::colSums(exprs)/UMI_count
    }
    MM.f<-data.frame(score=MM)
    colnames(MM.f)<-set.name
    Genes.score<-cbind(Genes.score,MM.f)
  }
  Genes.score<-Genes.score[-c(1)]
  ####Assign colors#####################
  if(length(custom_color)==0){
    multi.colors<-rainbow(ncol(Genes.score))
  }else{
    multi.colors<-custom_color
  }
  ######################################
  my.colors.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    x.color<-multi.colors[j]
    p1<-ggplot(my.dat,aes(X1,X2)) + geom_point(aes(colour = Genes.score[,j])) + scale_colour_gradientn(colours = c("gray85",x.color,x.color),values=c(0,0.1,1))
    My_color<-ggplot_build(p1)$data[[1]]$colour
    my.colors.list[[set.name]]<-My_color
  }
  my.alpha.list<-list()
  for(j in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[j])
    my.alpha<-as.numeric(Genes.score[,j]/rowSums(Genes.score))
    my.alpha.list[[set.name]]<-my.alpha
  }
  ###################################
  rownames(fa2.layout) <- 1:nrow(fa2.layout)
  ####################################
  #edgelist <- get.edgelist(get_igraph.graph(object))
  #edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
  #colnames(edges) <- c("X1","Y1","X2","Y2")
  ###################################
  p<-ggplot(my.dat,aes(X1,X2))
    #geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)
  for(k in 1:ncol(Genes.score)){
    set.name<-colnames(Genes.score[k])
    p<-p+geom_point(color=my.colors.list[[set.name]], size = vertex.size,alpha=my.alpha.list[[set.name]],stroke = 0)
  }
  
  p<-p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+ggtitle(title)+
    theme(plot.title = element_text(color=title.color, size=title.font.size, face="bold.italic",hjust = 0.5))+
    theme(plot.background = element_rect(fill = background.color))
  
  if(showLegend==T){
    df<-data.frame(Name=show_gene_sets)
    df<-as.matrix(df)
  
    tp<-ggtexttable(df,rows = NULL,
                  theme = ttheme(
                    colnames.style = colnames_style(fill = "white"),
                    tbody.style = tbody_style(fill = multi.colors)
                  )
    )
    return(ggarrange(p,tp,ncol = 2, nrow = 1,widths = c(1,0.2)))
  }else{
    return(p)
  }
}

#' Plot force-direct graph with AUCell score
#' @param  object The SingCellaR object.
#' @param  AUCell_gene_set_name The vector of gene set names.
#' @param  AUCell_score The dataframe of AUCell score per gene set.
#' @param  AUCell_cutoff The cutoff score. Default 0
#' @param  IsShowOnlySampleIDs is logical. If TRUE, only selected sample IDs will be shown.
#' @param  selected.sampleID The selected sample IDs.
#' @param  IsLimitedAUCscoreByClusters is logical. If TRUE, AUCell score will be limited by selected clusters.
#' @param  selected.limited.clusters The selected clusters of interest.
#' @param  clustering_method The clustering method.
#' @param  vertex.size The node size. Default 1.5
#' @param  edge.size The edge size. Default 0.2
#' @param  vertex.color1 The color for a low AUCell score. Default black
#' @param  vertex.color2 The color for a high AUCell score. Default red
#' @param  edge.color The edge color.
#' @param  background.color The background color.
#' @param  showEdge is logical. If TRUE, the edge will be displayed.
#' @export

plot_forceDirectedGraph_label_by_AUCell_score<-function(object,AUCell_gene_set_name=c(),AUCell_score,AUCell_cutoff=0,
                                                        IsShowOnlySampleIDs=FALSE,selected.sampleID=c(),
                                                        IsLimitedAUCscoreByClusters=FALSE,selected.limited.clusters=c(),
                                                        clustering_method="louvain",vertex.size=1.5,edge.size=0.2,
                                                        vertex.color1="black",vertex.color2="red",
                                                        edge.color="gray",background.color="white",showEdge=T){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(AUCell_gene_set_name)==0){
    stop("Required a name of gene set!")
  }
  #####
  fa2.layout<-as.data.frame(get_fa2_graph.layout(object))
  colnames(fa2.layout)<-c("X1","X2")
  fa2.layout$Cell<-rownames(fa2.layout)
  #####
  my.AUCell<-AUCell_score[as.character(fa2.layout$Cell),]
  my.AUCell.score<-my.AUCell[,which(colnames(my.AUCell)==AUCell_gene_set_name)]
  
  #####
  Genes.score<-data.frame(AUCell_Score=my.AUCell.score)
  f.umap<-cbind(fa2.layout,Genes.score)
  
  cell.anno<-get_cells_annotation(object)
  rownames(cell.anno)<-cell.anno$Cell
  cell.anno<-cell.anno[,-c(1)]
  cell.anno<-cell.anno[as.character(fa2.layout$Cell),]
  f.umap<-cbind(f.umap,cell.anno)
  ######
  plot.fa2<-f.umap
  plot.fa2$AUCell_Score[plot.fa2$AUCell_Score < AUCell_cutoff]<-0
  AUCell_gene_set_name<-paste(AUCell_gene_set_name," (Cells with AUC >",AUCell_cutoff, ")",sep="")
  ######
  rownames(fa2.layout) <- 1:nrow(fa2.layout)
  fa2.layout<-fa2.layout[,c(1:2)]
  ####################################
  edgelist <- get.edgelist(get_igraph.graph(object))
  edges <- data.frame(fa2.layout[edgelist[,1],], fa2.layout[edgelist[,2],])
  colnames(edges) <- c("X1","Y1","X2","Y2")
  
  
  if(IsShowOnlySampleIDs==T & IsLimitedAUCscoreByClusters==F){
    sample.index<-which(plot.fa2$sampleID %in% selected.sampleID)
    plot.fa2<-plot.fa2[sample.index,]
    names.sampleIDs<-paste(names(table(plot.fa2$sampleID)), collapse=",") 
    AUCell_gene_set_name<-paste(AUCell_gene_set_name," [show only: ",names.sampleIDs,"]",sep="")
    
  }else if(IsShowOnlySampleIDs==F & IsLimitedAUCscoreByClusters==T){
    
    clusters.info<-get_clusters(object)
    
    if(clustering_method=="walktrap"){
      clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="louvain"){
      clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_walktrap"){
      clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_louvain"){
      clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else{
      stop("Need community detection or clustering-method names!")
    }
    
    plot.fa2<-merge(plot.fa2,clusters.info)
    sample.index<-which(plot.fa2$cluster %in% selected.limited.clusters)
    edited_score<-rep(0,nrow(plot.fa2))
    edited_score[sample.index]<-plot.fa2$AUCell_Score[sample.index]
    plot.fa2$AUCell_Score<-edited_score
    total.n.cells<-nrow(plot.fa2)
    #################
    z<-plot.fa2[sample.index,]
    z<-subset(z,AUCell_Score > AUCell_cutoff)
    AUC.n.cells<-nrow(z)
    #################
    txt<-paste("total cells=",total.n.cells," ; cells with AUCell_score=",AUC.n.cells,sep="")
    AUCell_gene_set_name<-paste(AUCell_gene_set_name,txt,sep="\n")
    
  }else if(IsShowOnlySampleIDs==T & IsLimitedAUCscoreByClusters==T){
    sample.index<-which(plot.fa2$sampleID %in% selected.sampleID)
    plot.fa2<-plot.fa2[sample.index,]
    ##########
    names.sampleIDs<-paste(names(table(plot.fa2$sampleID)), collapse=",") 
    AUCell_gene_set_name<-paste(AUCell_gene_set_name," [show only: ",names.sampleIDs,"]",sep="")
    ##########
    clusters.info<-get_clusters(object)
    
    if(clustering_method=="walktrap"){
      clusters.info<-clusters.info[,c("Cell","walktrap_cluster","walktrap_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="louvain"){
      clusters.info<-clusters.info[,c("Cell","louvain_cluster","louvain_cluster_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_walktrap"){
      clusters.info<-clusters.info[,c("Cell","merged_walktrap","merged_walktrap_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_louvain"){
      clusters.info<-clusters.info[,c("Cell","merged_louvain","merged_louvain_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else if(clustering_method=="merged_kmeans"){
      clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans","merged_kmeans_color")]
      colnames(clusters.info)<-c("Cell","cluster","cluster_color")
    }else{
      stop("Need community detection or clustering-method names!")
    }
    plot.fa2<-merge(plot.fa2,clusters.info)
    sample.index<-which(plot.fa2$cluster %in% selected.limited.clusters)
    edited_score<-rep(0,nrow(plot.fa2))
    edited_score[sample.index]<-plot.fa2$AUCell_Score[sample.index]
    plot.fa2$AUCell_Score<-edited_score
    total.n.cells<-nrow(plot.fa2)
    #################
    z<-plot.fa2[sample.index,]
    z<-subset(z,AUCell_Score > AUCell_cutoff)
    AUC.n.cells<-nrow(z)
    #################
    txt<-paste("total cells=",total.n.cells," ; cells with AUCell_score=",AUC.n.cells,sep="")
    AUCell_gene_set_name<-paste(AUCell_gene_set_name,txt,sep="\n")
  }
  ################# 
  if(showEdge==T){
      p<-qplot(X1,X2, data=plot.fa2,colour=AUCell_Score)+
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = edge.size, colour=edge.color)
  }else{
      p<-qplot(X1,X2, data=plot.fa2,colour=AUCell_Score)
  }
    p<-p+geom_point(size=vertex.size)+ggtitle(AUCell_gene_set_name)+
    scale_colour_gradientn(colours = c("gray85",vertex.color1,vertex.color2),values=c(0,0.1,1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank())+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.title.y=element_blank())+
    theme(plot.background = element_rect(fill = background.color))
   
    return(p)
}

#' Plot bubble-plot for gene expression per cluster
#' @param  object The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  gene_list The vector of gene names.
#' @param  point.color1 The color for low expression level. Default orange
#' @param  point.color2 The color for high expression level. Default firebrick
#' @param  buble.scale The scale size of the bubble. Default 6
#' @param  gene.font.size The font size of gene name. Default 12
#' @param  axist.x.font.size The label on the x-axis font size. Default 12
#' @param  show.percent is logical. If TRUE, the percentage will be displayed.
#' @param  log is logical. If TRUE, the log expression will be applied.
#' @param  font.size.in.bubble The font size inside the bubble. Default 2
#' @export

plot_bubble_for_genes_per_cluster<-function(object,cluster.type=c("louvain","walktrap","kmeans"),
                                            gene_list=c(),point.color1="orange",point.color2="firebrick",
                                            buble.scale=6,gene.font.size=12,axist.x.font.size=12,
                                            show.percent=FALSE,log=T,font.size.in.bubble=2){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  ####################################
  if(cluster.type=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster")]
  }else if(cluster.type=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)
  }else if(cluster.type=="merged_walktrap"){
    clusters.info<-clusters.info[,c("Cell","merged_walktrap")]
  }else if(cluster.type=="merged_louvain"){
    clusters.info<-clusters.info[,c("Cell","merged_louvain")]
  }else if(cluster.type=="merged_kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans")]
  }else{
    stop("Need a clustering-method name!")
  }
  clusters.ids<-1:length(unique(clusters.info[,2]))
  ####################################
  umi.dat<-get_normalized_umi(object)
  ####################################
  my.new.dat<-data.frame()
  for(i in 1:length(gene_list)){
    my.gene<-as.character(gene_list)[i]
    my.sub.dat<-umi.dat[rownames(umi.dat)==my.gene,]
    for(k in 1:length(clusters.ids)){
      cluster.id<-paste("cl",k,sep="")
      x.cell.names<-subset(clusters.info,clusters.info[,2]==cluster.id)
      x.expr<-my.sub.dat[names(my.sub.dat) %in% as.character(x.cell.names$Cell)]
      ##################
      exp.val<-as.numeric(x.expr)
      n.exp<-length(exp.val[exp.val>0])
      n.total<-length(exp.val)
      exp.fraction<-round((n.exp/n.total)*100)
      exp.mean<-mean(exp.val)
      if(log==T){
        exp.mean<-log1p(exp.mean)
      }
      ##################
      my.f<-data.frame(gene=my.gene,cluster=cluster.id,percent_express=exp.fraction,Average_expression=exp.mean)
      my.new.dat<-rbind(my.new.dat,my.f)
    }
  }
    p<-ggplot(my.new.dat, aes(x=cluster, y=gene)) + 
    geom_point(aes(size=percent_express,colour = Average_expression)) + 
    scale_size(range = c(0, buble.scale))+ 
    scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.y = element_text(face="italic",size=gene.font.size),axis.title.y = element_blank())+
    theme(axis.text.x = element_text(face="bold",size=axist.x.font.size),axis.title.x = element_blank())+
    guides(size = guide_legend(title = 'Expressing cells (%)',order=2))
    
    if(show.percent==TRUE){
      p<-p+geom_text(aes(label=percent_express),size=font.size.in.bubble)
    }
    return(p)
}

#' Plot bubble-plot for gene expression per custom group of cells
#' @param  object The SingCellaR object.
#' @param  custom_group_of_cells The list of custome group of cells.
#' @param  gene_list The vector of gene names.
#' @param  point.color1 The color for low expression level. Default orange
#' @param  point.color2 The color for high expression level. Default firebrick
#' @param  buble.scale The scale size of the bubble. Default 6
#' @param  gene.font.size The font size of gene name. Default 12
#' @param  axist.x.font.size The label on the x-axis font size. Default 12
#' @param  show.percent is logical. If TRUE, the percentage will be displayed.
#' @param  log is logical. If TRUE, the log expression will be applied.
#' @param  font.size.in.bubble The font size inside the bubble. Default 2
#' @export

plot_bubble_for_genes_per_custom_group_of_cells<-function(object,custom_group_of_cells=list(),
                                            gene_list=c(),point.color1="orange",point.color2="firebrick",
                                            buble.scale=6,gene.font.size=12,axist.x.font.size=12,font.size.in.bubble=2,
                                            show.percent=FALSE,log=T){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  
  ####################################
  umi.dat<-get_normalized_umi(object)
  ####################################
  my.new.dat<-data.frame()
  for(i in 1:length(gene_list)){
    my.gene<-as.character(gene_list)[i]
    my.sub.dat<-umi.dat[rownames(umi.dat)==my.gene,]
    for(k in 1:length(custom_group_of_cells)){
      x.cell.names<-custom_group_of_cells[[k]]
      x.expr<-my.sub.dat[names(my.sub.dat) %in% as.character(x.cell.names)]
      ##################
      exp.val<-as.numeric(x.expr)
      n.exp<-length(exp.val[exp.val>0])
      n.total<-length(exp.val)
      exp.fraction<-round((n.exp/n.total)*100)
      exp.mean<-mean(exp.val)
      if(log==TRUE){
        exp.mean<-log1p(exp.mean)
      }
      ##################
      my.f<-data.frame(gene=my.gene,cluster=names(custom_group_of_cells)[k],percent_express=exp.fraction,Average_expression=exp.mean)
      my.new.dat<-rbind(my.new.dat,my.f)
    }
  }
  p<-ggplot(my.new.dat, aes(x=cluster, y=gene)) + 
    geom_point(aes(size=percent_express,colour = Average_expression)) + 
    scale_size(range = c(0, buble.scale))+ 
    scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text.y = element_text(face="italic",size=gene.font.size),axis.title.y = element_blank())+
    theme(axis.text.x = element_text(face="bold",size=axist.x.font.size),axis.title.x = element_blank())+
    guides(size = guide_legend(title = 'Expressing cells (%)',order=2))
  
  if(show.percent==TRUE){
    p<-p+geom_text(aes(label=percent_express),size=font.size.in.bubble)
  }
  return(p)
}

#' Plot violin-plot for gene expression per custom group of cells
#' @param  object The SingCellaR object.
#' @param  custom_group_of_cells The list of custome group of cells.
#' @param  gene_list The vector of gene names.
#' @param  take_log2 is logical. If TRUE, log2 expression will be applied.
#' @param  xlab.text.size The font size of the label on the x-axis. Default 5
#' @param  point.size The size of data point. Default 0.2
#' @param  point.alpha The alpha parameter of the data point. Default 0.1
#' @param  grid.ncol the column number of the grid of the containing plots to be displayed. Default 3
#' @param  grid.nrow the row number of the grid of the containing plots to be displayed. Default 3
#' @export

plot_violin_for_genes_per_custom_group_of_cells<-function(object,custom_group_of_cells=list(),gene_list=c(),
                                                          take_log2=T,xlab.text.size=5,point.size=0.2,point.alpha=0.1,
                                                          grid.ncol = 3,grid.nrow=3){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gene_list)==0){
    stop("Required list of genes!")
  }
  
  ####################################
  umi.dat<-get_normalized_umi(object)
  my.p<-list()
  for(i in 1:length(gene_list)){
    my.gene<-as.character(gene_list)[i]
    my.sub.dat<-umi.dat[rownames(umi.dat)==my.gene,]
    my.new.dat<-data.frame()
    for(k in 1:length(custom_group_of_cells)){
      x.cell.names<-custom_group_of_cells[[k]]
      x.expr<-my.sub.dat[names(my.sub.dat) %in% as.character(x.cell.names)]
      ##################
      exp.val<-as.numeric(x.expr)
      n.exp<-length(exp.val[exp.val>0])
      n.total<-length(exp.val)
      my.label<-paste(names(custom_group_of_cells)[k],"\n",n.exp,"/",n.total,sep="")
      ##################
      y.legend<-""
      ##################
      if(take_log2==T){
        my.f<-data.frame(cluster=my.label,normUMI=log2(as.numeric(x.expr)+1))
        y.legend<-"log2(normalized UMI)"
      }else{
        my.f<-data.frame(cluster=my.label,normUMI=as.numeric(x.expr))
        y.legend<-"normalized UMI"
      }
      my.new.dat<-rbind(my.new.dat,my.f)
    }
    my.p[[i]]<-ggplot(
      data = my.new.dat,aes(x=cluster, y=normUMI,fill=cluster))+labs(y=y.legend)+
      geom_violin(trim=TRUE,scale="width")+
      geom_jitter(height = 0,size=point.size,alpha = point.alpha)+
      ggtitle(my.gene)+
      guides(fill=FALSE)+
      theme(axis.text=element_text(size=xlab.text.size))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } 
  grid.arrange(grobs = my.p,nrow=grid.nrow,ncol=grid.ncol)
}