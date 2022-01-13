###############################################################################
#getNumberOfDetectedCellsPerGene <- function(expression_matrix, min_expr){
#
#	FM <- expression_matrix
#	FM_genes <- do.call(rbind, apply(FM, 1,
#					function(x) {
#						return(data.frame(
#										num_detected_cells=sum(unlist(as.list(x)) >= min_expr)
#								)
#						)
#					})
#	)
#	return(FM_genes)
#}

#getNumberOfDetectedGenesPerCell <- function(expression_matrix, min_expr){
#
#	FM <- expression_matrix
#	FM_cells <- do.call(rbind, apply(FM, 2,
#					function(x) {
#						return(data.frame(
#										num_genes_detected=sum(unlist(as.list(x)) >= min_expr)
#								)
#						)
#					})
#	)
#	return(FM_cells)
#}

#' List of available gene sets
#' @param  gmt.file The GMT file name.
#' @export 

listAvailableGeneSets <- function(gmt.file) {
  Lines <- strsplit(readLines(gmt.file), "\t")
  genesets <- lapply(Lines, utils::tail, -2)
  names(genesets) <- sapply(Lines, head, 1)
  print(names(genesets))
}

#' Get gene sets from GMT file
#' @param  gmt.file The GMT file name.
#' @export 
#' @return List of gene sets.

get_gmtGeneSets <- function(gmt.file) {
	Lines <- strsplit(readLines(gmt.file), "\t")
	genesets <- lapply(Lines, utils::tail, -2)
	names(genesets) <- sapply(Lines, head, 1)
	return(genesets)
}


CustomPalette <- function(low = "white",high = "red",mid = NULL,k = 50) {
  low <- col2rgb(col = low) / 255
  high <- col2rgb(col = high) / 255
  if (is.null(x = mid)) {
    r <- seq(from = low[1], to = high[1], len = k)
    g <- seq(from = low[2], to = high[2], len = k)
    b <- seq(from = low[3], to = high[3], len = k)
  } else {
    k2 <- round(x = k / 2)
    mid <- col2rgb(col = mid) / 255
    r <- c(
      seq(from = low[1], to = mid[1], len = k2),
      seq(from = mid[1], to = high[1], len = k2)
    )
    g <- c(
      seq(from = low[2], to = mid[2], len = k2),
      seq(from = mid[2], to = high[2],len = k2)
    )
    b <- c(
      seq(from = low[3], to = mid[3], len = k2),
      seq(from = mid[3], to = high[3], len = k2)
    )
  }
  return(rgb(red = r, green = g, blue = b))
}

identifyKNN <- function(Xmat,metric=c("euclidean","cosine"),n.neighbors=30,n_trees=50,n_threads =1){
  
  print(paste("Building Annoy index with metric = ", metric, ", n_trees = ", n_trees),sep=" ")
  ####Build annoy######################
  metric<-match.arg(metric)
  if(metric=="euclidean"){
    my.annoy <- new(RcppAnnoy::AnnoyEuclidean, ncol(Xmat))
    search_knn <- new(RcppAnnoy::AnnoyEuclidean, ncol(Xmat))  
  }else if(metric=="cosine"){
    my.annoy <- new(RcppAnnoy::AnnoyAngular, ncol(Xmat))
    search_knn <- new(RcppAnnoy::AnnoyAngular, ncol(Xmat))  
  }else{
    stop("Please input the type of metric!")
  }
  # Add items
  progress <- Progress_bar$new(max = nrow(Xmat), display = TRUE)
  for (i in 1:nrow(Xmat)) {
    my.annoy$addItem(i - 1, Xmat[i, ])
    progress$increment()
  }
  # Build index
  my.annoy$build(n_trees)
  index_file <- tempfile(tmpdir = tempdir())
  my.annoy$save(index_file)
  #####################################
  ##Search knn#########################
  search_knn$load(index_file)
  #####################################
  if(n_threads > 1){
    
    print(paste("Searching Annoy index using ",n_threads," threads search_k  =",n.neighbors*100,sep=""))
    res <- annoy_search_parallel_cpp(index_file,Xmat,n.neighbors, n.neighbors*100,metric,grain_size = 1,n_threads=n_threads)
    knn<-res$item+1
  }else{
    print(paste("Searching Annoy index, search_k = ", n.neighbors*100),sep="")
    search_progress <- Progress_bar$new(max = nrow(Xmat), display = TRUE)
    idx <- matrix(nrow = nrow(Xmat), ncol = n.neighbors)
    for (i in 1:nrow(Xmat)) {
      res <- search_knn$getNNsByVectorList(Xmat[i, ], n.neighbors, n.neighbors*100, FALSE)
      idx[i, ] <- res$item
      search_progress$increment()
    }
    knn<-idx+1
  }
  #########remove index_file###
  unlink(index_file)
  #############################
  return(knn)
}

Progress_bar <- setRefClass("Progress_bar",
                        fields = list(
                          value = "numeric",
                          max = "numeric",
                          curr_stars = "numeric",
                          max_stars = "numeric",
                          display = "logical"
                        ),
                        methods = list(
                          initialize = function(max, display = TRUE) {
                            max_stars <<- 51 # length of the progress bar
                            value <<- 0
                            curr_stars <<- 0
                            max <<- max
                            display <<- display
                            if (display) {
                              message("0%   10   20   30   40   50   60   70   80   90   100%")
                              message("[----|----|----|----|----|----|----|----|----|----|")
                            }
                          },
                          increment = function() {
                            if (display && curr_stars < max_stars) {
                              value <<- value + 1
                              num_stars <- round(max_stars * value / max)
                              if (num_stars > curr_stars) {
                                # Number of new stars to print
                                num_new_stars <- num_stars - curr_stars
                                
                                # If we are going to reach the end of the progress bar
                                # save space for the terminal "|"
                                if (curr_stars + num_new_stars >= max_stars) {
                                  num_new_stars <- num_new_stars - 1
                                }
                                new_stars <- paste(rep("*", num_new_stars), collapse = "")
                                
                                message(new_stars, appendLF = FALSE)
                                flush.console()
                                curr_stars <<- num_stars
                              }
                              if (curr_stars >= max_stars) {
                                # The terminal "|" character that appears instead of a *
                                message("|")
                              }
                            }
                          }
                        )
)

#Binarize Sparse Matrix
binarizeMat <- function(mat){
	mat@x[mat@x > 0] <- 1
	mat
}

#'  LISI analysis 
#' @param  lisi_label1 The batch effect variable of interest such as 'batch' and 'donor'. Default 'donor'
#' @param  lisi_label2 The variable name that represents ground truth or high AUC score cell type'. Default 'CellType' 
#' @param  reference.celltypes.rds.file The RDS file name that contains cell type information.
#' @param  integrative.umap.rds.files The vector of file names that contain UMAP coordinate information generated by different integrative methods.
#' @param  integrative.method.names The vector of integrative method names represented in the same order as in 'integrative.umap.rds.files'.
#' @param  point.size The size of data points in the plot.
#' @param  IsShowPlot If TRUE, the function will show the plot. If not, it will return the lisi's scores.
#' @export 
#'

runLISI <- function(lisi_label1="donor",lisi_label2="CellType",
                    reference.celltypes.rds.file="",
                    integrative.umap.rds.files= "",
                    integrative.method.names="",
                    point.size=5,IsShowPlot=TRUE){
  
  reference_f<-readRDS(reference.celltypes)
  reference_f<-subset(reference_f,CellType !="")
  reference_f$Cell<-rownames(reference_f)
  lisi.score<-data.frame()
  
  for(i in 1:length(integrative.umaps)){
    my.umap<-readRDS(integrative.umaps[i])
    my.name<-method.names[i]
    my.dat<-merge(reference_f,my.umap,all.x = T)
    x.umap<-my.dat[,c("UMAP1","UMAP2")]
    my_lisi <- lisi::compute_lisi(x.umap,my.dat,c(lisi_label1,lisi_label2))
    z<-t(colMeans(my_lisi))
    rownames(z)<-my.name
    lisi.score<-rbind(lisi.score,z)
  }
  z<-lisi.score
  colnames(z)<-c("iLISI","cLISI")
  z$method=rownames(z)
  
  if(IsShowPlot==TRUE){
    ggplot(z,aes(x=cLISI,y=iLISI,color = method,shape=method)) +
      geom_point(size=point.size) +theme_bw()
  }else{
    return(z)
  }
}

#'  kBET analysis 
#' @param  covariate_variable_name The batch effect variable of interest such as 'batch' and 'donor'. Default 'donor'
#' @param  reference.celltypes.rds.file The RDS file name that contains cell type information.
#' @param  integrative.umap.rds.files The vector of file names that contain UMAP coordinate information generated by different integrative methods.
#' @param  integrative.method.names The vector of integrative method names represented in the same order as in 'integrative.umap.rds.files'.
#' @param  n.sample The downsample size of data points used in kBET analysis.
#' @export 
#'

runKBET <- function(covariate_variable_name="donor",
                    reference.celltypes.rds.file="",
                    integrative.umap.rds.files= "",
                    integrative.method.names="",
                    n.sample=1000){
  
  reference_f<-readRDS(reference.celltypes)
  reference_f<-subset(reference_f,CellType !="")
  reference_f$Cell<-rownames(reference_f)
  ###get the ground truth 'CellType' as cluster###
  clusters <- as.character(reference_f$CellType)
  ##############################################
  kBET_result_all<-data.frame()
  ####Iterate to each integrative method  RDS file###
  for(i in 1:length(integrative.umaps)){
    
    my.umap<-readRDS(integrative.umaps[i])
    my.name<-method.names[i]
    my.dat<-merge(reference_f,my.umap,all.x = T)
    
    kBET_result_list <- list()
    ####Iterate to each 'CellType'##############
    for (cluster_level in unique(clusters)){
      
      my.f<-subset(my.dat,CellType==cluster_level)
      my.umap<-as.matrix(my.f[,c("UMAP1","UMAP2")])
      my.batch<-as.character(my.f[,c(covariate_variable_name)])
      
      if(nrow(my.f)  > n.sample){
        s.index<-sample(nrow(my.f),size = n.sample,replace=FALSE)
        batch_tmp<-my.batch[s.index]
        data_tmp <-my.umap[s.index,]
      }else{
        batch_tmp<-my.batch
        data_tmp <-my.umap
      }
      kBET_tmp <- kBET(df=data_tmp, batch=batch_tmp, plot=FALSE)
      kBET_result_list[[cluster_level]] <- kBET_tmp
      kBET_result_list[[cluster_level]] <- mean(kBET_tmp$stats$kBET.observed)
      out.p<-paste(my.name,":",cluster_level,"is done!",sept="")
      print(out.p)
    }
    
    kBET_result <- melt(kBET_result_list)
    colnames(kBET_result) <- c("RejectRate","Cluster")
    kBET_result$AcceptanceRate <- 1 - kBET_result$RejectRate
    kBET_result$Method<-my.name
    kBET_result_all<-rbind(kBET_result_all,kBET_result)
  }
  return(kBET_result_all)
}
