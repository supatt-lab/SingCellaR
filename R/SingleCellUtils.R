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
    if(metric=="euclidean"){
      res <- annoy_euclidean_nns(index_file,Xmat,n.neighbors, n.neighbors*100, grain_size = 1,verbose = FALSE)
    }
    if(metric=="cosine"){
      res <- annoy_cosine_nns(index_file,Xmat,n.neighbors, n.neighbors*100, grain_size = 1,verbose = FALSE)
    }
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
