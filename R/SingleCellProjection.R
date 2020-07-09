runLSI <- function(object,use.components=50, binarize = FALSE, use.regressout.data=FALSE){
	
	set.seed(1)
	
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	cells.info<-get_cells_annotation(object)
	cells.used<-subset(cells.info,IsPassed==TRUE)
	##############################
	genes_info<-get_genes_metadata(object)
	selected.genes<-subset(genes_info,IsVarGenes==TRUE)
	##############################
	if(binarize==TRUE){
		my.select.data<-get_umi_count(object)
		mat<-my.select.data[rownames(my.select.data) %in% rownames(selected.genes),as.character(cells.used$Cell)]
		message(paste0("Binarizing matrix..."))
		mat@x[mat@x > 0] <- 1 
	}else if(binarize==FALSE){
		if(use.regressout.data==T){
			my.select.data<-get_regressout_log_data(object)
			if(is.null(my.select.data)==T){
				stop("Please run 'remove_unwanted_confounders' function!")
			}
			#my.use.dat<-my.select.data[rownames(my.select.data) %in% rownames(selected.genes),as.character(cells.used$Cell)]
			S.m<-my.select.data[rownames(my.select.data) %in% rownames(selected.genes),as.character(cells.used$Cell)]
		}else if (use.regressout.data==F){
			my.select.data<-get_normalized_umi(object)
			#my.use.dat<-my.select.data[rownames(my.select.data) %in% rownames(selected.genes),as.character(cells.used$Cell)]
			S.m<-log1p(my.select.data[rownames(my.select.data) %in% rownames(selected.genes),as.character(cells.used$Cell)])
		}else{
			stop("Please set the parameter 'use.regressout.data' to T or F")
		}
		mat<-scale(S.m,center = FALSE,scale = TRUE)##scale with no center
	}
	
	#Calc RowSums and ColSums
	colSm <- Matrix::colSums(mat)
	rowSm <- Matrix::rowSums(mat)
	
	#Calc TF IDF
	message("Computing Term Frequency IDF...")
	freqs <- t(t(mat)/colSm)
	idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
	tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
	
	#Calc SVD then LSI
	message("Computing SVD using irlba...")
	svd <- irlba::irlba(tfidf, use.components, use.components)
	svdDiag <- matrix(0, nrow=use.components, ncol=use.components)
	diag(svdDiag) <- svd$d
	matSVD <- t(svdDiag %*% t(svd$v))
	rownames(matSVD) <- colnames(mat)
	colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
	###########################################
	#Return Object
	out <- list(
			matSVD = matSVD, 
			rowSm = rowSm, 
			colSm = colSm, 
			svd = svd, 
			binarize = binarize, 
			use.components = use.components,
			date = Sys.Date(),
			seed = 1)
	
	get_lsi.result(object)<-out
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("LSI analysis is done!.")
}