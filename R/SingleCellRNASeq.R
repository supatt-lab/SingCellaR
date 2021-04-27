#' Load 10X sparse data matrices provided by the Cell Ranger software
#' @param  object The SingCellaR object.
#' @param  cellranger.version is the version of the Cell Ranger software used to generate input matrices.
#' @param  isMultiFeatures is logical
#' @param  selectedFeature SingCellaR supports two features 'Gene Expression' and 'Antibody Capture'
#' @export
#' @importFrom Matrix readMM
load_matrices_from_cellranger <- function(object,cellranger.version=3,
		isMultiFeatures=FALSE,selectedFeature="Gene Expression"){
	
	objName <- deparse(substitute(object))
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	data.dir<-object@dir_path_10x_matrix
	
	if (! dir.exists(data.dir)){
		stop("Directory provided does not exist")
	}
	input.files<-list.files(data.dir)
	if (length(input.files) !=3){
		stop("Required three files to be shown in the folder 'barcodes.tsv', 'genes.tsv or features.tsv', and 'matrix.mtx'")
	}
	if(cellranger.version==2){
		barcode.file <- paste(data.dir,"barcodes.tsv",sep="/")
		gene.file    <- paste(data.dir,"genes.tsv",sep="/")
		matrix.file  <- paste(data.dir,"matrix.mtx",sep="/")
		#if (!file.exists(barcode.file)){
		#  stop("Missing barcodes.tsv")
		#}
		#if (! file.exists(gene.file)){
		#  stop("Missing genes.tsv or features.tsv")
		#}
		#if (! file.exists(matrix.file )){
		#  stop("Missing matrix.mtx")
		#}
	}else if(cellranger.version >=3){
		barcode.file <- paste(data.dir,"barcodes.tsv.gz",sep="/")
		gene.file    <- paste(data.dir,"features.tsv.gz",sep="/")
		matrix.file  <- paste(data.dir,"matrix.mtx.gz",sep="/")
		#if (!file.exists(barcode.file)){
		#  stop("Missing barcodes.tsv.gz")
		#}
		#if (! file.exists(gene.file)){
		#  stop("Missing features.tsv.gz")
		#}
		#if (! file.exists(matrix.file )){
		#  stop("Missing matrix.mtx.gz")
		#}
	}else{
		stop("Required the version number of the cellranger2 or 3")
	}
	
	if(isMultiFeatures==FALSE){
		
		mat <- readMM(file = matrix.file)
		gene_info <- read.delim(gene.file, stringsAsFactors = FALSE,
				sep = "\t", header = FALSE)
		
		if (dim(mat)[1] != length(gene_info[, 1])) {
				stop("Mismatch dimension between gene file and the matrix file")
		}else {
				rownames(mat) <- make.unique(gene_info[, 2])
		}
		barcodes <- read.delim(barcode.file, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
		if (dim(mat)[2] != length(barcodes[, 1])) {
			stop("Mismatch dimension between barcodes file and the matrix file")
		}else {
			colnames(mat) <- paste(barcodes[, 1],"_",object@sample_uniq_id,sep="")
		}
		#my.counts = DelayedArray(mat)
		#gene.expressing.cells<-DelayedArray::rowSums(mat)
		#gene.expressing.cells<-Matrix::rowSums(mat)
		#sce<-SingleCellExperiment(assays = list(counts = mat[gene.expressing.cells > 0,]))
		sce<-SingleCellExperiment(assays = list(counts = mat))
		if(class(object)[1]=="SingCellaR"){
			object <- new("SingCellaR", sce, dir_path_10x_matrix=data.dir,sample_uniq_id=object@sample_uniq_id)
		}else if(class(object)[1]=="SingCellaR_Int"){
			object <- new("SingCellaR_Int", sce, dir_path_10x_matrix=data.dir,sample_uniq_id=object@sample_uniq_id)
		}else{
			stop(paste("There is no initiated class!"))
		}
	}else if(isMultiFeatures==TRUE){
		
		if(selectedFeature %in% c("Gene Expression","Antibody Capture") == FALSE){
			stop("Please input selected feature 'Gene Expression' or 'Antibody Capture'")
		}
		mat <- readMM(file = matrix.file)
		gene_info <- read.delim(gene.file, stringsAsFactors = FALSE,
				sep = "\t", header = FALSE)
		#features<-unique(gene_info[,3])
		my.index<-which(gene_info[,3]==selectedFeature)
		my.mat<-mat[my.index,]
		my.genes<-gene_info[,2][my.index]
		
		if (dim(my.mat)[1] != length(my.genes)) {
			stop("Mismatch dimension between gene file and the matrix file")
		}else {
			rownames(my.mat) <- make.unique(my.genes)
		}
		barcodes <- read.delim(barcode.file, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
		if (dim(my.mat)[2] != length(barcodes[, 1])) {
			stop("Mismatch dimension between barcodes file and the matrix file")
		}else {
			colnames(my.mat) <- paste(barcodes[, 1],"_",object@sample_uniq_id,sep="")
		} 
		sce<-SingleCellExperiment(assays = list(counts = my.mat))
		if(class(object)[1]=="SingCellaR"){
			object <- new("SingCellaR", sce, dir_path_10x_matrix=data.dir,sample_uniq_id=object@sample_uniq_id)
		}else if(class(object)[1]=="SingCellaR_Int"){
			object <- new("SingCellaR_Int", sce, dir_path_10x_matrix=data.dir,sample_uniq_id=object@sample_uniq_id)
		}else{
			stop(paste("There is no initiated class!"))
		}
	}
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("The sparse matrix is created.")
}


load_custom_MEX_matrices <- function(object,barcode.file="barcodes.tsv.gz",
		gene.file="features.tsv.gz",matrix.file="matrix.mtx.gz",
		gene_name_column=1){
	
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	data.dir<-object@dir_path_10x_matrix
	if (! dir.exists(data.dir)){
		stop("Directory provided does not exist")
	}
	barcode.file <- paste(data.dir,barcode.file,sep="/")
	gene.file    <- paste(data.dir,gene.file,sep="/")
	matrix.file  <- paste(data.dir,matrix.file,sep="/")
	
	mat <- readMM(file = matrix.file)
	gene_info <- read.delim(gene.file, stringsAsFactors = FALSE,
			sep = "\t", header = FALSE)
	
	if (dim(mat)[1] != length(gene_info[, 1])) {
		stop("Mismatch dimension between gene file and the matrix file")
	}else {
		rownames(mat) <- make.unique(gene_info[, gene_name_column])
	}
	barcodes <- read.delim(barcode.file, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
	if (dim(mat)[2] != length(barcodes[, 1])) {
		stop("Mismatch dimension between barcodes file and the matrix file")
	}else {
		colnames(mat) <- paste(barcodes[, 1],"_",object@sample_uniq_id,sep="")
	}
	sce<-SingleCellExperiment(assays = list(counts = mat))
	if(class(object)[1]=="SingCellaR"){
		object <- new("SingCellaR", sce, dir_path_10x_matrix=data.dir,sample_uniq_id=object@sample_uniq_id)
	}else if(class(object)[1]=="SingCellaR_Int"){
		object <- new("SingCellaR_Int", sce, dir_path_10x_matrix=data.dir,sample_uniq_id=object@sample_uniq_id)
	}else{
		stop(paste("There is no initiated class!"))
	}
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("The sparse matrix is created.")
}

#' Load gene expression from a file
#' @param  object The SingCellaR object.
#' @param  input_file The input gene expression file name.
#' @param  isTargetSeq is the logical input. default = TRUE
#' @param  sep The field separator character of the gene expression table.
#' @export 

load_gene_expression_from_a_file <- function(object,input_file="",isTargetSeq=T,sep=","){
  
  objName <- deparse(substitute(object))
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  #######Input file###########################
  if(isTargetSeq==T){
    input_file<-object@GenesExpressionMatrixFile
  }else if(input_file !=""){
    input_file<-input_file
  }else{
    stop("Requires input file!")
  }
  ############################################
  expM<-fread(input_file,sep = sep)
  g.names<-as.data.frame(expM[,1])
  g.names<-g.names[,1]
  expM<-expM[,-c(1)]
  expM.matrix<-as.matrix(expM)
  rownames(expM.matrix)<-make.unique(g.names)
  #######convert matrix to sparse matrix#####
  Sp <- Matrix(expM.matrix,sparse = T)
  sce<-SingleCellExperiment(assays = list(counts = Sp))
  ###########################################
  if(isTargetSeq==T){
    object <- new("TargetSeq", sce,GenesExpressionMatrixFile=object@GenesExpressionMatrixFile,
                  CellsMetaDataFile=object@CellsMetaDataFile)
  }else{
    object <- new("SingCellaR", sce)
  }
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("The sparse matrix is created.")
}

#' Subsample cells from the whole data set
#' @param  object The SingCellaR object.
#' @param  n.subsample The number of subsample cells to be selected.
#' @param  n.seed The seed number for the subsampling. This number is used for reproducing the same subsampling results.
#' @export
#' 

subsample_cells <- function(object,n.subsample=10000,n.seed = 123){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  ######################
  set.seed(n.seed)
  ######################
  s.index<-sample(1:nrow(object@meta.data), n.subsample, replace = FALSE)
  ######################
  object<-object[,s.index]
  s.meta.data<-object@meta.data[s.index,]
  object@meta.data<-s.meta.data
  ######################
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Downstream analysis will be performed using the subsampling cells!.")
  
}

#' Process all cells annotation
#' @param  object The SingCellaR object.
#' @param  mitochondiral_genes_start_with The unique prefix alphabets for mitocondrial genes such as 'MT-' for human and 'mt-' for mouse genes.
#' @export
#' @importFrom Matrix colSums
#' 
process_cells_annotation <- function(object,mitochondiral_genes_start_with="MT-"){

	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
  	umi.dat<-counts(object)
	##get meta data
	meta.data<-data.frame(Cell=colnames(umi.dat))
	meta.data$sampleID<-""
	cell.splitted<-strsplit(as.character(meta.data$Cell), "-", fixed = TRUE)
	ids <-do.call(rbind,cell.splitted)
	if(ncol(ids) > 1){
	  meta.data$sampleID<-ids[,2]
	}else{
	  meta.data$sampleID<-""
	}
	#meta.data$UMI_count<-DelayedArray::colSums(umi.dat)
	meta.data$UMI_count<-Matrix::colSums(umi.dat)
	##get number of detected genes per cell
	##n.genePerCell<-getNumberOfDetectedGenesPerCell(umi.dat,1)
	##meta.data$detectedGenesPerCell<-as.numeric(n.genePerCell$num_genes_detected)
	##binary.mat<-umi.dat
	##binary.mat[binary.mat > 1]<-1
	##meta.data$detectedGenesPerCell<-getDetectedGenesPerCell(as(umi.dat, "sparseMatrix"),1)
	meta.data$detectedGenesPerCell<-getDetectedGenesPerCell(umi.dat,1)
	
	##check % of mitocondrial genes
	searching.string<-paste("^",mitochondiral_genes_start_with,sep="")
	mito.genes <- grep(searching.string, rownames(umi.dat), value = T)
	print("List of mitochondrial genes:")
	print(mito.genes)
	#percent.mito.info<-data.frame(percent.mito=DelayedArray::colSums(umi.dat[mito.genes, ])/DelayedArray::colSums(umi.dat))
	percent.mito.info<-data.frame(percent.mito=Matrix::colSums(umi.dat[mito.genes, ])/Matrix::colSums(umi.dat))
	meta.data$percent_mito<-percent.mito.info$percent.mito*100
	##
	get_cells_annotation(object)<-meta.data
	##
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("The meta data is processed.")
	
}

#' Target-Seq processing for cell annotation
#' @param  object The TargetSeq object.
#' @param  mitochondiral_genes_start_with The unique prefix alphabets for mitocondrial genes such as 'MT-' for human and 'mt-' for mouse genes.
#' @param  ERCC_genes_start_with The unique prefix alphabets for ERCC genes.
#' @export
#' @importFrom Matrix colSums

TargetSeq_process_cells_annotation <- function(object,mitochondiral_genes_start_with="MT-",ERCC_genes_start_with="ERCC-"){
  
  objName <- deparse(substitute(object))
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  umi.dat<-counts(object)
  ##get meta data
  meta.data<-data.frame(Cell=colnames(umi.dat))
  meta.data$sampleID<-""
  cell.splitted<-strsplit(as.character(meta.data$Cell), "-", fixed = TRUE)
  ids <-do.call(rbind,cell.splitted)
  if(ncol(ids) > 1){
    meta.data$sampleID<-ids[,2]
  }else{
    meta.data$sampleID<-""
  }
  meta.data$UMI_count<-Matrix::colSums(umi.dat)
  #################################################
  normalized.umi<-(t(t(umi.dat)/meta.data$UMI_count))*round(mean(meta.data$UMI_count))
  #################################################
  meta.data$detectedGenesPerCell<-getDetectedGenesPerCell(normalized.umi,1)
  ##check % of mitocondrial genes
  searching.string<-paste("^",mitochondiral_genes_start_with,sep="")
  mito.genes <- grep(searching.string, rownames(umi.dat), value = T)
  print("List of mitochondrial genes:")
  print(mito.genes)
  #percent.mito.info<-data.frame(percent.mito=DelayedArray::colSums(umi.dat[mito.genes, ])/DelayedArray::colSums(umi.dat))
  percent.mito.info<-data.frame(percent.mito=Matrix::colSums(umi.dat[mito.genes, ])/Matrix::colSums(umi.dat))
  meta.data$percent_mito<-percent.mito.info$percent.mito*100
  ##
  ##check % of ERCC
  if(ERCC_genes_start_with!=""){
    searching.string<-paste("^",ERCC_genes_start_with,sep="")
    ERCC.genes <- grep(searching.string, rownames(umi.dat), value = T)
    print("List of ERCC:")
    print(ERCC.genes)
    #percent.mito.info<-data.frame(percent.mito=DelayedArray::colSums(umi.dat[mito.genes, ])/DelayedArray::colSums(umi.dat))
    percent.ERCC.info<-data.frame(percent.ERCC=Matrix::colSums(umi.dat[ERCC.genes, ])/Matrix::colSums(umi.dat))
    meta.data$percent_ERCC<- percent.ERCC.info$percent.ERCC*100
  }
  get_cells_annotation(object)<-meta.data
  ##
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("The meta data is processed.")
}

#' Filter out cells and genes
#' @param  object The SingCellaR object.
#' @param  min_UMIs The minimum UMI cutoff.
#' @param  max_UMIs The maximum UMI cutoff.
#' @param  min_detected_genes The minimum number of genes detected per cell.
#' @param  max_detected_genes The maximum number of genes detected per cell.
#' @param  min_percent_mito The minimum percent of UMIs from the mitocondrial genes.
#' @param  max_percent_mito The maximum percent of UMIs from the mitocondrial genes.
#' @param  genes_with_expressing_cells The cutoff for genes expressing in at least a number of cells.
#' @export
#' @importFrom Matrix rowSums

filter_cells_and_genes <- function(object,min_UMIs=1000,max_UMIs=30000,min_detected_genes=1000,max_detected_genes=8000,
                                   min_percent_mito=0,max_percent_mito=10,
                                   genes_with_expressing_cells=10){
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  meta.data<-get_cells_annotation(object)
  meta.data$IsPassed<-FALSE
  filtered.cells<-subset(meta.data,UMI_count >=min_UMIs & UMI_count <=max_UMIs)
  filtered.cells<-subset(filtered.cells,detectedGenesPerCell >=min_detected_genes & detectedGenesPerCell <=max_detected_genes)
  filtered.cells<-subset(filtered.cells,percent_mito >=min_percent_mito & percent_mito <=max_percent_mito)
  f.cells<-as.character(filtered.cells$Cell)
  meta.data$IsPassed[meta.data$Cell %in% f.cells]<-TRUE
  n.cells<-nrow(meta.data)
  n.filter<-nrow(meta.data[meta.data$IsPassed=="FALSE",])
  ###
  #umi.dat<-get_umi_count(object)
  #binary.mat<-umi.dat
  #binary.mat[binary.mat > 1]<-1
  binary.mat<-makeBinaryMatrix(get_umi_count(object),1)
  ###
  CountCellPerGene<-data.frame(num_detected_cells=Matrix::rowSums(binary.mat))
  rownames(CountCellPerGene)<-rownames(get_umi_count(object))
  CountCellPerGene$IsExpress<-FALSE
  CountCellPerGene$IsExpress[CountCellPerGene$num_detected_cells >=genes_with_expressing_cells]<-TRUE
  ##Update the metadata
  get_cells_annotation(object)<-meta.data
  get_genes_metadata(object)<-CountCellPerGene
  ##Update the object
  #cell.used<-subset(meta.data,IsPassed==TRUE)
  #object<-object[,as.character(cell.used$Cell)]
  ###################
  print("The cells and genes metadata are updated by adding the filtering status.")
  print(paste(n.filter,"/",n.cells," cells will be filtered out from the downstream analyses!.",sep=""))
  #print("The count matrix in the object is updated due to filtering out cells!.")
  #print(object)
  assign(objName,object,envir=parent.frame())
  invisible(1)
}

#' Loading cell metdata from TARGET-Seq analysis
#' @param  object The SingCellaR object.
#' @param  sep  The field separator character.
#' @param  a_column_of_the_cell_names The column that contains cell names.
#' @export
#' 
TargetSeq_load_cell_metadata <- function(object,sep="\t",
                                         a_column_of_the_cell_names=1){
  
  objName <- deparse(substitute(object))
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(a_column_of_the_cell_names==""){
    stop("Please input the number of column for the cell names!")
  }
  #if(a_column_of_the_cell_QC==""){
  #  stop("Please input the number of column for the cell QC")
  #}
  #######Input file###########################
  if(object@CellsMetaDataFile==""){
    stop("Please input the cell metadata file in the TargetSeq's object")
  }else{
    input_file<-object@CellsMetaDataFile
    CellsMetaData<-data.table::fread(input_file,sep = sep)
    CellsMetaData<-as.data.frame(CellsMetaData)
    ####check QC column###
    #CellsMetaData<-CellsMetaData[CellsMetaData[,a_column_of_the_cell_QC]==TRUE,]
  }
  ############################################
  ori.meta.data<-get_cells_annotation(object)
  ori.colnames<-colnames(ori.meta.data)
  ############################################
  cellMetaData.colnames<-colnames(CellsMetaData)
  matched.index<-which(cellMetaData.colnames %in% ori.colnames)
  if(length(matched.index) > 0){
    stop("The column names in the cell meta data file must not contain 'Cell','sampleID',
         'UMI_count','detectedGenesPerCell','percent_mito', 'percent_ERCC', and 'IsPassed'")
  }else{
    column_of_cell_name<-colnames(CellsMetaData[a_column_of_the_cell_names])
    updated.cell.metadata<-merge(ori.meta.data,CellsMetaData,by.x="Cell",by.y=column_of_cell_name)
  }
  ##Update the metadata
  get_cells_annotation(object)<-updated.cell.metadata
  ###################
  print("The cell metadata is updated!")
  ####################
  assign(objName,object,envir=parent.frame())
  invisible(1)
}


#' TargetSeq processing for filtering out cells and genes
#' @param  object TargetSeq object.
#' @param  min_Reads The minimum number of reads cutoff.
#' @param  min_detected_genes The minimum number of genes detected per cell.
#' @param  min_percent_mito The minimum percent of reads from mitocondrial genes.
#' @param  max_percent_mito The maximum percent of reads from mitocondrial genes.
#' @param  min_percent_ERCC The minimum percent of reads from ERCC.
#' @param  max_percent_ERCC The maximum percent of reads from ERCC.
#' @param  genes_with_expressing_cells The cutoff for genes expressing in expressing cells.
#' @export
#' @importFrom Matrix rowSums


TargetSeq_filter_cells_and_genes <- function(object,min_Reads=2000,min_detected_genes=500,
                                             min_percent_mito=0,max_percent_mito=10,min_percent_ERCC=0,
                                             max_percent_ERCC=50,genes_with_expressing_cells=10){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  meta.data<-get_cells_annotation(object)
  meta.data$IsPassed<-FALSE
  filtered.cells<-subset(meta.data,UMI_count >=min_Reads)
  filtered.cells<-subset(filtered.cells,detectedGenesPerCell >=min_detected_genes)
  filtered.cells<-subset(filtered.cells,percent_mito >=min_percent_mito & percent_mito <=max_percent_mito)
  filtered.cells<-subset(filtered.cells,percent_ERCC >=min_percent_ERCC & percent_ERCC <=max_percent_ERCC)
  f.cells<-as.character(filtered.cells$Cell)
  meta.data$IsPassed[meta.data$Cell %in% f.cells]<-TRUE
  n.cells<-nrow(meta.data)
  n.filter<-nrow(meta.data[meta.data$IsPassed=="FALSE",])
  ###
  umi.dat<-get_umi_count(object)
  binary.mat<-umi.dat
  binary.mat[binary.mat > 1]<-1
  ###
  CountCellPerGene<-data.frame(num_detected_cells=Matrix::rowSums(binary.mat))
  CountCellPerGene$IsExpress<-FALSE
  CountCellPerGene$IsExpress[CountCellPerGene$num_detected_cells >=genes_with_expressing_cells]<-TRUE
  ##Update the metadata
  get_cells_annotation(object)<-meta.data
  get_genes_metadata(object)<-CountCellPerGene
  ##Update the object
  #cell.used<-subset(meta.data,IsPassed==TRUE)
  #object<-object[,as.character(cell.used$Cell)]
  ###################
  print("The cells and genes metadata are updated by adding the filtering status.")
  print(paste(n.filter,"/",n.cells," cells will be filtered out from the downstream analyses!.",sep=""))
  #print("The count matrix in the object is updated due to filtering out cells!.")
  #print(object)
  assign(objName,object,envir=parent.frame())
  invisible(1)
}


#' Normalisation
#' @param  object The SingCellaR object.
#' @param  scale.factor The normalised scale factor which is usually set up at 10000.
#' @param  use.scaled.factor default is TRUE if it is FALSE, the function will use the mean of library size as the scale factor.
#' @export
#' @importFrom Matrix colSums

normalize_UMIs <- function(object,scale.factor = 1e4,use.scaled.factor=T){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  #######################
  umi.dat<-get_umi_count(object)
  #######################
  if(use.scaled.factor==T){
    
    if(ncol(umi.dat) < 50000){
      totalUMI.lib<-Matrix::colSums(umi.dat)
      normalized.umi<-(t(t(umi.dat)/totalUMI.lib))*scale.factor
    }else{
      normalized.umi<-count_divided_by_libsize(umi.dat)*scale.factor
      rownames(normalized.umi)<-rownames(umi.dat)
      colnames(normalized.umi)<-colnames(umi.dat)
    }
    assay(object, "normalized.umi")<-normalized.umi
  }else{
    if(ncol(umi.dat) < 50000){
      totalUMI.lib<-Matrix::colSums(umi.dat)
      normalized.umi<-(t(t(umi.dat)/totalUMI.lib))*round(mean(totalUMI.lib))
    }else{
      totalUMI.lib<-Matrix::colSums(umi.dat)
      normalized.umi<-count_divided_by_libsize(umi.dat)*round(mean(totalUMI.lib))
      rownames(normalized.umi)<-rownames(umi.dat)
      colnames(normalized.umi)<-colnames(umi.dat)
    }
    assay(object, "normalized.umi")<-normalized.umi
  }
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Normalization is completed!.")
}

normalize_ADTs <- function(object,method="CLR"){
	
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	#######################
	data<-counts(object)
	#######################
	if(method=="CLR"){
		
		clr_function = function(x) {
			return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
		}
		myapply <- ifelse(test = T, yes = pbapply, no = apply)
		norm.data <- myapply(X = data,FUN = clr_function,MARGIN=1)
		assay(object, "normalized.umi")<-t(norm.data)
	}
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("Normalization is completed!.")
}


#' Add the cell cycle gene scores into the cell metadata
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file that contains singature gene sets.
#' @param  gene_sets The vector of gene set names.
#' @export
#' @importFrom Matrix colMeans

add_cell_cycle_genes_score<-function(object,gmt.file=c(),gene_sets=c("G2M_Core","S_phase_Core")){
	
	objName <- deparse(substitute(object))
	
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(gmt.file)==0){
		stop("Required the path to GMT file!")
	}
	if(length(gene_sets)==0){
		stop("Required gene set names!")
	}
	#####
	signature.sets<-get_gmtGeneSets(gmt.file)
	#####
	check.1<-match(gene_sets,names(signature.sets))
	check.2<-sum(length(which(is.na(check.1))))
	if(check.2 > 0){
		stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
	}
	normalized.umi<-get_normalized_umi(object)
	#####
	Genes.score<-data.frame(score=rep(0,ncol(normalized.umi)))
	for(i in 1:length(gene_sets)){
		set.name<-gene_sets[i]
		genes.x<-signature.sets[[set.name]]
		exprs<-as.matrix(normalized.umi[rownames(normalized.umi) %in% genes.x,])
		MM<-Matrix::colMeans(exprs)
		MM.f<-data.frame(score=MM)
		colnames(MM.f)<-set.name
		Genes.score<-cbind(Genes.score,MM.f)
	}
	cells.anno<-get_cells_annotation(object)
	r.index<-which(colnames(cells.anno) %in% gene_sets)
	if(length(r.index) > 0){
		cells.anno<-cells.anno[,-c(r.index)]
		get_cells_annotation(object)<-cbind(cells.anno,Genes.score[-c(1)])
	}else{
		get_cells_annotation(object)<-cbind(get_cells_annotation(object),Genes.score[-c(1)])
	}
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("The cell meta data is updated!.")
}

#' Removing unwanted confounders using limma
#' @param  object The SingCellaR object.
#' @param  residualModelFormulaStr The model for removing confounder variables.
#' @param  preserved_feature is a defined preserved variable such as a cell genotype. 
#' @param  block.gene.size is the number of genes in each processing block.
#' @export 
#' @importFrom Matrix sparse.model.matrix
#' @importFrom limma lmFit

remove_unwanted_confounders<-function(object,residualModelFormulaStr="~UMI_count+percent_mito",
                                      preserved_feature="",block.gene.size=2000){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  
  FM<-log1p(get_normalized_umi(object)[,as.character(cells.used$Cell)])
  
  if(preserved_feature==""){
    set.seed(1)
    X.model_mat <- Matrix::sparse.model.matrix(as.formula(residualModelFormulaStr),
                                               data = cells.used, drop.unused.levels = TRUE)
    step.size = block.gene.size
    ends <- c(seq(from=block.gene.size,to=nrow(FM),by=step.size),as.integer(nrow(FM)))
    starts <- seq(from=1,to=nrow(FM),by=step.size)
    window.genes<-data.frame(start=starts,end=ends)
    
    pb <- txtProgressBar(max = nrow(window.genes), style = 3)
    datalist = vector("list",nrow(window.genes))
    
    for(i in 1:nrow(window.genes)){
      Sys.sleep(0.5);
      x.start<-window.genes$start[i]
      x.end <-window.genes$end[i]
      fit <- limma::lmFit(FM[x.start:x.end,], X.model_mat)
      beta <- fit$coefficients[, -1, drop = FALSE]
      beta[is.na(beta)] <- 0
      
      FM2 <- FM[x.start:x.end,] - beta %*% t(X.model_mat[, -1])
      #FM2 <- FM2[!is.na(row.names(FM2)), ]
      datalist[[i]] <- as.big.matrix(as.matrix(FM2))
      setTxtProgressBar(pb, pb$getVal()+1)
    }
    object@regressout_matrix<-as.RowLinkedMatrix(datalist)
    close(pb)
  }else{
    
    set.seed(1)
    X.model_mat <- Matrix::sparse.model.matrix(as.formula(residualModelFormulaStr),data = cells.used, drop.unused.levels = TRUE)
    design<-Matrix::sparse.model.matrix(as.formula(preserved_feature),data = cells.used, drop.unused.levels = TRUE)
    
    step.size = block.gene.size
    ends <- c(seq(from=block.gene.size,to=nrow(FM),by=step.size),as.integer(nrow(FM)))
    starts <- seq(from=1,to=nrow(FM),by=step.size)
    window.genes<-data.frame(start=starts,end=ends)
    
    pb <- txtProgressBar(max = nrow(window.genes), style = 3)
    datalist = vector("list",nrow(window.genes))
    
    for(i in 1:nrow(window.genes)){
      
      Sys.sleep(0.5);
      x.start<-window.genes$start[i]
      x.end <-window.genes$end[i]
      fit <- limma::lmFit(FM[x.start:x.end,], cbind(design,X.model_mat))
      beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
      beta[is.na(beta)] <- 0
      
      FM2 <- FM[x.start:x.end,] - beta %*% t(X.model_mat)
      datalist[[i]] <- as.big.matrix(as.matrix(FM2))
      setTxtProgressBar(pb, pb$getVal()+1)
    }
    object@regressout_matrix<-as.RowLinkedMatrix(datalist)
    close(pb)
  }
  
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Removing unwanted sources of variation is done!")
}

#' Identified variable genes by fitting the GLM model using mean vs coefficient of variation
#' @param  object The SingCellaR object.
#' @param  mean_expr_cutoff The mean expression cutoff.
#' @param  disp_zscore_cutoff The dispersion of z-score cutoff. 
#' @param  quantile_genes_expr_for_fitting The quantile gene expression used for fitting the model.
#' @param  quantile_genes_cv2_for_fitting The quantile coefficient of variation used for fitting the model.
#' @export 
#'

get_variable_genes_by_fitting_GLM_model <- function(object,mean_expr_cutoff=0.1,disp_zscore_cutoff=0.1,
		quantile_genes_expr_for_fitting=0.20,quantile_genes_cv2_for_fitting=0.80){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	cells.info<-get_cells_annotation(object)
	cells.used<-subset(cells.info,IsPassed==TRUE)
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
	###
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
	#vars  <- matrixStats::rowVars(as.matrix(normalized.umi.used))
	###
	cv2 <- vars/means^2
	
	na.index<-which(is.na(cv2)==T)
	if(length(na.index) > 0){
		means<-means[-c(na.index)]
		cv2<-cv2[-c(na.index)]
	}
	###############################################
	####Get genes for fitting######################
	minMeanForFit <- unname(quantile( means,quantile_genes_expr_for_fitting))
	maxVarForFit <- unname(quantile(cv2,quantile_genes_cv2_for_fitting))
	###############################################
	genesForFit<-which(means >=minMeanForFit & cv2 <=maxVarForFit)
	print(paste("Using :",length(genesForFit)," genes for fitting the GLM model!",sep=""))
	################################################
	fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[genesForFit] ),cv2[genesForFit] )
	a0 <- unname( fit$coefficients["a0"] )
	a1 <- unname( fit$coefficients["a1tilde"])
	fit$coefficients
	
	afit <- a1/means+a0
	n.m<-data.frame(gene=names(means),means=means,cv2=cv2,varFit=afit,log.cv2=log(cv2),log.varFit=log(afit))
	y.m<-merge(n.m,selected.genes,all=T)
	y.m$percent_detect_cells<-y.m$num_detected_cells/nrow(cells.info)
	y.m<-y.m[order(y.m$means,decreasing=T),]
	y.m<-na.omit(y.m)
	###
	y.m$z<-(y.m$log.cv2-y.m$log.varFit)/sd(y.m$log.cv2-y.m$log.varFit)
	###
	y.m<-subset(y.m,cv2 > varFit)
	y.m<-subset(y.m,means >=mean_expr_cutoff & z>disp_zscore_cutoff)
	###
	y.m<-subset(y.m,percent_detect_cells < 0.95)
	#z.m<-subset(y.m,percent_detect_cells >= 0.90)
	###
	#if(nrow(z.m) > 10){
	#  y.m<-subset(y.m,percent_detect_cells < mean(z.m$percent_detect_cells))
	#}else{
	#  y.m<-subset(y.m,percent_detect_cells < 0.96)
	#}
	###
	var.genes<-as.character(y.m$gene)
	###
	genes_info$IsVarGenes<-FALSE
	genes_info$IsVarGenes[rownames(genes_info) %in% var.genes]<-TRUE
	###
	get_genes_metadata(object)<-genes_info
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print(paste("Identified :",length(var.genes)," variable genes",sep=""))
}

#' Identified variable genes by using the standard mean of gene expression and coefficient of variation
#' @param  object The SingCellaR object.
#' @param  mean_expr_cutoff The mean expression cutoff. Default 0.1
#' @param  max_of_mean_expr_cutoff The maximum number of gene expression cutoff. This is to filter out housekeeping and ribosomal genes. Default 5
#' @param  cv.max The maximum CV cutoff. Default 3
#' @param  cv.min The minimum CV cutoff. Default 0.5
#' @export 
#'


get_variable_genes_basic <- function(object,mean_expr_cutoff=0.1,max_of_mean_expr_cutoff=5,cv.max=3,cv.min=0.5){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	cells.info<-get_cells_annotation(object)
	cells.used<-subset(cells.info,IsPassed==TRUE)
	###
	normalized.umi<-get_normalized_umi(object)
	normalized.umi.used<-normalized.umi[,as.character(cells.used$Cell)]
	###
	genes_info<-get_genes_metadata(object)
	selected.genes<-subset(genes_info,IsExpress==TRUE)
	selected.genes$gene<-rownames(selected.genes)
	###
	normalized.umi.used<-normalized.umi.used[rownames(normalized.umi.used) %in% rownames(selected.genes),]
	###
	my.mean<-Matrix::rowMeans(normalized.umi.used)
	vars.S <- rowSds(as.matrix(normalized.umi.used))
	my.cv <- vars.S / my.mean
	###
	n.m<-data.frame(mean=my.mean,cv=my.cv)
	n.m<-n.m[order(n.m$mean,decreasing=T),]
	###
	n.m<-subset(n.m,mean >= mean_expr_cutoff & cv >= cv.min & cv <= cv.max)
	###
	var.genes<-as.character(rownames(n.m))
	###
	genes_info$IsVarGenes<-FALSE
	genes_info$IsVarGenes[rownames(genes_info) %in% var.genes]<-TRUE
	###
	get_genes_metadata(object)<-genes_info
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print(paste("Identified :",length(var.genes)," variable genes",sep=""))
}

#' Identified variable genes for the integrative samples. The gene list is combined from individual sample.
#' @param  object The SingCellaR object.
#' @param  mean_expr_cutoff The mean expression cutoff. Default 0.1
#' @param  disp_zscore_cutoff The dispersion of z-score cutoff. Default 0.1
#' @param  quantile_genes_expr_for_fitting The quantile gene expression used for fitting the model. Default 0.2
#' @param  quantile_genes_cv2_for_fitting The quantile coefficient of variation used for fitting the model. Default 0.8
#' @export 
#'

get_variable_genes_for_integrative_data_by_fitting_GLM_model <- function(object,mean_expr_cutoff=0.1,disp_zscore_cutoff=0.1,
                                                                         quantile_genes_expr_for_fitting=0.20,quantile_genes_cv2_for_fitting=0.80){
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  normalized.umi<-get_normalized_umi(object)
  normalized.umi.used<-normalized.umi[,as.character(cells.used$Cell)]
  rm(normalized.umi)
  ###use identified variable genes from each individual sample###
  variable.genes<-unique(unlist(object@Variable.genes))
  ######################################################
  genes_info<-get_genes_metadata(object)
  selected.genes<-genes_info[rownames(genes_info) %in% variable.genes,]
  selected.genes<-selected.genes[c("num_detected_cells")]
  selected.genes$gene=rownames(selected.genes)
  ###
  normalized.umi.used<-normalized.umi.used[rownames(normalized.umi.used) %in% variable.genes,]
  ###
  means <- Matrix::rowMeans(normalized.umi.used)
  ###
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
  #vars  <- matrixStats::rowVars(as.matrix(normalized.umi.used))
  ###
  cv2 <- vars/means^2
  
  na.index<-which(is.na(cv2)==T)
  if(length(na.index) > 0){
    means<-means[-c(na.index)]
    cv2<-cv2[-c(na.index)]
  }
  ###############################################
  ####Get genes for fitting######################
  minMeanForFit <- unname(quantile( means,quantile_genes_expr_for_fitting))
  maxVarForFit <- unname(quantile(cv2,quantile_genes_cv2_for_fitting))
  ###############################################
  genesForFit<-which(means >=minMeanForFit & cv2 <=maxVarForFit)
  print(paste("Using :",length(genesForFit)," genes for fitting the GLM model!",sep=""))
  ################################################
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[genesForFit] ),cv2[genesForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  fit$coefficients
  
  afit <- a1/means+a0
  n.m<-data.frame(gene=names(means),means=means,cv2=cv2,varFit=afit,log.cv2=log(cv2),log.varFit=log(afit))
  y.m<-merge(n.m,selected.genes,all=T)
  y.m$percent_detect_cells<-y.m$num_detected_cells/nrow(cells.info)
  y.m<-y.m[order(y.m$means,decreasing=T),]
  y.m<-na.omit(y.m)
  ###
  y.m$z<-(y.m$log.cv2-y.m$log.varFit)/sd(y.m$log.cv2-y.m$log.varFit)
  ###
  y.m<-subset(y.m,cv2 > varFit)
  y.m<-subset(y.m,means >=mean_expr_cutoff & z>disp_zscore_cutoff)
  ###
  y.m<-subset(y.m,percent_detect_cells < 0.95)
  #z.m<-subset(y.m,percent_detect_cells >= 0.90)
  ###
  #if(nrow(z.m) > 10){
  #  y.m<-subset(y.m,percent_detect_cells < mean(z.m$percent_detect_cells))
  #}else{
  #  y.m<-subset(y.m,percent_detect_cells < 0.96)
  #}
  ###
  var.genes<-as.character(y.m$gene)
  ###
  genes_info$IsVarGenes<-FALSE
  genes_info$IsVarGenes[rownames(genes_info) %in% var.genes]<-TRUE
  ###
  get_genes_metadata(object)<-genes_info
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print(paste("Identified :",length(var.genes)," variable genes",sep=""))
}

#' Remove unwanted genes from identified list of variable genes
#' @param  object The SingCellaR object.
#' @param  gmt.file The GMT file that containes the signature gene sets.
#' @param  removed_gene_sets The vector of gene sets to be removed from the list of variable genes.
#' @export 
#'

remove_unwanted_genes_from_variable_gene_set<-function(object,gmt.file=c(),removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene")){
  
  objName <- deparse(substitute(object))
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(gmt.file)==0){
    stop("Required the path to GMT file!")
  }
  if(length(removed_gene_sets)==0){
    stop("Required gene set names!")
  }
  #####
  signature.sets<-get_gmtGeneSets(gmt.file)
  #####
  check.1<-match(removed_gene_sets,names(signature.sets))
  check.2<-sum(length(which(is.na(check.1))))
  if(check.2 > 0){
    stop(paste(check.2," input gene names do not present in the GMT file.",sep=""))
  }
  r.genes<-c()
  for(i in 1:length(removed_gene_sets)){
    set.name<-removed_gene_sets[i]
    genes.x<-signature.sets[[set.name]]
    r.genes<-append(r.genes,genes.x)
  }
  r.genes.uniq<-unique(r.genes)
  ##
  genes_info<-get_genes_metadata(object)
  n.before<-nrow(subset(genes_info,IsVarGenes==T))
  genes_info$IsVarGenes[rownames(genes_info) %in%  r.genes.uniq]<-FALSE
  ###
  n.after<-nrow(subset(genes_info,IsVarGenes==T))
  n.diff<-n.before-n.after
  get_genes_metadata(object)<-genes_info
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print(paste(n.diff, "genes are removed from the variable gene set."),sep="")
}

#' PCA analysis
#' @param  object The SingCellaR object.
#' @param  use.components The number of principal components. Default 50
#' @param  regressout.data is logical, if TRUE, PCA uses the regressout data. If FALSE, PCA uses the normalized data without removing any umwanted source of variations.
#' @param  use.scanorama.integrative.matrix is the logical, if TRUE scannorama integrative matrix will be used. 
#' @export 
#'

runPCA <- function(object,use.components=50,use.regressout.data=T,use.scanorama.integrative.matrix=F){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	set.seed(seed = 1)
	cells.info<-get_cells_annotation(object)
	cells.used<-subset(cells.info,IsPassed==TRUE)
	##############################
	if(use.scanorama.integrative.matrix==TRUE){
	  print("The scanorama matrix will be used for PCA!")
	   S.m<-object@Scanorama.integrative.matrix
	}else{
	  genes_info<-get_genes_metadata(object)
	  selected.genes<-subset(genes_info,IsVarGenes==TRUE)
	  ##############################
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
	}
	###PCA using irlba ###########
	#S.m<-scale(my.use.dat,center = FALSE,scale = TRUE)##This line would be able to skip.
	#S.m<-scale(S.m,center = FALSE,scale = TRUE)
	#rm(my.use.dat)
	center.m <- Matrix::colMeans(t(S.m))
	##############################
	irlba.res <- irlba(A = t(S.m), nv = use.components,tol=1e-10,center = center.m)
	gene.loadings <- irlba.res$v	
	sdev <- irlba.res$d/sqrt(max(1, ncol(S.m) - 1))
	cell.embeddings <- irlba.res$u %*% diag(irlba.res$d)
	
	rownames(gene.loadings) <- rownames(S.m)
	colnames(gene.loadings) <- paste0('PC', 1:use.components)
	rownames(cell.embeddings) <- colnames(S.m)
	colnames(cell.embeddings) <- colnames(gene.loadings)
	
	res.pca <- list(gene.loadings = gene.loadings, x=cell.embeddings, sdev=sdev)
	
	get_pca.result(object)<-res.pca
	#assay(object,"regressedout.log.data")<-c()
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("PCA analysis is done!.")
}


#' NNMF analysis
#' @param  object The SingCellaR object.
#' @param  use.regressout.data is logical, if TRUE, NNMF uses the regressout data. If FALSE, NNMF uses the normalized data without removing any umwanted source of variations.
#' @param  use.scanorama.integrative.matrix is logical, if TRUE scannorama integrative matrix will be used. 
#' @param  n.threads The number of threads for running. Default 4
#' @param  max.iter The number of NNMF interations. Default 1000
#' @param  rel.tol The nnmf parameter. Default 1e-05
#' @export 
#'

runNNMF <- function(object,k=30,use.regressout.data=T,use.scanorama.integrative.matrix=F,n.threads=4,max.iter=1000,rel.tol=1e-05){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  set.seed(seed = 1)
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  ##############################
  if(use.scanorama.integrative.matrix==TRUE){
    print("The scanorama matrix will be used for NNMF!")
    S.m<-object@Scanorama.integrative.matrix
  }else{
    genes_info<-get_genes_metadata(object)
    selected.genes<-subset(genes_info,IsVarGenes==TRUE)
    ##############################
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
  }
  S.m<-scale(S.m,center = FALSE,scale = TRUE)##scale with no center
  #########RUN NNMF##########################
  decomp <- nnmf(S.m,k,n.threads = n.threads,max.iter = max.iter,rel.tol=rel.tol)
  get_nnmf.result(object)<-decomp
  ###########################################
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Nonnegative matrix factorization analysis is done!.")
}

#'  TSNE analysis
#' @param  object The SingCellaR object.
#' @param  dim_reduction_method The dimensional reduction method name that specifies the embedding matrix.
#' @param  useIntegrativeEmbeddings is logical, if TRUE the embedding matrix genereated from integrative analysis will be used.
#' @param  integrative_method The name of an integrative method.
#' @param  n.dims.use The number of dimensions as the input.
#' @param  n.dims.embed The number of output dimension. Default 2
#' @param  n.perplexity The TSNE perplexity. Default 30
#' @param  n.seed The number of set seed.
#' @export 
#'

runTSNE <- function(object,dim_reduction_method=c("pca","nnmf"),useIntegrativeEmbeddings=F,
                    integrative_method=c("combat","seurat","harmony","supervised_harmony"),
                    n.dims.use=50,n.dims.embed=2,n.perplexity=30, n.seed = 1){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	###
	if(n.dims.embed >3){
		stop("Please use 'n.dims.embed' = 2 or 3")
	}
	dim_reduction_method<-match.arg(dim_reduction_method)
	integrative_method<-match.arg(integrative_method)
	###
	cells.info<-get_cells_annotation(object)
	cells.used<-subset(cells.info,IsPassed==TRUE)
	###
	my.select.data<-get_normalized_umi(object)
	my.use.dat<-my.select.data[,as.character(cells.used$Cell)]
	###
	if(useIntegrativeEmbeddings==FALSE){
	  if(dim_reduction_method=="pca"){
	     res.pca<-get_pca.result(object)
	     my.pca<-as.matrix(res.pca$x[,1:n.dims.use])
	  }else if(dim_reduction_method=="nnmf"){
	     res.pca<-get_nnmf.result(object)
	     my.pca<-t(res.pca$H)[,1:n.dims.use]
	  }
	}else if(useIntegrativeEmbeddings==TRUE & integrative_method=="harmony"){
	  my.pca<-object@Harmony.embeddings[,1:n.dims.use]
	}else if(useIntegrativeEmbeddings==TRUE & integrative_method=="supervised_harmony"){
	    my.pca<-object@SupervisedHarmony.embeddings[,1:n.dims.use]
	}else if(useIntegrativeEmbeddings==TRUE & integrative_method=="seurat"){
	  my.pca<-object@Seurat.embeddings[,1:n.dims.use]
	}else if(useIntegrativeEmbeddings==TRUE & integrative_method=="combat"){
		my.pca<-object@Combat.embeddings[,1:n.dims.use]
	}
	###
	set.seed(seed = n.seed)
	tsne_out <- Rtsne(my.pca,dims = n.dims.embed,check_duplicates = FALSE,perplexity = n.perplexity)
	my.results<-data.frame(Cell=colnames(my.use.dat),tsne_out$Y)
	
	if(n.dims.embed==2){
		tsne.dim.names<-c("Cell","TSNE1","TSNE2")
	}else if(n.dims.embed==3){
		tsne.dim.names<-c("Cell","TSNE1","TSNE2","TSNE3")
	}else{
		stop("Please use 'n.dims.embed' = 2 or 3")
	}
	colnames(my.results)<-tsne.dim.names
	res.tsne<-merge(my.results,cells.used,by.x="Cell",by.y="Cell",all=T)
	get_tsne.result(object)<-res.tsne
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("TSNE analysis is done!.")
}

#'  UMAP analysis
#' @param  object The SingCellaR object.
#' @param  dim_reduction_method Dimensional reduction method name that specifies the embedding matrix.
#' @param  useIntegrativeEmbeddings is logical, if TRUE the embedding matrix genereated from integrative analysis will be used.
#' @param  integrative_method The name of an integrative method.
#' @param  umap_method The seleced UMAP analysis method.
#' @param  n.dims.use The number of dimensions as the input.
#' @param  n.neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Default 30
#' @param  n.seed The number of seed.
#' @param  uwot.metric uwot's distance metric. Default 'cosine'.
#' @param  uwot.spread  uwot spread parameter. Default 1.0
#' @param  uwot.set.op.mix.ratio uwot set.op.mix.ratio parameter. Default 1.0
#' @param  uwot.local.connectivity uwot local connectivity parameter. Default 1L
#' @param  uwot.repulsion.strength uwot repulsion.strength parameter. Default 1
#' @param  uwot.negative.sample.rate uwot negative.sample.rate parameter. Default 5
#' @param  uwot.a uwot a parameter.
#' @param  uwot.b uwot b parameter.
#' @param  uwot.metric.kwds uwot metric.kwds parameter.
#' @param  uwot.angular.rp.forest is the logical of uwot angular.rp.forest. 
#' @param  uwot.verbose is logical for uwot verbose parameter.
#' @param  uwot.save.model is the logical. If TRUE, the uwot model will be saved.
#' @param  uwot.save.file The file name of uwot model to be saved.
#' @export 
#'

runUMAP <- function(object,dim_reduction_method=c("pca","nnmf","lsi"),useIntegrativeEmbeddings=FALSE,integrative_method=c("combat","seurat","harmony","supervised_harmony"),
                    umap_method=c("uwot"),n.dims.use=50,n.neighbors=30,n.seed = 1,uwot.metric = 'cosine',
                    uwot.n.epochs = NULL,uwot.learning.rate = 1.0,uwot.min.dist = 0.25,uwot.spread = 1.0,uwot.set.op.mix.ratio = 1.0,
                    uwot.local.connectivity = 1L,uwot.repulsion.strength = 1,uwot.negative.sample.rate = 5,uwot.a = NULL,uwot.b = NULL,
                    uwot.metric.kwds = NULL,uwot.angular.rp.forest = FALSE,uwot.verbose = TRUE,uwot.save.model = FALSE,uwot.save.file="uwot.out.rdata"){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  dim_reduction_method<-match.arg(dim_reduction_method)
  integrative_method<-match.arg(integrative_method)
  umap_method <- match.arg(umap_method)
  cells.info  <-get_cells_annotation(object)
  cells.used  <-subset(cells.info,IsPassed==TRUE)
  ###
  if(useIntegrativeEmbeddings==FALSE){
    if(dim_reduction_method=="pca"){
      res.pca<-get_pca.result(object)
      my.pca<-as.matrix(res.pca$x[,1:n.dims.use])
    }else if(dim_reduction_method=="nnmf"){
      res.pca<-get_nnmf.result(object)
      my.pca<-t(res.pca$H)[,1:n.dims.use]
    }else if(dim_reduction_method=="lsi"){
      res.pca<-get_lsi.result(object)
      my.pca<-res.pca$matSVD[,1:n.dims.use]
    }
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="harmony"){
    my.pca<-object@Harmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="supervised_harmony"){
    my.pca<-object@SupervisedHarmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="seurat"){
    my.pca<-object@Seurat.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="combat"){
	  my.pca<-object@Combat.embeddings[,1:n.dims.use]
  }
  ###
  #if(umap_method=="umap-learn"){
  #  ###custom settings#############
  #  custom.settings = umap.defaults
  #  custom.settings$n_neighbors = n.neighbors
  #  custom.settings$min_dist= min_dist
  #  custom.settings$random_state= n.seed
  #  ###############################
  #  if(py_module_available("umap")==TRUE){
  #    umap_out <- umap(my.pca,config = custom.settings,method="umap-learn",metric=metric)
  #  }else{
  #    stop("The umap-learn python module is not installed!. To use the efficient umap, please install using pip ('pip install umap-learn').")
  #  }
  #  res.umap<-data.frame(Cell=rownames(umap_out$layout),umap_out$layout)
  #}else if(umap_method=="uwot"){
  if(umap_method=="uwot"){ 
    set.seed(seed = n.seed)
    umap.out <- uwot::umap(my.pca,
                           n_neighbors = as.integer(x = n.neighbors),
                           n_components = as.integer(x = 2L),
                           metric = uwot.metric,
                           n_epochs = uwot.n.epochs,
                           learning_rate = uwot.learning.rate,
                           min_dist = uwot.min.dist,
                           spread = uwot.spread,
                           set_op_mix_ratio = uwot.set.op.mix.ratio,
                           local_connectivity = uwot.local.connectivity,
                           repulsion_strength = uwot.repulsion.strength,
                           negative_sample_rate = uwot.negative.sample.rate,
                           a = uwot.a,
                           b = uwot.b,
                           verbose = uwot.verbose, 
                           ret_nn = TRUE,
                           ret_model = TRUE)
                           
    if(uwot.save.model==TRUE){
    	save_uwot(umap.out, uwot.save.file)
    }
    res.umap<-data.frame(Cell=rownames(my.pca),umap.out$embedding)
    
  }else{
    stop("Please choose the umap methods : uwot or umap-learn.")
  }
  colnames(res.umap)<-c("Cell","UMAP1","UMAP2")
  res.umap<-merge(cells.used,res.umap,by.x="Cell",by.y="Cell",all=T)
  get_umap.result(object)<-res.umap
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("UMAP analysis is done!.")
}

#'  Diffusion map analysis
#' @param  object The SingCellaR object.
#' @param  dim_reduction_method Dimensional reduction method name that specifies the embedding matrix.
#' @param  useIntegrativeEmbeddings is logical, if TRUE the embedding matrix genereated from integrative analysis will be used.
#' @param  integrative_method The name of an integrative method.
#' @param  n.dims.use The number of dimensions as the input.
#' @param  n.dims.embed is the number of output dimension. Default 3
#' @param  n.neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Default 30
#' @param  distance The distance metric.
#' @param  n.seed The number of seed.
#' @export 
#'

runDiffusionMap <- function(object,dim_reduction_method=c("pca","nnmf"),useIntegrativeEmbeddings=F,
                    integrative_method=c("combat","seurat","harmony","supervised_harmony"),
                    n.dims.use=50,n.dims.embed=3,n.neighbors=30,distance=c("euclidean","cosine"), 
                    n.seed = 1){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  ###
  dim_reduction_method<-match.arg(dim_reduction_method)
  integrative_method<-match.arg(integrative_method)
  distance<-match.arg(distance)
  ###
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  ###
  my.select.data<-get_normalized_umi(object)
  my.use.dat<-my.select.data[,as.character(cells.used$Cell)]
  ###
  if(useIntegrativeEmbeddings==FALSE){
    if(dim_reduction_method=="pca"){
      res.pca<-get_pca.result(object)
      my.pca<-as.matrix(res.pca$x[,1:n.dims.use])
    }else if(dim_reduction_method=="nnmf"){
      res.pca<-get_nnmf.result(object)
      my.pca<-t(res.pca$H)[,1:n.dims.use]
    }
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="harmony"){
    my.pca<-object@Harmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="supervised_harmony"){
    my.pca<-object@SupervisedHarmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="seurat"){
    my.pca<-object@Seurat.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="combat"){
    my.pca<-object@Combat.embeddings[,1:n.dims.use]
  }
  ###
  set.seed(seed = n.seed)
  print("This process will take time and requires large RAM depending on the number of cells in your analysis!")
  dfm_out <- DiffusionMap(my.pca,n_eigs = n.dims.embed,k=n.neighbors,distance = distance)
  my.results<-data.frame(Cell=colnames(my.use.dat),eigenvectors(dfm_out))
  ###
  res.dfm<-merge(my.results,cells.used,by.x="Cell",by.y="Cell",all=T)
  get_dfm.result(object)<-res.dfm
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("DiffusionMap analysis is done!.")
}

#'  ForceAtlas2  analysis
#' @param  object The SingCellaR object.
#' @param  useIntegrativeEmbeddings is logical, if TRUE the embedding matrix genereated from integrative analysis will be used.
#' @param  integrative_method The name of an integrative method.
#' @param  knn.metric Type of distance metric to use to find nearest neighbors.
#' @param  knn.n_trees The number of trees for building the annoy index. Default 50
#' @param  knn.n_threads The number of threads. Default 1
#' @param  n.dims.use The number of dimensions as the input. Default 30
#' @param  n.neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Default 30
#' @param  fa2_n_iter The FA2 number of iterations. Default 1000
#' @param  fa2_edgeWeightInfluence The FA2 edgeWeightInfluence parameter. Default 1
#' @param  fa2_barnesHutTheta The FA2 barnesHutTheta paramenter. Default 2
#' @param  fa2_scalingRatio The FA2 scalingRatio parameter. Default 1
#' @param  fa2_gravity The FA2 gravity parameter. Default 0.05
#' @param  fa2_jitterTolerance The FA2 jitterTolerance parameter. Default 1
#' @param  n.seed The number of set seed.
#' @export 
#'

runFA2_ForceDirectedGraph <- function(object,dim_reduction_method=c("pca","nnmf","lsi"),useIntegrativeEmbeddings=F,
                                      integrative_method=c("combat","seurat","harmony","supervised_harmony"),
                                      knn.metric=c("euclidean","cosine"),
                                      n.dims.use=30,n.neighbors=5,n.seed = 1,knn.n_trees=50,knn.n_threads =1,
                                      fa2_n_iter=1000,fa2_edgeWeightInfluence=1,fa2_barnesHutTheta=2, 
                                      fa2_scalingRatio=1, fa2_gravity=0.05, fa2_jitterTolerance=1){
  
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  dim_reduction_method<-match.arg(dim_reduction_method)
  integrative_method<-match.arg(integrative_method)
  knn.metric<-match.arg(knn.metric)
  ###
  if(useIntegrativeEmbeddings==FALSE){
    if(dim_reduction_method=="pca"){
      res.pca<-get_pca.result(object)
      my.pca<-as.matrix(res.pca$x[,1:n.dims.use])
    }else if(dim_reduction_method=="nnmf"){
      res.pca<-get_nnmf.result(object)
      my.pca<-t(res.pca$H)[,1:n.dims.use]
    }else if(dim_reduction_method=="lsi"){
      res.pca<-get_lsi.result(object)
      my.pca<-res.pca$matSVD[,1:n.dims.use]
    }
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="harmony"){
    my.pca<-object@Harmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="supervised_harmony"){
    my.pca<-object@SupervisedHarmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="seurat"){
    my.pca<-object@Seurat.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="combat"){
    my.pca<-object@Combat.embeddings[,1:n.dims.use]
  }
  
  knn <- identifyKNN(as.matrix(my.pca),n.neighbors=n.neighbors,
                     metric = knn.metric,n_trees = knn.n_trees,
                     n_threads = knn.n_threads)
  ###############################
  links<-compute_weights(knn)
  links <- links[links[,1]>0, ]
  ####make igraph#########
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  G.igraph <- graph.data.frame(relations, directed=FALSE)
  G.igraph <- igraph::simplify(G.igraph)
  ###############################
  print("Processing fa2..")
  ###############################
  if(py_module_available("fa2")==FALSE){
    stop("The fa2 python module is not installed!. Please install using pip ('pip install fa2')")
  }
  if(py_module_available("networkx")==FALSE){
    stop("The networkx python module is not installed!. Please install using pip ('pip install networkx')")
  }
  if(py_module_available("numpy")==FALSE){
    stop("The numpy python module is not installed!. Please install using pip ('pip install numpy')")
  }
  
  if(py_module_available("fa2")==TRUE){
    
    ###3 lines below were for reading python in the ".py" functions found in the "src" folder.
	### These 3 lines will be commented when apply to the package###
    python.fa2 <<- reticulate::import("fa2", delay_load=TRUE)
    python.nx <<- reticulate::import("networkx", delay_load=TRUE)
    python.np <<- reticulate::import("numpy", convert=FALSE,delay_load=TRUE)
    
    python.links <-python.np$array(links,dtype="object")
    #######################
    G<-python.nx$Graph()
    #G$add_nodes_from(nodes)
    G$add_weighted_edges_from(python.links)
    
    fa2.config = list(
      # Behavior alternatives
      outboundAttractionDistribution=F,  # Dissuade hubs
      linLogMode=F,  # NOT IMPLEMENTED
      adjustSizes=F,  # Prevent overlap (NOT IMPLEMENTED)
      edgeWeightInfluence=fa2_edgeWeightInfluence,
      
      # Performance
      jitterTolerance=fa2_jitterTolerance,  # Tolerance
      barnesHutOptimize=T,
      barnesHutTheta=fa2_barnesHutTheta,
      multiThreaded=F,  # NOT IMPLEMENTED
      # Tuning
      scalingRatio=fa2_scalingRatio,
      strongGravityMode=F,
      gravity=fa2_gravity,
      # Log
      verbose=T
    )
    
    FA2 = do.call(python.fa2$ForceAtlas2,fa2.config)
    py_set_seed(n.seed, disable_hash_randomization = TRUE)
    
    positions = FA2$forceatlas2_networkx_layout(G, iterations=as.integer(fa2_n_iter))
    positions  <-positions[order(as.numeric(names(positions)))]
    fa2.layout <- matrix(unlist(positions), ncol = 2, byrow = TRUE)
    #######################
    rownames(fa2.layout)<-rownames(my.pca)
    
  }
  get_igraph.graph(object)<-G.igraph
  get_fa2_graph.layout(object)<-fa2.layout
  
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Force directed graph analysis is done!.")
}

#'  K-nearest neighbor graph analysis
#' @param  object The SingCellaR object.
#' @param  useIntegrativeEmbeddings is logical, if TRUE the embedding matrix genereated from integrative analysis will be used.
#' @param  integrative_method The name of an integrative method.
#' @param  knn.metric Type of distance metric to use to find nearest neighbors.
#' @param  knn.n_trees The number of trees for building the annoy index. Default 50
#' @param  knn.n_threads The number of threads. Default 1
#' @param  n.dims.use The number of dimensions as the input. Default 30
#' @param  n.neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Default 30
#' @param  niter Integer scalar, the number of iterations to perform.
#' @param  n.seed The number of set seed.
#' @export 
#'


runKNN_Graph <- function(object,dim_reduction_method=c("pca","nnmf"),
                         useIntegrativeEmbeddings=F,integrative_method=c("combat","seurat","harmony","supervised_harmony"),
                         n.dims.use=30,n.neighbors=5,knn.metric=c("euclidean","cosine"),knn.n_trees=50,
                         knn.n_threads=1,niter=500L,n.seed=1){
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  dim_reduction_method<-match.arg(dim_reduction_method)
  integrative_method<-match.arg(integrative_method)
  knn.metric<-match.arg(knn.metric)
  #####################################
  if(useIntegrativeEmbeddings==FALSE){
    if(dim_reduction_method=="pca"){
      res.pca<-get_pca.result(object)
      my.pca<-as.matrix(res.pca$x[,1:n.dims.use])
    }else if(dim_reduction_method=="nnmf"){
      res.pca<-get_nnmf.result(object)
      my.pca<-t(res.pca$H)[,1:n.dims.use]
    }
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="harmony"){
    my.pca<-object@Harmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="supervised_harmony"){
      my.pca<-object@SupervisedHarmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="seurat"){
    my.pca<-object@Seurat.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="combat"){
    my.pca<-object@Combat.embeddings[,1:n.dims.use]
  }
  #####################################
  knn <- identifyKNN(as.matrix(my.pca),n.neighbors=n.neighbors,
                     metric = knn.metric,n_trees = knn.n_trees,
                     n_threads = knn.n_threads)
  #####################################
  links<-compute_weights(knn)
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  G <- graph.data.frame(relations, directed=FALSE)
  G <- igraph::simplify(G)
  ######################################
  print("Processing graph-layout..")
  
  set.seed(n.seed)
  my.layout=layout_with_fr(G,dim=3,niter=niter)
  rownames(my.layout)<-rownames(my.pca)
  #####################################
  get_knn_graph.graph(object)<-G
  get_knn_graph.layout(object)<-my.layout
  #######################################
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("KNN-Graph analysis is done!.")
}

#'  Cell cluster identification 
#' @param  object The SingCellaR object.
#' @param  useIntegrativeEmbeddings is logical, if TRUE the embedding matrix genereated from integrative analysis will be used.
#' @param  integrative_method The name of an integrative method.
#' @param  n.dims.use The number of dimensions as the input. Default 30
#' @param  n.neighbors The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. Default 30
#' @param  knn.n_trees The number of trees for building the annoy index. Default 50
#' @param  knn.metric Type of distance metric to use to find nearest neighbors.
#' @param  knn.n_threads The number of threads. Default 1
#' @export 
#'

identifyClusters <- function(object,dim_reduction_method=c("pca","nnmf","lsi"),useIntegrativeEmbeddings=FALSE,
                             integrative_method=c("combat","seurat","harmony","supervised_harmony"),n.dims.use=30,n.neighbors=30,
                             knn.n_trees=50,knn.metric=c("euclidean","cosine"),knn.n_threads =1){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  dim_reduction_method<-match.arg(dim_reduction_method)
  integrative_method<-match.arg(integrative_method)
  knn.metric<-match.arg(knn.metric)
  #####################################
  if(useIntegrativeEmbeddings==FALSE){
    if(dim_reduction_method=="pca"){
      res.pca<-get_pca.result(object)
      my.pca<-as.matrix(res.pca$x[,1:n.dims.use])
    }else if(dim_reduction_method=="nnmf"){
      res.pca<-get_nnmf.result(object)
      my.pca<-t(res.pca$H)[,1:n.dims.use]
    }else if(dim_reduction_method=="lsi"){
      res.pca<-get_lsi.result(object)
      my.pca<-res.pca$matSVD[,1:n.dims.use]
    }
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="harmony"){
    my.pca<-object@Harmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="supervised_harmony"){
    my.pca<-object@SupervisedHarmony.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="seurat"){
    my.pca<-object@Seurat.embeddings[,1:n.dims.use]
  }else if(useIntegrativeEmbeddings==TRUE & integrative_method=="combat"){
    my.pca<-object@Combat.embeddings[,1:n.dims.use]
  }
  ######################################
  mat<-as.matrix(my.pca)
  ######################################
  knn <- identifyKNN(as.matrix(my.pca),n.neighbors=n.neighbors,
                     metric = knn.metric,n_trees = knn.n_trees,
                     n_threads = knn.n_threads)
  #####################################
  print("Building graph network..")
  links<-compute_weights(knn)
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  G <- graph.data.frame(relations, directed=FALSE)
  G <- igraph::simplify(G)
  print("Building graph network is done.")
  #####################################
  #res.tsne<-get_tsne.result(object)
  ######method:walk_trap###############
  #km_walktrap <- igraph::cluster_walktrap(G)
  #com_walktrap <- km_walktrap$membership
  #com_walktrap.t<-table(com_walktrap)
  #com_walktrap.t<-com_walktrap.t[order(com_walktrap.t,decreasing = T)]
  ########ranking clusters#############
  #com_walktrap.l<-list()
  #for( k in 1:length(com_walktrap.t)){
  #  cl.x<-names(com_walktrap.t)[k]
  #  com_walktrap.l[cl.x]<-k
  #}
  #com_walktrap.ranked<-as.numeric(com_walktrap.l[as.character(com_walktrap)])
  ######################################
  #walktrap.cluster<-data.frame(Cell=rownames(mat),walktrap_cluster=paste("cl",com_walktrap.ranked,sep=""))
  #p<- qplot(0,0, data=walktrap.cluster,colour=walktrap_cluster)+geom_point(size=2)
  #walktrap.cluster$walktrap_cluster_color<-ggplot_build(p)$data[[1]]$colour
  #print("Walktrap analysis is done!.")
  #######method:louvain#################
  print("Process community detection..")
  km_louvain <- igraph::cluster_louvain(G)
  com_louvain <- km_louvain$membership
  com_louvain.t<-table(com_louvain)
  com_louvain.t<-com_louvain.t[order(com_louvain.t,decreasing = T)]
  #########ranking clusters#############
  com_louvain.l<-list()
  for( y in 1:length(com_louvain.t)){
    cl.x<-names(com_louvain.t)[y]
    com_louvain.l[cl.x]<-y
  }
  com_louvain.ranked<-as.numeric(com_louvain.l[as.character(com_louvain)])
  #######################################
  louvain.cluster<-data.frame(Cell=rownames(mat),louvain_cluster=paste("cl",com_louvain.ranked,sep=""))
  p<- qplot(0,0, data=louvain.cluster,colour=louvain_cluster)+geom_point(size=2)
  louvain.cluster$louvain_cluster_color<-ggplot_build(p)$data[[1]]$colour
  print("Louvain analysis is done!.")
  #######method:infomap#################
  #km_infomap <- igraph::cluster_infomap(G,nb.trials = 30)
  #com_infomap <- km_infomap$membership
  #infomap.cluster<-data.frame(Cell=rownames(mat),infomap_cluster=paste("cl",com_infomap,sep=""))
  #infomap.tsne<-merge(res.tsne[,1:3],infomap.cluster,by.x="Cell",by.y="Cell",all=T)
  #p<- qplot(TSNE1,TSNE2, data=infomap.tsne,colour=infomap_cluster)+geom_point(size=2)
  #infomap.tsne$infomap_cluster_color<-ggplot_build(p)$data[[1]]$colour
  #print("Infomap analysis is done!.")
  #########################################
  #walktrap.res<-walktrap.tsne[,c(1,4,5)]
  #louvain.res<-louvain.tsne[,c(1,4,5)]
  #infomap.res<-infomap.tsne[,c(1,4,5)]
  #########combined clusters###############
  #combined.res<-merge(walktrap.cluster,louvain.cluster)
  #combined.res<-merge(combined.res,infomap.res)
  #########################################
  get_clusters(object)<-louvain.cluster
  assign(objName,object,envir=parent.frame())
  invisible(1)
}

#' Identify marker genes for each cluster.
#' @param  object The SingCellaR object.
#' @param  cluster.type The type of clustering method.
#' @param  min.log2FC The minimum cutoff of the log2FC. Default 0.5
#' @param  min.expFraction The minimum expressing cells fraction. Default 0.3
#' @param  limit.top.genes The limited number of differential genes to be shown. Default 500
#' @param  use.sampling.cells is logical. If TRUE, the downsample cells will be processed.
#' @param  use.fraction.of.cells The fraction of downsample. Default 0.2
#' @param  n.seed The number of set seed.
#' @export 
#'

findMarkerGenes <- function(object,cluster.type=c("louvain","walktrap","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                            min.log2FC=0.5,min.expFraction=0.3,limit.top.genes = 500,
                            use.sampling.cells=FALSE,use.fraction.of.cells=0.2,n.seed=1){
  
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
  #######################################
  umi.dat<-get_normalized_umi(object)
  genes.info<-get_genes_metadata(object)
  genes.info<-subset(genes.info,IsExpress==TRUE)
  gene.index<-rownames(umi.dat) %in% rownames(genes.info)
  ######################################
  print(paste("Number of genes for analysis: ",nrow(genes.info),sep=""))
  umi.dat<-umi.dat[gene.index,]
  #########################
  print("Creating binary matrix..")
  binary.mat<-makeBinaryMatrix_for_float(umi.dat)
  
  print("Finding differential genes..")
  MarkerGenes<-list()
  #########################
  for(k in 1:length(clusters.ids)){
    
    cluster.id<-paste("cl",k,sep="")
    start.message<-paste("Identifying marker genes for:",cluster.id)
    print(start.message)  
    
    cols.A<-clusters.info$Cell[clusters.info[,2] ==cluster.id]
    cols.B<-clusters.info$Cell[clusters.info[,2] !=cluster.id]
    
    a.index<-which(colnames(umi.dat) %in% cols.A)
    b.index<-which(colnames(umi.dat) %in% cols.B)
    set.seed(n.seed)
    if(use.sampling.cells==TRUE){
      a.index<-sample(a.index,length(a.index)*use.fraction.of.cells,replace = FALSE)
      b.index<-sample(b.index,length(b.index)*use.fraction.of.cells,replace = FALSE)
    }
    
    cellsA.m<-umi.dat[,a.index]
    cellsB.m<-umi.dat[,b.index]
    
    cellsA.bi<-binary.mat[,a.index]
    #cellsA.bi[cellsA.bi > 0]<-1
    A.freq<-Matrix::rowSums(cellsA.bi)
    
    cellsB.bi<-binary.mat[,b.index]
    #cellsB.bi[cellsB.bi > 0]<-1
    B.freq<-Matrix::rowSums(cellsB.bi)
    
    MeanA<-Matrix::rowMeans(cellsA.m)
    MeanB<-Matrix::rowMeans(cellsB.m)
    Foldchange<-MeanA/MeanB
    
    ExpFractionA=A.freq/ncol(cellsA.m)
    ExpFractionB=B.freq/ncol(cellsB.m)
    
    z0<-data.frame(Gene=rownames(cellsA.m),
                   ExpA=MeanA,
                   ExpB=MeanB,
                   FoldChange=Foldchange,
                   log2FC=log2(Foldchange),
                   ExpFreqA=A.freq,
                   ExpFreqB=B.freq,
                   TotalA=ncol(cellsA.m),
                   TotalB=ncol(cellsB.m),
                   ExpFractionA=ExpFractionA,
                   ExpFractionB=ExpFractionB)
    
    z0<-subset(z0,ExpFractionA >=min.expFraction & ExpFractionA > ExpFractionB)
    z0<-subset(z0,log2FC >=min.log2FC)
    
    if(nrow(z0) > 0){
      cont.t<-data.frame(x1=z0$ExpFreqA,x2=ncol(cellsA.m)-z0$ExpFreqA,y1=z0$ExpFreqB,y2=ncol(cellsB.m)-z0$ExpFreqB)
      print("Processing Fisher's exact test!")
      fishers.p = pbsapply (1:nrow(cont.t),
                            function (x) fisher.test(matrix(c(cont.t[x,1],
                                                              cont.t[x,2],
                                                              cont.t[x,3],
                                                              cont.t[x,4]),
                                                            nrow=2
                            ))$p.value
      )
      
      z0$fishers.pval<-fishers.p
      z0<-subset(z0,fishers.pval < 0.05)
      
      if(nrow(z0) > 1){
        z0<-z0[order(z0$FoldChange,decreasing = T),]
        if(nrow(z0) >limit.top.genes){
          z0<-z0[1:limit.top.genes,]
        }
        cellsA.f.m<-cellsA.m[rownames(cellsA.m) %in% as.character(z0$Gene),]
        cellsB.f.m<-cellsB.m[rownames(cellsB.m) %in% as.character(z0$Gene),]
        cellsA.f.m<-log1p(cellsA.f.m)
        cellsB.f.m<-log1p(cellsB.f.m)
        
        print("Processing Wilcoxon Rank-sum test!")
        wilcoxon.p = pbsapply (1:nrow(cellsA.f.m),
                               function (x) wilcox.test(cellsA.f.m[x,],cellsB.f.m[x,])$p.value
        )
        z0$wilcoxon.pval<-wilcoxon.p
        
        fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
        
        print("Combining p-values using the fisher's method!")
        combined.pvalues = pbsapply (1:nrow(z0),
                                     function (x) fishersMethod(c(z0[x,12],z0[x,13]))
        )
        z0$combined.pval<-combined.pvalues
        z0<-z0[order(z0$wilcoxon.pval),]
        MarkerGenes[[cluster.id]]<-as.data.table(z0)
      }
    }else{
      print(paste(cluster.id,": No gene has passed the criteria!",sep=""))
    }
  }
  if(cluster.type=="walktrap"){
    object@marker.genes$walktrap<-MarkerGenes
  }
  if(cluster.type=="louvain"){
    object@marker.genes$louvain<-MarkerGenes
  }
  if(cluster.type=="kmeans"){
    object@marker.genes$kmeans<-MarkerGenes
  }
  if(cluster.type=="merged_walktrap"){
    object@marker.genes$merged_walktrap<-MarkerGenes
  }
  if(cluster.type=="merged_louvain"){
    object@marker.genes$merged_louvain<-MarkerGenes
  }
  if(cluster.type=="merged_kmeans"){
    object@marker.genes$merged_kmeans<-MarkerGenes
  }
  
  assign(objName,object,envir=parent.frame())
  invisible(1)
}

#' Get marker genes information for any selected cluster.
#' @param  object The SingCellaR object.
#' @param  cluster.type The type of clustering method.
#' @param  cluster_id The specified cluster ID (e.g. cluster1)
#' @export 
#' @return a dataframe of marker information

getMarkerGenesInfo <- function(object,cluster.type=c("louvain","walktrap","kmeans","merged_louvain","merged_walktrap","merged_kmeans"),cluster_id="cl1"){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	cluster.type <- match.arg(cluster.type)
	
	if(cluster.type=="walktrap"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$walktrap
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])	
		return(cluster.genes)
	}else if(cluster.type=="louvain"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$louvain
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else if(cluster.type=="kmeans"){
		markers_genes<-get_marker_genes(object)
		markers_genes<-markers_genes$kmeans
		cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		return(cluster.genes)
	}else if(cluster.type=="merged_louvain"){
		  markers_genes<-get_marker_genes(object)
		  markers_genes<-markers_genes$merged_louvain
		  cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		  return(cluster.genes)
	 }else if(cluster.type=="merged_walktrap"){
		    markers_genes<-get_marker_genes(object)
		    markers_genes<-markers_genes$merged_walktrap
		    cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
		    return(cluster.genes)
	 }else if(cluster.type=="merged_kmeans"){
	   markers_genes<-get_marker_genes(object)
	   markers_genes<-markers_genes$merged_kmeans
	   cluster.genes<-as.data.frame(markers_genes[[cluster_id]])
	   return(cluster.genes)
	}else{
		stop("Need a clustering-method name!")
	}	
}

#' Get cell names from selected clusters.
#' @param  object The SingCellaR object.
#' @param  cluster.type The type of clustering method.
#' @param  cluster_id The specified cluster ID (e.g. cluster1)
#' @export 
#' @return  a vector of cell names.

getCellsFromClusters <- function(object,cluster.type=c("louvain","walktrap","kmeans"),cluster_id=c("cl1","cl2")){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	cluster.type <- match.arg(cluster.type)
	
	if(length(cluster_id)==0){
		stop("Please add cluster ids!")
	}
	
	if(cluster.type=="walktrap"){
		
		cell_names<-c()
		cluster.info<-get_clusters(object)
		for( k in 1:length(cluster_id)){
			cl.x<-cluster_id[k]
			cl.sub<-subset(cluster.info,walktrap_cluster==cl.x)
			cl.name<-as.character(cl.sub$Cell)
			cell_names<-append(cell_names,cl.name)
		}
		return(cell_names)
		
	}else if(cluster.type=="louvain"){
		cell_names<-c()
		cluster.info<-get_clusters(object)
		for( k in 1:length(cluster_id)){
			cl.x<-cluster_id[k]
			cl.sub<-subset(cluster.info,louvain_cluster==cl.x)
			cl.name<-as.character(cl.sub$Cell)
			cell_names<-append(cell_names,cl.name)
		}
		return(cell_names)
		
	}else if(cluster.type=="kmeans"){
		cell_names<-c()
		cluster.info<-get_knn_graph.kmeans.cluster(object)
		for( k in 1:length(cluster_id)){
			cl.x<-cluster_id[k]
			cl.sub<-subset(cluster.info,kmeans_cluster==cl.x)
			cl.name<-as.character(cl.sub$Cell)
			cell_names<-append(cell_names,cl.name)
		}
		return(cell_names)
		
	}else{
		stop("Need a clustering-method name!")
	}	
}

#' Identify differentially expressed genes between a selected cluster vs the rest of clusters
#' @param  object The SingCellaR object.
#' @param  cluster.type The type of clustering method.
#' @param  cluster_id The specified cluster ID (e.g. cluster1)
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.25
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.3
#' @param  write.to.file The output file name.
#' @export 
#' @return a dataframe of differentially expressed genes.
#'

identifyDifferentialGenes_in_a_cluster_vs_the_rest_of_clusters <- function(object,cluster.id="",
                                                                                    method=c("wilcoxon"),
                                                                                    cluster.type=c("louvain","walktrap","kmeans",
                                                                                                   "merged_walktrap","merged_louvain",
                                                                                                   "merged_kmeans"),
                                                                                    min.log2FC=0.25,min.expFraction=0.3,write.to.file=""){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(cluster.id)==0){
    stop("Please input cluster id of interest (e.g. cluster.id='cl1')")
  }
  
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  
  if(cluster.type=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster")]
  }else if(cluster.type=="louvain"){
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
  #######################################
  umi.dat<-get_normalized_umi(object)
  genes.info<-get_genes_metadata(object)
  genes.info<-subset(genes.info,IsExpress==TRUE)
  umi.dat<-umi.dat[rownames(umi.dat) %in% rownames(genes.info),]
  #########################
  start.message<-paste("Identifying differentially expressed genes for:",cluster.id)
  print(start.message)  
  
  cols.A<-clusters.info$Cell[clusters.info[,2] ==cluster.id]
  cols.B<-clusters.info$Cell[clusters.info[,2] !=cluster.id]
  
  cellsA.m<-umi.dat[,colnames(umi.dat) %in% cols.A]
  cellsB.m<-umi.dat[,colnames(umi.dat) %in% cols.B]
  
  if(method=="wilcoxon"){
    
    cellsA.bi<-makeBinaryMatrix_for_float(cellsA.m)
    #cellsA.bi[cellsA.bi > 0]<-1
    A.freq<-Matrix::rowSums(cellsA.bi)
    
    cellsB.bi<-makeBinaryMatrix_for_float(cellsB.m)
    #cellsB.bi[cellsB.bi > 0]<-1
    B.freq<-Matrix::rowSums(cellsB.bi)
    
    
    MeanA<-Matrix::rowMeans(cellsA.m)
    MeanB<-Matrix::rowMeans(cellsB.m)
    Foldchange<-MeanA/MeanB
    
    ExpFractionA=A.freq/ncol(cellsA.m)
    ExpFractionB=B.freq/ncol(cellsB.m)
    
    z0<-data.frame(Gene=rownames(cellsA.m),
                   ExpA=MeanA,
                   ExpB=MeanB,
                   FoldChange=Foldchange,
                   log2FC=log2(Foldchange),
                   ExpFreqA=A.freq,
                   ExpFreqB=B.freq,
                   TotalA=ncol(cellsA.m),
                   TotalB=ncol(cellsB.m),
                   ExpFractionA=ExpFractionA,
                   ExpFractionB=ExpFractionB)
    
    z0<-subset(z0,ExpFractionA >=min.expFraction | ExpFractionB >=min.expFraction)
    z0<-subset(z0,abs(log2FC) >=min.log2FC)
    
    cont.t<-data.frame(x1=z0$ExpFreqA,x2=ncol(cellsA.m)-z0$ExpFreqA,y1=z0$ExpFreqB,y2=ncol(cellsB.m)-z0$ExpFreqB)
    print("Processing Fisher's exact test!")
    fishers.p = pbsapply (1:nrow(cont.t),
                          function (x) fisher.test(matrix(c(cont.t[x,1],
                                                            cont.t[x,2],
                                                            cont.t[x,3],
                                                            cont.t[x,4]),
                                                          nrow=2
                          ))$p.value
    )
    z0$fishers.pval<-fishers.p
    z0<-subset(z0,fishers.pval < 0.1)
    
    cellsA.f.m<-cellsA.m[rownames(cellsA.m) %in% as.character(z0$Gene),]
    cellsB.f.m<-cellsB.m[rownames(cellsB.m) %in% as.character(z0$Gene),]
    cellsA.f.m<-log1p(cellsA.f.m)
    cellsB.f.m<-log1p(cellsB.f.m)
    
    print("Processing Wilcoxon Rank-sum test!")
    wilcoxon.p = pbsapply (1:nrow(cellsA.f.m),
                           function (x) wilcox.test(cellsA.f.m[x,],cellsB.f.m[x,])$p.value
    )
    z0$wilcoxon.pval<-wilcoxon.p
    
    fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
    print("Combining p-values using the fisher's method!")
    combined.pvalues = pbsapply (1:nrow(z0),
                                 function (x) fishersMethod(c(z0[x,12],z0[x,13]))
    )
    z0$combined.pval<-combined.pvalues
    z0$adjusted.pval<-p.adjust(z0$combined.pval,method = "BH")
    z0<-z0[order(z0$adjusted.pval),]
    #####write to a file############
    if(write.to.file!=""){
      write.table(z0,file=write.to.file,
                  append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
    }
    ################################
    return(z0)
  }else{
    stop("Please input the statistical method name!.")
  }
}

#' Identify differentially expressed genes between multiple selected clusters vs the rest of clusters
#' @param  object The SingCellaR object.
#' @param  method  The differential gene expression testing method. Default 'wilcoxon'  test
#' @param  cluster.type The type of clustering method.
#' @param  cluster_id The specified cluster ID (e.g. cluster1)
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.25
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.3
#' @param  write.to.file The output file name.
#' @export 
#' @return a dataframe of differentially expressed genes.
#'

identifyDifferentialGenes_in_multiple_clusters_vs_the_rest_of_clusters <- function(object,cluster.ids=c(),
                                                                                            method=c("wilcoxon"),
                                                                                            cluster.type=c("louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                                                                            min.log2FC=0.25,min.expFraction=0.3,write.to.file=""){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(cluster.ids)==0){
    stop("Please input cluster id of interest (e.g. cluster.ids=c('cl1','cl2')")
  }
  
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  
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
  #######################################
  umi.dat<-get_normalized_umi(object)
  genes.info<-get_genes_metadata(object)
  genes.info<-subset(genes.info,IsExpress==TRUE)
  umi.dat<-umi.dat[rownames(umi.dat) %in% rownames(genes.info),]
  #########################
  start.message<-paste("Identifying differentially expressed genes for:",toString(cluster.ids))
  print(start.message)  
  
  cols.A<-clusters.info$Cell[clusters.info[,2] %in% cluster.ids==T]
  cols.B<-clusters.info$Cell[clusters.info[,2] %in% cluster.ids==F]
  
  cellsA.m<-umi.dat[,colnames(umi.dat) %in% cols.A]
  cellsB.m<-umi.dat[,colnames(umi.dat) %in% cols.B]
  
  if(method=="wilcoxon"){
    
    cellsA.bi<-makeBinaryMatrix_for_float(cellsA.m)
    #cellsA.bi[cellsA.bi > 0]<-1
    A.freq<-Matrix::rowSums(cellsA.bi)
    
    cellsB.bi<-makeBinaryMatrix_for_float(cellsB.m)
    #cellsB.bi[cellsB.bi > 0]<-1
    B.freq<-Matrix::rowSums(cellsB.bi)
    
    
    MeanA<-Matrix::rowMeans(cellsA.m)
    MeanB<-Matrix::rowMeans(cellsB.m)
    Foldchange<-MeanA/MeanB
    
    ExpFractionA=A.freq/ncol(cellsA.m)
    ExpFractionB=B.freq/ncol(cellsB.m)
    
    z0<-data.frame(Gene=rownames(cellsA.m),
                   ExpA=MeanA,
                   ExpB=MeanB,
                   FoldChange=Foldchange,
                   log2FC=log2(Foldchange),
                   ExpFreqA=A.freq,
                   ExpFreqB=B.freq,
                   TotalA=ncol(cellsA.m),
                   TotalB=ncol(cellsB.m),
                   ExpFractionA=ExpFractionA,
                   ExpFractionB=ExpFractionB)
    
    z0<-subset(z0,ExpFractionA >=min.expFraction | ExpFractionB >=min.expFraction)
    z0<-subset(z0,abs(log2FC) >=min.log2FC)
    
    cont.t<-data.frame(x1=z0$ExpFreqA,x2=ncol(cellsA.m)-z0$ExpFreqA,y1=z0$ExpFreqB,y2=ncol(cellsB.m)-z0$ExpFreqB)
    print("Processing Fisher's exact test!")
    fishers.p = pbsapply (1:nrow(cont.t),
                          function (x) fisher.test(matrix(c(cont.t[x,1],
                                                            cont.t[x,2],
                                                            cont.t[x,3],
                                                            cont.t[x,4]),
                                                          nrow=2
                          ))$p.value
    )
    z0$fishers.pval<-fishers.p
    z0<-subset(z0,fishers.pval < 0.1)
    
    cellsA.f.m<-cellsA.m[rownames(cellsA.m) %in% as.character(z0$Gene),]
    cellsB.f.m<-cellsB.m[rownames(cellsB.m) %in% as.character(z0$Gene),]
    cellsA.f.m<-log1p(cellsA.f.m)
    cellsB.f.m<-log1p(cellsB.f.m)
    
    print("Processing Wilcoxon Rank-sum test!")
    wilcoxon.p = pbsapply (1:nrow(cellsA.f.m),
                           function (x) wilcox.test(cellsA.f.m[x,],cellsB.f.m[x,])$p.value
    )
    z0$wilcoxon.pval<-wilcoxon.p
    
    fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
    print("Combining p-values using the fisher's method!")
    combined.pvalues = pbsapply (1:nrow(z0),
                                 function (x) fishersMethod(c(z0[x,12],z0[x,13]))
    )
    z0$combined.pval<-combined.pvalues
    z0$adjusted.pval<-p.adjust(z0$combined.pval,method = "BH")
    z0<-z0[order(z0$adjusted.pval),]
    #####write to a file############
    if(write.to.file!=""){
      write.table(z0,file=write.to.file,
                  append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
    }
    ################################
    return(z0)
  }else{
    stop("Please input the statistical method name!.")
  }
}

#' Identify differentially expressed genes for all identified clusters
#' @param  object The SingCellaR object.
#' @param  method  The differential gene expression testing method. Default 'wilcoxon'  test
#' @param  cluster.type The type of clustering method.
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.20
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.2
#' @export 


identifyDifferentialGenes_for_all_clusters <- function(object,method=c("wilcoxon"),
                                                       cluster.type=c("louvain","walktrap","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                                       min.log2FC=0.20,min.expFraction=0.20){
  
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
  }else if(cluster.type=="merged_walktrap"){
    clusters.info<-clusters.info[,c("Cell","merged_walktrap")]
  }else if(cluster.type=="merged_louvain"){
    clusters.info<-clusters.info[,c("Cell","merged_louvain")]
  }else if(cluster.type=="merged_kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans")]
  }else{
    stop("Need a clustering-method name!")
  }
  #######################################
  clusters.ids<-1:length(unique(clusters.info[,2]))
  #######################################
  DiffGenes<-list()
  for(k in 1:length(clusters.ids)){
    cluster.id<-paste("cl",k,sep="")
    my.diff.genes<-identifyDifferentialGenes_in_a_selected_cluster_vs_the_rest_of_clusters(object,cluster.id = cluster.id,method = method,cluster.type = cluster.type,min.log2FC = min.log2FC,min.expFraction = min.expFraction)
    DiffGenes[[cluster.id]]<-as.data.table(my.diff.genes)
  }
  return(DiffGenes)
}

#' Identify GSEA's preranked genes for a selected cluster vs the rest of clusters
#' @param  object The SingCellaR object.
#' @param  cluster.id The specified cluster id (e.g. 'cluster1')
#' @param  cluster.type The type of clustering method.
#' @param  fishers_exact_test The fisher's exact test cutoff p-value. Default 0.1
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.01
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.1
#' @param  write.to.file The output file name.
#' @export 
#' @return a dataframe of preranked genes.

identifyGSEAPrerankedGenes_in_a_cluster_vs_the_rest_of_clusters <- function(object,cluster.id="",
                                                                                     cluster.type=c("louvain","kmeans","merged_walktrap",
                                                                                                    "merged_louvain","merged_kmeans"),
                                                                                     fishers_exact_test=0.1,min.expFraction=0.01,
                                                                                     min.log2FC=0.1,write.to.file=""){
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(length(cluster.id)==0){
		stop("Please input cluster id of interest (e.g. cluster.id='cl1')")
	}
	
	clusters.info<-get_clusters(object)
	cluster.type <- match.arg(cluster.type)
	
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
	#######################################
	umi.dat<-get_normalized_umi(object)
	genes.info<-get_genes_metadata(object)
	genes.info<-subset(genes.info,IsExpress==TRUE)
	umi.dat<-umi.dat[rownames(umi.dat) %in% rownames(genes.info),]
	#######################################
	start.message<-paste("Identifying differentially expressed genes for:",cluster.id)
	print(start.message)  
	
	cols.A<-clusters.info$Cell[clusters.info[,2] ==cluster.id]
	cols.B<-clusters.info$Cell[clusters.info[,2] !=cluster.id]
	
	cellsA.m<-umi.dat[,colnames(umi.dat) %in% cols.A]
	cellsB.m<-umi.dat[,colnames(umi.dat) %in% cols.B]
	
		
	  cellsA.bi<-makeBinaryMatrix_for_float(cellsA.m)
		#cellsA.bi[cellsA.bi > 0]<-1
	  A.freq<-Matrix::rowSums(cellsA.bi)
		
	  cellsB.bi<-makeBinaryMatrix_for_float(cellsB.m)
		#cellsB.bi[cellsB.bi > 0]<-1
	  B.freq<-Matrix::rowSums(cellsB.bi)
		
		
		MeanA<-Matrix::rowMeans(cellsA.m)
		MeanB<-Matrix::rowMeans(cellsB.m)
		Foldchange<-MeanA/MeanB
		
		ExpFractionA=A.freq/ncol(cellsA.m)
		ExpFractionB=B.freq/ncol(cellsB.m)
		
		z0<-data.frame(Gene=rownames(cellsA.m),
				ExpA=MeanA,
				ExpB=MeanB,
				FoldChange=Foldchange,
				log2FC=log2(Foldchange),
				ExpFreqA=A.freq,
				ExpFreqB=B.freq,
				TotalA=ncol(cellsA.m),
				TotalB=ncol(cellsB.m),
				ExpFractionA=ExpFractionA,
				ExpFractionB=ExpFractionB)
		
		z0<-subset(z0,(ExpFractionA >=min.expFraction | ExpFractionB >=min.expFraction))
		z0<-subset(z0,abs(log2FC) >= min.log2FC)
		
		cont.t<-data.frame(x1=z0$ExpFreqA,x2=ncol(cellsA.m)-z0$ExpFreqA,y1=z0$ExpFreqB,y2=ncol(cellsB.m)-z0$ExpFreqB)
		print("Processing Fisher's exact test!")
		fishers.p = pbsapply (1:nrow(cont.t),
				function (x) fisher.test(matrix(c(cont.t[x,1],
											cont.t[x,2],
											cont.t[x,3],
											cont.t[x,4]),
									nrow=2
							))$p.value
		)
		z0$fishers.pval<-fishers.p
		z0<-subset(z0,fishers.pval < fishers_exact_test)
		z0$adjusted.pval<-p.adjust(z0$fishers.pval,method = "BH")
		z0<-z0[order(z0$log2FC,decreasing = T),]
		#####write to a file############
		if(write.to.file!=""){
			write.table(z0,file=write.to.file,
					append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
		}
		################################
		return(z0)
}

#' Identify GSEA's preranked genes for all identified clusters
#' @param  object The SingCellaR object.
#' @param  cluster.type The type of clustering method.
#' @param  fishers_exact_test The fisher's exact test cutoff p-value. Default 0.1
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.01
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.1
#' @export 

identifyGSEAPrerankedGenes_for_all_clusters <- function(object,cluster.type=c("louvain","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                                       fishers_exact_test=0.1,min.expFraction=0.01,min.log2FC=0.1){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  
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
  #######################################
  clusters.ids<-1:length(unique(clusters.info[,2]))
  #######################################
  DiffGenes<-list()
  for(k in 1:length(clusters.ids)){
    cluster.id<-paste("cl",k,sep="")
    my.diff.genes<-identifyGSEAPrerankedGenes_in_a_cluster_vs_the_rest_of_clusters(object,cluster.id = cluster.id,
                                                                                            cluster.type = cluster.type,
                                                                                            fishers_exact_test=fishers_exact_test,
                                                                                            min.expFraction=min.expFraction,
                                                                                            min.log2FC=min.expFraction)
    DiffGenes[[cluster.id]]<-as.data.table(my.diff.genes)
  }
  return(DiffGenes)
}

#' Differential gene expression analysis
#' @param  method  The differential gene expression testing method. Default 'wilcoxon'  test
#' @param  objectA The SingCellaR object A
#' @param  objectB The SingCellaR object B
#' @param  cellsA  The vector of cell names from object A
#' @param  cellsB  The vector of cell names from object B
#' @param  write.to.file The output file name.
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.3
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.25
#' @export 

identifyDifferentialGenes <- function(method=c("wilcoxon"),objectA=objectA,objectB=objectB,
                                      cellsA=c(),cellsB=c(),write.to.file="",min.log2FC=0.25,min.expFraction=0.3){
	
	if(!is(objectA,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	if(!is(objectB,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	
	if(method=="wilcoxon"){
		genesA.info<-get_genes_metadata(objectA)
		genesA.info<-subset(genesA.info,IsExpress==TRUE)
		
		genesB.info<-get_genes_metadata(objectB)
		genesB.info<-subset(genesB.info,IsExpress==TRUE)
		
		expr.genes<-unique(c(rownames(genesA.info),rownames(genesB.info)))
		########################
		umiA.dat<-get_normalized_umi(objectA)
		umiB.dat<-get_normalized_umi(objectB)
		
		umiA.dat<-umiA.dat[rownames(umiA.dat) %in% expr.genes,]
		umiB.dat<-umiB.dat[rownames(umiB.dat) %in% expr.genes,]
		#########################
		cellsA.m<-umiA.dat[,colnames(umiA.dat) %in% cellsA]
		cellsB.m<-umiB.dat[,colnames(umiB.dat) %in% cellsB]
		
		cellsA.bi<-makeBinaryMatrix_for_float(cellsA.m)
		#cellsA.bi[cellsA.bi > 0]<-1
		A.freq<-Matrix::rowSums(cellsA.bi)
		
		cellsB.bi<-makeBinaryMatrix_for_float(cellsB.m)
		#cellsB.bi[cellsB.bi > 0]<-1
		B.freq<-Matrix::rowSums(cellsB.bi)
		
		
		MeanA<-Matrix::rowMeans(cellsA.m)
		MeanB<-Matrix::rowMeans(cellsB.m)
		Foldchange<-MeanA/MeanB
		
		ExpFractionA=A.freq/ncol(cellsA.m)
		ExpFractionB=B.freq/ncol(cellsB.m)
		
		z0<-data.frame(Gene=rownames(cellsA.m),
				ExpA=MeanA,
				ExpB=MeanB,
				FoldChange=Foldchange,
				log2FC=log2(Foldchange),
				ExpFreqA=A.freq,
				ExpFreqB=B.freq,
				TotalA=ncol(cellsA.m),
				TotalB=ncol(cellsB.m),
				ExpFractionA=ExpFractionA,
				ExpFractionB=ExpFractionB)
		
		z0<-subset(z0,ExpFractionA >=min.expFraction | ExpFractionB >=min.expFraction)
		z0<-subset(z0,abs(log2FC) >=min.log2FC)
		
		cont.t<-data.frame(x1=z0$ExpFreqA,x2=ncol(cellsA.m)-z0$ExpFreqA,y1=z0$ExpFreqB,y2=ncol(cellsB.m)-z0$ExpFreqB)
		print("Processing Fisher's exact test!")
		fishers.p = pbsapply (1:nrow(cont.t),
				function (x) fisher.test(matrix(c(cont.t[x,1],
											cont.t[x,2],
											cont.t[x,3],
											cont.t[x,4]),
									nrow=2
							))$p.value
		)
		z0$fishers.pval<-fishers.p
		z0<-subset(z0,fishers.pval < 0.1)
		
		cellsA.f.m<-cellsA.m[rownames(cellsA.m) %in% as.character(z0$Gene),]
		cellsB.f.m<-cellsB.m[rownames(cellsB.m) %in% as.character(z0$Gene),]
		cellsA.f.m<-log1p(cellsA.f.m)
		cellsB.f.m<-log1p(cellsB.f.m)
		
		print("Processing Wilcoxon Rank-sum test!")
		wilcoxon.p = pbsapply (1:nrow(cellsA.f.m),
				function (x) wilcox.test(cellsA.f.m[x,],cellsB.f.m[x,])$p.value
		)
		z0$wilcoxon.pval<-wilcoxon.p
		
		fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
		print("Combining p-values using the fisher's method!")
		combined.pvalues = pbsapply (1:nrow(z0),
				function (x) fishersMethod(c(z0[x,12],z0[x,13]))
		)
		z0$combined.pval<-combined.pvalues
		z0$adjusted.pval<-p.adjust(z0$combined.pval,method = "BH")
		z0<-z0[order(z0$adjusted.pval),]
		#####write to a file############
		if(write.to.file!=""){
			write.table(z0,file=write.to.file,
					append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
		}
		################################
		return(z0)
	}else{
		stop("Please input the statistical method name!.")
	}
}

#' Identify GSEA preranked genes 
#' @param  method  The differential gene expression testing method. Default 'wilcoxon'  test
#' @param  objectA The SingCellaR object A
#' @param  objectB The SingCellaR object B
#' @param  cellsA  The vector of cell names from object A
#' @param  cellsB  The vector of cell names from object B
#' @param  write.to.file The output file name.
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.01
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.1
#' @export 

identifyGSEAPrerankedGenes <- function(method=c("wilcoxon"),objectA=objectA,
                                       objectB=objectB,cellsA=c(),cellsB=c(),
                                       write.to.file="",min.log2FC=0.1,min.expFraction=0.01){
  
  if(!is(objectA,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(!is(objectB,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  if(method=="wilcoxon"){
    genesA.info<-get_genes_metadata(objectA)
    genesA.info<-subset(genesA.info,IsExpress==TRUE)
    
    genesB.info<-get_genes_metadata(objectB)
    genesB.info<-subset(genesB.info,IsExpress==TRUE)
    
    expr.genes<-unique(c(rownames(genesA.info),rownames(genesB.info)))
    ########################
    umiA.dat<-get_normalized_umi(objectA)
    umiB.dat<-get_normalized_umi(objectB)
    
    umiA.dat<-umiA.dat[rownames(umiA.dat) %in% expr.genes,]
    umiB.dat<-umiB.dat[rownames(umiB.dat) %in% expr.genes,]
    #########################
    cellsA.m<-umiA.dat[,colnames(umiA.dat) %in% cellsA]
    cellsB.m<-umiB.dat[,colnames(umiB.dat) %in% cellsB]
    
    cellsA.bi<-makeBinaryMatrix_for_float(cellsA.m)
    #cellsA.bi[cellsA.bi > 0]<-1
    A.freq<-Matrix::rowSums(cellsA.bi)
    
    cellsB.bi<-makeBinaryMatrix_for_float(cellsB.m)
    #cellsB.bi[cellsB.bi > 0]<-1
    B.freq<-Matrix::rowSums(cellsB.bi)
    
    
    MeanA<-Matrix::rowMeans(cellsA.m)
    MeanB<-Matrix::rowMeans(cellsB.m)
    Foldchange<-MeanA/MeanB
    
    ExpFractionA=A.freq/ncol(cellsA.m)
    ExpFractionB=B.freq/ncol(cellsB.m)
    
    z0<-data.frame(Gene=rownames(cellsA.m),
                   ExpA=MeanA,
                   ExpB=MeanB,
                   FoldChange=Foldchange,
                   log2FC=log2(Foldchange),
                   ExpFreqA=A.freq,
                   ExpFreqB=B.freq,
                   TotalA=ncol(cellsA.m),
                   TotalB=ncol(cellsB.m),
                   ExpFractionA=ExpFractionA,
                   ExpFractionB=ExpFractionB)
    
    z0<-subset(z0,ExpFractionA >=min.expFraction | ExpFractionB >=min.expFraction)
    z0<-subset(z0,abs(log2FC) >=min.log2FC)
    
    cont.t<-data.frame(x1=z0$ExpFreqA,x2=ncol(cellsA.m)-z0$ExpFreqA,y1=z0$ExpFreqB,y2=ncol(cellsB.m)-z0$ExpFreqB)
    print("Processing Fisher's exact test!")
    fishers.p = pbsapply (1:nrow(cont.t),
                          function (x) fisher.test(matrix(c(cont.t[x,1],
                                                            cont.t[x,2],
                                                            cont.t[x,3],
                                                            cont.t[x,4]),
                                                          nrow=2
                          ))$p.value
    )
    z0$fishers.pval<-fishers.p
    z0<-subset(z0,fishers.pval < 0.1)
    
    cellsA.f.m<-cellsA.m[rownames(cellsA.m) %in% as.character(z0$Gene),]
    cellsB.f.m<-cellsB.m[rownames(cellsB.m) %in% as.character(z0$Gene),]
    cellsA.f.m<-log1p(cellsA.f.m)
    cellsB.f.m<-log1p(cellsB.f.m)
    
    print("Processing Wilcoxon Rank-sum test!")
    wilcoxon.p = pbsapply (1:nrow(cellsA.f.m),
                           function (x) wilcox.test(cellsA.f.m[x,],cellsB.f.m[x,])$p.value
    )
    z0$wilcoxon.pval<-wilcoxon.p
    
    fishersMethod <-function(x) pchisq(-2 * sum(log(x)),df=2*length(x),lower.tail=FALSE)
    print("Combining p-values using the fisher's method!")
    combined.pvalues = pbsapply (1:nrow(z0),
                                 function (x) fishersMethod(c(z0[x,12],z0[x,13]))
    )
    z0$combined.pval<-combined.pvalues
    z0$adjusted.pval<-p.adjust(z0$combined.pval,method = "BH")
    z0<-z0[order(z0$adjusted.pval),]
    #####write to a file############
    if(write.to.file!=""){
      write.table(z0,file=write.to.file,
                  append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
    }
    ################################
    return(z0)
  }else{
    stop("Please input the statistical method name!.")
  }
}

#' Run GSEA analysis 
#' @param  objectA The SingCellaR object A
#' @param  objectB The SingCellaR object B
#' @param  cellsA  The vector of cell names from object A
#' @param  cellsB  The vector of cell names from object B
#' @param  GSEAPrerankedGenes_info The dataframe object of preranked genes.
#' @param  gmt.file The GMT file name
#' @param  plotGSEA  is logical. If TRUE, this function will plot GSEA enrichment.
#' @param  nTopGenesets The number of showing top gene sets.
#' @param  minSize The cutoff minimum number of genes in each gene set. Gene set that contains the number of genes lower than this number will be excluded. Default 5
#' @param  maxSize The cutoff maxiumum number of genes in each gene set. Gene set that contains the number of genes higher than this number will be excluded. Default 1000
#' @param  eps The eps paramenter for fgsea, this parameter sets the boundary for calculating the p value. Default 1e-10
#' @export 

Run_fGSEA_analysis <- function(objectA=objectA,objectB=objectB,cellsA=c(),cellsB=c(),
                               GSEAPrerankedGenes_info="",gmt.file="",plotGSEA=T,
                               nTopGenesets=10,minSize=5,maxSize=1000,eps = 1e-10){
  
  
  if(!is(objectA,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(!is(objectB,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(gmt.file==""){
    stop("Need the GMT file!")
  }
  ###################################
  diff_cl<-GSEAPrerankedGenes_info
  ###################################
  i.index<-which(diff_cl$adjusted.pval==0)
  if(length(i.index) > 0){
    diff_cl$adjusted.pval[i.index]<-min(diff_cl$adjusted.pval[-c(i.index)])
  }
  diff_cl$prerank.genes<- (-log10(diff_cl$adjusted.pval))*diff_cl$log2FC
  #########remove -inf genes#########
  inf.index<-which(is.infinite(diff_cl$prerank.genes)==T)
  if(length(inf.index) > 0){
    diff_cl<-diff_cl[-c(inf.index),]
  } 
  ###################################
  diff_cl<-diff_cl[order(diff_cl$prerank.genes,decreasing = T),]
  prerank.genes<- diff_cl$prerank.genes
  names(prerank.genes)<-diff_cl$Gene
  ###################################
  my.db<- gmtPathways(gmt.file)
  
  print("Processing GSEA!")
  fgseaRes <- fgseaMultilevel(pathways = my.db, 
                              stats = prerank.genes,
                              minSize=minSize,
                              maxSize=maxSize,
                              eps = eps)
  
  fgseaRes<-fgseaRes[order(fgseaRes$padj), ]
  ################################
  if(plotGSEA==T){
    topUp <- fgseaRes[ES > 0][head(order(pval), n=nTopGenesets), pathway]
    topDown <- fgseaRes[ES < 0][head(order(pval), n=nTopGenesets), pathway]
    topPathways <- c(topUp, rev(topDown))
    plotGseaTable(my.db[topPathways], prerank.genes, fgseaRes,gseaParam=0.5)
  }
  fgseaRes.return<-list()
  fgseaRes.return[["fgseaResult"]]<-fgseaRes
  fgseaRes.return[["preRankedGenes"]]<-prerank.genes
  fgseaRes.return[["pathways"]]<-my.db
  return(fgseaRes.return)
}
#' Run GSEA analysis for a selected cluster vs the rest of clusters 
#' @param  object  The SingCellaR object.
#' @param  cluster.id The specified cluster id (e.g. 'cluster1').
#' @param  cluster.type The type of clustering method.
#' @param  diff.gene.method The method for differential gene expression analysis. Default 'wilcoxon'.
#' @param  gmt.file The GMT file name.
#' @param  fishers_exact_test The fisher's exact test cutoff p-value.
#' @param  plotGSEA  is logical. If TRUE, this function will plot GSEA enrichment.
#' @param  nTopGenesets The number of showing top gene sets.
#' @param  minSize The cutoff minimum number of genes in each gene set. Gene set that contains the number of genes lower than this number will be excluded. Default 5
#' @param  maxSize The cutoff maxiumum number of genes in each gene set. Gene set that contains the number of genes higher than this number will be excluded. Default 2500
#' @param  eps The eps paramenter for fgsea, this parameter sets the boundary for calculating the p value. Default 1e-10
#' @export 

Run_fGSEA_for_a_selected_cluster_vs_the_rest_of_clusters <- function(object,cluster.id="",cluster.type=c("louvain","walktrap","kmeans","merged_walktrap","merged_louvain","merged_kmeans"),
                                                                     diff.gene.method="wilcoxon",gmt.file="",fishers_exact_test=0.1,plotGSEA=T,
                                                                     nTopGenesets=10,minSize=5,maxSize=2500,eps = 1e-10){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(gmt.file==""){
    stop("Need the GMT file!")
  }
  ########################
  cluster.type <- match.arg(cluster.type)
  diff_cl<-identifyGSEAPrerankedGenes_in_a_selected_cluster_vs_the_rest_of_clusters(object,cluster.id = cluster.id,
                                                                                   method=diff.gene.method,
                                                                                   cluster.type = cluster.type)
 
  i.index<-which(diff_cl$adjusted.pval==0)
  if(length(i.index) > 0){
    diff_cl$adjusted.pval[i.index]<-min(diff_cl$adjusted.pval[-c(i.index)])
  }
  diff_cl$prerank.genes<- (-log10(diff_cl$adjusted.pval))*diff_cl$log2FC
  #########remove -inf genes#########
  inf.index<-which(is.infinite(diff_cl$prerank.genes)==T)
  if(length(inf.index) > 0){
    diff_cl<-diff_cl[-c(inf.index),]
  } 
  ###################################
  prerank.genes<- diff_cl$prerank.genes
  names(prerank.genes)<-diff_cl$Gene
  ###################################
  my.db<- gmtPathways(gmt.file)
  
  print("Processing GSEA!")
  fgseaRes <- fgseaMultilevel(pathways = my.db, 
                    stats = prerank.genes,
                    minSize=minSize,
                    maxSize=maxSize,
                    eps = eps)
  
  fgseaRes<-fgseaRes[order(fgseaRes$padj), ]
  ################################
  if(plotGSEA==T){
    topUp <- fgseaRes[ES > 0][head(order(pval), n=nTopGenesets), pathway]
    topDown <- fgseaRes[ES < 0][head(order(pval), n=nTopGenesets), pathway]
    topPathways <- c(topUp, rev(topDown))
    plotGseaTable(my.db[topPathways], prerank.genes, fgseaRes,gseaParam=0.5)
  }
  fgseaRes.return<-list()
  fgseaRes.return[["fgseaResult"]]<-fgseaRes
  fgseaRes.return[["preRankedGenes"]]<-prerank.genes
  fgseaRes.return[["pathways"]]<-my.db
  return(fgseaRes.return)
}


#' Run GSEA analysis for multiple comparisons
#' @param  object  The SingCellaR object.
#' @param  GSEAPrerankedGenes_list The list of preranked genes.
#' @param  gmt.file The GMT file name.
#' @param  minSize The cutoff minimum number of genes in each gene set. Gene set that contains the number of genes lower than this number will be excluded. Default 5
#' @param  maxSize The cutoff maxiumum number of genes in each gene set. Gene set that contains the number of genes higher than this number will be excluded. Default 2500
#' @param  eps The eps paramenter for fgsea, this parameter sets the boundary for calculating the p value. Default 1e-10
#' @param  n.seed The set seed number.
#' @export 
#' @return a dataframe of GSEA analysis results.


Run_fGSEA_for_multiple_comparisons <- function(GSEAPrerankedGenes_list,gmt.file="",minSize=5,maxSize=2500,eps = 1e-10, n.seed=1){
  
  if(gmt.file==""){
    stop("Need the GMT file!")
  }
  ########################
  cluster.ids<-names(GSEAPrerankedGenes_list)
  ########################
  gsea.results<-data.table()
  ########################
  set.seed(n.seed)
  for(k in 1:length(cluster.ids)){
    
    my.cluster.id<-cluster.ids[k]
    diff_cl<-GSEAPrerankedGenes_list[[my.cluster.id]]
    
    i.index<-which(diff_cl$adjusted.pval==0)
    
    if(length(i.index) > 0){
      diff_cl$adjusted.pval[i.index]<-min(diff_cl$adjusted.pval[-c(i.index)])
    }
    diff_cl$prerank.genes<- (-log10(diff_cl$adjusted.pval))*diff_cl$log2FC
    #########remove -inf genes#########
    inf.index<-which(is.infinite(diff_cl$prerank.genes)==T)
    if(length(inf.index) > 0){
      diff_cl<-diff_cl[-c(inf.index)]
    } 
    ###################################
    prerank.genes<- diff_cl$prerank.genes
    names(prerank.genes)<-diff_cl$Gene
    ###################################
    my.db<- gmtPathways(gmt.file)
    ###################################
    print(paste("Processing fGSEA for:",my.cluster.id,sep=""))
    fgseaRes <- fgseaMultilevel(pathways = my.db, 
                      stats = prerank.genes,
                      minSize=minSize,
                      maxSize=maxSize,
                      eps = eps)
    fgseaRes$EnrichedIn[fgseaRes$ES > 0]<-my.cluster.id
    fgseaRes$EnrichedIn[fgseaRes$ES < 0]<-"the_rest_of_clusters"
    fgseaRes$cluster<-my.cluster.id
    fgseaRes.f<-fgseaRes[,c("pathway","pval","padj","ES","NES","EnrichedIn","cluster")]
    gsea.results<-rbind(gsea.results,fgseaRes.f)
  }
  return(gsea.results)
}

#DensityClusteringOnTSNE <- function(object,n.max.cluster=10){
#	objName <- deparse(substitute(object))
#	if(!is(object,"SingCellaR")){
#		stop("Need to initialize the SingCellaR object")
#	}
#	##################################
#	res.tsne<-get_tsne.result(object)
#	TSNE.mat<-res.tsne[c("TSNE1","TSNE2")]
#	rownames(TSNE.mat)<-res.tsne$Cell
#	##################################
#	dataClust <- densityClust::densityClust(TSNE.mat, gaussian = gaussian)
#	delta_rho_df <- data.frame("delta" = dataClust$delta, "rho" = dataClust$rho)
#	rho_threshold <- 0
#	num_clusters<-n.max.cluster
#	delta_threshold <- sort(delta_rho_df$delta, decreasing = T)[num_clusters]-.Machine$double.eps
#
#	dataClust <- densityClust::findClusters(dataClust,rho = rho_threshold, delta = delta_threshold)
#	density.cluster<-data.frame(Cell=rownames(TSNE.mat),density_cluster=paste("cl",dataClust$clusters,sep="_"))
#	#######################################
#	#########UPDATED TSNE WITH CLUSTERS#####
#	res.updated.tsne<-merge(res.tsne[,1:8],density.cluster,by.x="Cell",by.y="Cell",all=T)
#
#	p<- qplot(TSNE1,TSNE2, data=res.updated.tsne,colour=density_cluster)+geom_point(size=2)
#	res.updated.tsne$density_cluster_color<-ggplot_build(p)$data[[1]]$colour
#	#######################################
#	get_tsne.result(object)<-res.updated.tsne
#	assign(objName,object,envir=parent.frame())
#	invisible(1)
#	print("Clustering analysis is done!.")
#}

#' Run K-Means clustering on the KNN-graph
#' @param  object  The SingCellaR object.
#' @param  k The number of k-means clusters.
#' @param  iter.max The maximum number of k-mean's iterations.
#' @export 

clustering_KMeansOnKNN_Graph <- function(object,k=10,iter.max=100){
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR")){
		stop("Need to initialize the SingCellaR object")
	}
	#######################################
	res.tsne<-get_tsne.result(object)
	#######################################
	k.mat<-get_knn_graph.layout(object)
	my.layout<-as.data.frame(k.mat)
	my.layout$Cell<-rownames(my.layout)
	#######################################
	#######################################
	set.seed(1)
	k.cluster<-kmeans(k.mat,k,iter.max = iter.max)
	
	com.k <- k.cluster$cluster
	com.t<-table(com.k)
	com.t<-com.t[order(com.t,decreasing = T)]
	########ranking clusters#############
	com.l<-list()
	for( k in 1:length(com.t)){
		cl.x<-names(com.t)[k]
		com.l[cl.x]<-k
	}
	com.ranked<-as.numeric(com.l[as.character(com.k)])
	
	k.cluster<-data.frame(Cell=rownames(my.layout),kmeans_cluster=paste("cl",com.ranked,sep=""))
	k.tsne<-merge(res.tsne[,1:3], k.cluster,by.x="Cell",by.y="Cell",all=T)
	p<- qplot(TSNE1,TSNE2, data=k.tsne,colour=kmeans_cluster)+geom_point(size=2)
	k.tsne$kmeans_cluster_color<-ggplot_build(p)$data[[1]]$colour
	#####################################
	kmeans.res<-k.tsne[,c(1,4,5)]
	get_knn_graph.kmeans.cluster(object)<-kmeans.res
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("KMeans analysis on the KNN-graph is done!.")
}

#' Identify the clusters that potenntially can be merged.
#' @param  object  The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  min.expFraction The minimum fraction of expressing cell frequency cutoff. Default 0.3
#' @param  min.log2FC The minimum log2FC cutoff. Default 0.3
#' @param  min.jaccard_similarity The minimum jaccard similarity.
#' @export 

identify_potentially_merging_clusters <- function(object,cluster.type=c("louvain","walktrap","kmeans"),
                                                 min.log2FC=0.3,min.expFraction=0.3,min.jaccard_similarity=0.50){
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
  jaccard.score<-c()
  for(k in 1:nrow(clusters.combn)){
    
    cl.x1<-clusters.combn[k,1]
    cl.x2<-clusters.combn[k,2]
    id.x1<-clusters.combn.ids[k,1]
    id.x2<-clusters.combn.ids[k,2]
    markers.x1<-getMarkerGenesInfo(object,cluster.type=cluster.type,cluster_id=cl.x1)
    markers.x1<-subset(markers.x1,fishers.pval< 0.05 
                       & wilcoxon.pval < 0.05 
                       & ExpFractionA > ExpFractionB
                       & ExpFractionA > min.expFraction 
                       & log2FC >= min.log2FC)
    
    markers.x2<-getMarkerGenesInfo(object,cluster.type=cluster.type,cluster_id=cl.x2)
    markers.x2<-subset(markers.x2,fishers.pval< 0.05 
                       & wilcoxon.pval < 0.05 
                       & ExpFractionA > ExpFractionB
                       & ExpFractionA > min.expFraction 
                       & log2FC >= min.log2FC)
    
    a<-intersect(markers.x1$Gene,markers.x2$Gene)
    b<-union(markers.x1$Gene,markers.x2$Gene)
    my.jaccard.index<-(length(a)/length(b))
    jaccard.score<-append(jaccard.score,my.jaccard.index)
  }
  
  clusters.combn.f<-as.data.frame(clusters.combn)
  clusters.combn.f$jaccard_index<-jaccard.score
  passed.pairs<-subset(clusters.combn.f,jaccard_index >=min.jaccard_similarity)
  ###########run differentially expressed genes#######
  if(nrow(passed.pairs) > 0){
    n.diff.genes<-c()
    for(m in 1:nrow(passed.pairs)){
      cl.x1<-as.character(passed.pairs[m,1])
      cl.x2<-as.character(passed.pairs[m,2])
      cellsA<-getCellsFromClusters(object,cluster.type = cluster.type,cluster_id = cl.x1)
      cellsB<-getCellsFromClusters(object,cluster.type = cluster.type,cluster_id = cl.x2)
      
      my.diff.genes<-identifyDifferentialGenes(objectA = object,objectB = object,
                                               cellsA = cellsA,cellsB = cellsB,
                                               min.log2FC=min.log2FC,min.expFraction=min.expFraction) 
      if(nrow(my.diff.genes)==0){
        n.diff.genes<-append(n.diff.genes,0)
      }else{
        n.diff.genes<-append(n.diff.genes,nrow(my.diff.genes))
      }
    }
    passed.pairs$total_genes<-nrow(object)
    passed.pairs$number_of_diff_genes<-n.diff.genes
    passed.pairs$percent_of_diff_genes<-(n.diff.genes/nrow(object))*100
  }else{
    print("No cluster pairs has passed the jaccard_similarity cut-off. Please reduce the jaccard_similarity index cut-off!")
  }
  return(passed.pairs)
}

#' Merged clusters
#' @param  object  The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  merge_cluster_ids The vector of cluster ids that will be merged together e.g. c('cl1:cl2','cl9:cl10').
#' @export 

merge_clusters <- function(object,cluster.type=c("louvain","walktrap","kmeans"),merge_cluster_ids=c()){
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(length(merge_cluster_ids)==0){
    stop("Please add cluster names to be merged (e.g. merge_cluster_ids=c('cl1:cl2','cl9:cl10')")
  }
  
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  if(cluster.type=="walktrap"){
    clusters.info<-object@sc.clusters
    merged_walktrap_cluster<-as.character(clusters.info$walktrap_cluster)
    for(i in 1:length(merge_cluster_ids)){
      my.clusters<-merge_cluster_ids[i]
      cluster.names<-unlist(strsplit(my.clusters,":"))
      merged_walktrap_cluster[merged_walktrap_cluster %in% cluster.names]<-my.clusters
      merged_walktrap_cluster<-merged_walktrap_cluster
    }
    merge_walktrap.t<-table(merged_walktrap_cluster)
    merge_walktrap.t<-merge_walktrap.t[order(merge_walktrap.t,decreasing = T)]
    #########ranking clusters#############
    merge_walktrap.l<-list()
    for( y in 1:length(merge_walktrap.t)){
      cl.x<-names(merge_walktrap.t)[y]
      merge_walktrap.l[cl.x]<-y
    }
    merge_walktrap.ranked<-as.numeric(merge_walktrap.l[as.character(merged_walktrap_cluster)])
    #######################################
    merge.walktrap.cluster<-data.frame(merge_walktrap_cluster=paste("cl",merge_walktrap.ranked,sep=""))
    p<- qplot(0,0, data=merge.walktrap.cluster,colour=merge_walktrap_cluster)
    merge_walktrap_cluster_color<-ggplot_build(p)$data[[1]]$colour
    object@sc.clusters$merged_walktrap<-as.character(merge.walktrap.cluster[,1])
    object@sc.clusters$merged_walktrap_color<-merge_walktrap_cluster_color
    print("Done!")
  }else if(cluster.type=="louvain"){
    clusters.info<-object@sc.clusters
    merged_louvain_cluster<-as.character(clusters.info$louvain_cluster)
    for(i in 1:length(merge_cluster_ids)){
      my.clusters<-merge_cluster_ids[i]
      cluster.names<-unlist(strsplit(my.clusters,":"))
      merged_louvain_cluster[merged_louvain_cluster %in% cluster.names]<-my.clusters
      merged_louvain_cluster<-merged_louvain_cluster
    }
    merge_louvain.t<-table(merged_louvain_cluster)
    merge_louvain.t<-merge_louvain.t[order(merge_louvain.t,decreasing = T)]
    #########ranking clusters#############
    merge_louvain.l<-list()
    for( y in 1:length(merge_louvain.t)){
      cl.x<-names(merge_louvain.t)[y]
      merge_louvain.l[cl.x]<-y
    }
    merge_louvain.ranked<-as.numeric(merge_louvain.l[as.character(merged_louvain_cluster)])
    #######################################
    merge.louvain.cluster<-data.frame(merge_louvain_cluster=paste("cl",merge_louvain.ranked,sep=""))
    p<- qplot(0,0, data=merge.louvain.cluster,colour=merge_louvain_cluster)
    merge_louvain_cluster_color<-ggplot_build(p)$data[[1]]$colour
    object@sc.clusters$merged_louvain<-as.character(merge.louvain.cluster[,1])
    object@sc.clusters$merged_louvain_color<-merge_louvain_cluster_color
    print("Done!")
    
  }else if(cluster.type=="kmeans"){
    clusters.info<-object@knn_graph.kmeans.cluster
    merged_kmeans_cluster<-as.character(clusters.info$kmeans_cluster)
    for(i in 1:length(merge_cluster_ids)){
      my.clusters<-merge_cluster_ids[i]
      cluster.names<-unlist(strsplit(my.clusters,":"))
      merged_kmeans_cluster[merged_kmeans_cluster %in% cluster.names]<-my.clusters
      merged_kmeans_cluster<-merged_kmeans_cluster
    }
    merge_kmeans.t<-table(merged_kmeans_cluster)
    merge_kmeans.t<-merge_kmeans.t[order(merge_kmeans.t,decreasing = T)]
    #########ranking clusters#############
    merge_kmeans.l<-list()
    for( y in 1:length(merge_kmeans.t)){
      cl.x<-names(merge_kmeans.t)[y]
      merge_kmeans.l[cl.x]<-y
    }
    merge_kmeans.ranked<-as.numeric(merge_kmeans.l[as.character(merged_kmeans_cluster)])
    #######################################
    merge.kmeans.cluster<-data.frame(merge_kmeans_cluster=paste("cl",merge_kmeans.ranked,sep=""))
    p<- qplot(0,0, data=merge.kmeans.cluster,colour=merge_kmeans_cluster)
    merge_kmeans_cluster_color<-ggplot_build(p)$data[[1]]$colour
    object@knn_graph.kmeans.cluster$merged_kmeans<-as.character(merge.kmeans.cluster[,1])
    object@knn_graph.kmeans.cluster$merged_kmeans_color<-merge_kmeans_cluster_color
    print("Done!")
  }else{
    stop("Need a clustering-method name!")
  }
  ##########################
  print("Please rerun findMarkerGenes function.")
  
  assign(objName,object,envir=parent.frame())
  invisible(1)
}

#' Remove clusters
#' @param  object  The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  cluster_ids The vector of cluster ids that will be removed e.g. c('cl1:cl2','cl9:cl10').
#' @export 

remove_clusters <- function(object,cluster.type=c("louvain","walktrap","kmeans"),cluster_ids=c()){
  
  objName <- deparse(substitute(object))
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  if(length(cluster_ids)==0){
    stop("Please add cluster names (e.g. cluster_ids=c('cl1','cl2')")
  }
  
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  
  if(cluster.type=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster")]
  }else if(cluster.type=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster")]
  }else if(cluster.type=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","kmeans_cluster")]
  }else{
    stop("Need a clustering-method name!")
  }
  ##########################
  removed.clusters<-clusters.info[as.character(clusters.info[,2]) %in% cluster_ids,]
  meta.data<-get_cells_annotation(object)
  meta.data$IsPassed[meta.data$Cell %in% as.character(removed.clusters$Cell)]<-FALSE
  ##Update the metadata#####
  get_cells_annotation(object)<-meta.data
  ##Reset object############
  object@regressout_matrix<-""
  object@pca.result<-""
  object@tsne.result<-""
  object@umap.result<-""
  object@knn_graph.graph<-""
  object@knn_graph.layout<-""
  object@knn_graph.kmeans.cluster<-""
  object@igraph.graph<-""
  object@fa2_graph.layout<-""
  object@sc.clusters<-""
  object@marker.genes<-""
  ##########################
  print(paste(nrow(removed.clusters)," cells are removed. Please follow the pipeline steps starting from using the normalize_UMIs function!",sep=""))
  print("The cell's metadata is updated!.")
  assign(objName,object,envir=parent.frame())
  invisible(1)
}

#' Extracting selected clusters. 
#' @param  object  The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  cluster_ids The vector of cluster ids that will be removed e.g. c('cl1:cl2','cl9:cl10').
#' @export 
#' @return new SingCellaR object with selected clusters.

sub_clusters <- function(object,cluster.type=c("louvain","walktrap","kmeans","merged_louvain","merged_walktrap","merged_kmeans"),cluster_ids=c()){
  
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  
  if(length(cluster_ids)==0){
    stop("Please add cluster names (e.g. cluster_ids=c('cl1','cl2')")
  }
  
  clusters.info<-get_clusters(object)
  cluster.type <- match.arg(cluster.type)
  
  if(cluster.type=="walktrap"){
    clusters.info<-clusters.info[,c("Cell","walktrap_cluster")]
  }else if(cluster.type=="louvain"){
    clusters.info<-clusters.info[,c("Cell","louvain_cluster")]
  }else if(cluster.type=="kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","kmeans_cluster")]
  }else if(cluster.type=="merged_walktrap"){
    clusters.info<-clusters.info[,c("Cell","merged_walktrap")]
  }else if(cluster.type=="merged_louvain"){
    clusters.info<-clusters.info[,c("Cell","merged_louvain")]
  }else if(cluster.type=="merged_kmeans"){
    clusters.info<-get_knn_graph.kmeans.cluster(object)[,c("Cell","merged_kmeans")]
  }else{
    stop("Need a clustering-method name!")
  }
  ##########################
  selected.clusters<-clusters.info[as.character(clusters.info[,2]) %in% cluster_ids,]
  new.object<-object[,as.character(selected.clusters$Cell)]
  ##Reset object############
  new.object@dir_path_10x_matrix<-""
  new.object@sample_uniq_id<-""
  new.object@genes.info<-data.frame()
  new.object@meta.data<-data.frame()
  new.object@regressout_matrix<-""
  new.object@pca.result<-""
  new.object@tsne.result<-""
  new.object@umap.result<-""
  new.object@knn_graph.graph<-""
  new.object@knn_graph.layout<-""
  new.object@knn_graph.kmeans.cluster<-""
  new.object@igraph.graph<-""
  new.object@fa2_graph.layout<-""
  new.object@sc.clusters<-""
  new.object@marker.genes<-""
  ##########################
  print("Please redo the analysis steps starting from using the process_cells_annotation function!")
  return(new.object)
}

#' Export identified marker genes to table. 
#' @param  object  The SingCellaR object.
#' @param  cluster.type The clustering method name.
#' @param  n.TopGenes The number of top differential genes. Default 5
#' @param  min.log2FC The minimum log2FC.
#' @param  min.expFraction The minimum fraction of expressing cell frequency.
#' @param  write.to.file The output file name.
#' @export 
#' 
export_marker_genes_to_table<-function(object,cluster.type=c("walktrap","louvain","kmeans","merged_louvain","merged_walktrap","merged_kmeans"),
                                        n.TopGenes=5,min.log2FC=0.5,min.expFraction=0.3,write.to.file=""){
  
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
  if(write.to.file!=""){
    write.table(TopGenes.info,file=write.to.file,
                append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
  }else{
    return(TopGenes.info)
  }
}

#' Building AUCell gene rankings
#' @param  object  The SingCellaR object.
#' @param  AUCell_buildRankings.file The output file of AUCell rankings.
#' @param  nCores The number of cores for running AUCell. Default 1
#' @param  plotStats is logical. If TRUE, the AUCell stats will be displayed.
#' @export 
#' @importFrom AUCell AUCell_buildRankings

Build_AUCell_Rankings<-function(object,AUCell_buildRankings.file="",nCores=1, plotStats=FALSE){
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(AUCell_buildRankings.file==""){
    stop("Please add the output file name of AUCell_buildRankings!")
  }
  ########################################
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  my.select.data<-get_umi_count(object)
  my.use.dat<-my.select.data[,as.character(cells.used$Cell)]
  AUCell_rankings <- AUCell_buildRankings(as.matrix(my.use.dat), nCores=nCores, plotStats=plotStats)
  save(AUCell_rankings,file=AUCell_buildRankings.file)
  #########################################
}

#' Fast building AUCell gene rankings
#' @param  object  The SingCellaR object.
#' @param  AUCell_buildRankings.file The output file of AUCell rankings.
#' @param  nCores The number of cores for running AUCell. Default 1
#' @param  block.cell.size The number of genes in each block of AUCell processing. Default 5000
#' @param  plotStats is logical. If TRUE, the AUCell stats will be displayed.
#' @export 
#' 


Build_AUCell_Rankings_Fast<-function(object,AUCell_buildRankings.file="",nCores=1, block.cell.size=5000,plotStats=FALSE){
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(AUCell_buildRankings.file==""){
    stop("Please add the output file name of AUCell_buildRankings!")
  }
  ########################################
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  my.select.data<-get_umi_count(object)
  my.use.dat<-my.select.data[,as.character(cells.used$Cell)]
  ########################################
  step.size = block.cell.size
  ends <- c(seq(from=block.cell.size,to=ncol(my.use.dat),by=step.size),as.integer(ncol(my.use.dat)))
  starts <- seq(from=1,to=ncol(my.use.dat),by=step.size)
  window.cells<-data.frame(start=starts,end=ends)
  ########run AUCell per window of cells##
  s.start<-window.cells$start[1]
  s.end <-window.cells$end[1]
  AUCell_rankings.1 <- AUCell_buildRankings(as.matrix(my.use.dat[,s.start:s.end]), nCores=nCores, plotStats=plotStats,verbose=TRUE)
  Rank.1<-getRanking(AUCell_rankings.1)
  print("Completed ranking round : 1")
  for(i in 2:nrow(window.cells)){
    x.start<-window.cells$start[i]
    x.end <-window.cells$end[i]
    AUCell_rankings <- AUCell_buildRankings(as.matrix(my.use.dat[,x.start:x.end]), nCores=nCores, plotStats=plotStats,verbose=TRUE)
    Rank.i<-getRanking(AUCell_rankings)
    Rank.1<-cbind(Rank.1,Rank.i)
    print(paste("Completed ranking round : ",i,sep=""))
  }
  #########Create auCellResults object#####
  names(dimnames(Rank.1)) <- c("genes", "cells")
  AUCell_rankings<-new("aucellResults",SummarizedExperiment::SummarizedExperiment(assays=list(ranking=Rank.1)))
  #########################################
  save(AUCell_rankings,file=AUCell_buildRankings.file)
  #########################################
}

#' Run AUCell
#' @param  object  The SingCellaR object.
#' @param  AUCell_buildRankings.file The input file name of AUCell rankings.
#' @param  geneSets.gmt.file The GMT file name that contains gene sets.
#' @param  n.random.genes The number of a random gene set.
#' @export 
#' 


Run_AUCell<-function(object,AUCell_buildRankings.file="",geneSets.gmt.file="",n.random.genes=500){
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  if(AUCell_buildRankings.file==""){
    stop("Please input the AUCell_buildRankings file")
  }
  if(geneSets.gmt.file==""){
    stop("Please input the AUCell_buildRankings file")
  }
  genes<-get_genes_metadata(object)
  genes<-subset(genes,IsExpress==T)
  gene.names<-as.character(rownames(genes))
  random.genes<-list()
  for(i in 1:10){
    r.index<-sample(1:length(gene.names), n.random.genes, replace = FALSE)
    sample.genes<-gene.names[r.index]
    h<-paste("random.genes",i,sep="")
    random.genes[[h]]<-sample.genes
  }
  
  geneSets<-get_gmtGeneSets(geneSets.gmt.file)
  ############################################################
  AUCell_buildRankings<-local(get(load(AUCell_buildRankings.file)))
  cells_AUC <- AUCell_calcAUC(geneSets, AUCell_buildRankings)
  geneSets.AUC<-as.data.frame(t(getAUC(cells_AUC)))
  #################RUN AUCell for random genes################
  cells_AUC.random <- AUCell_calcAUC(random.genes, AUCell_buildRankings)
  geneSets.AUC.random<-as.data.frame(t(getAUC(cells_AUC.random)))
  ############################################################
  geneSets.AUC$random_gene<-rowMeans(geneSets.AUC.random)
  #############################################################
  return(geneSets.AUC)
 }