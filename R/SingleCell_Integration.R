#' Preprocessing for data integration
#' @param  object The SingCellaR object.
#' @param  mitochondiral_genes_start_with The unique prefix alphabets for mitocondrial genes such as 'MT-' for human and 'mt-' for mouse genes.
#' @export 
#' 

preprocess_integration <- function(object,mitochondiral_genes_start_with="MT-"){
	
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR_Int")){
		stop("Need to initialize the SingCellaR_Int object")
	}
	
	data.dir<-object@dir_path_SingCellR_object_files
	
	if (! dir.exists(data.dir)){
		stop("Directory provided does not exist")
	}
	input.files<-object@SingCellR_object_files
	if(length(input.files)==1){
		stop("Need more SingCellaR object files as the input!")
	}else{
		#print("Combining Genes..")
		
		pb <- txtProgressBar(max = length(input.files), style = 3)
		
		#EXP.genes<-list()
		#for(n in 1:length(input.files)){
		#  Sys.sleep(0.5);
		#  file.name  <-input.files[n]
		#  local.obj <-local(get(load(file=file.name)))
		#  local.exp.genes<-get_genes_metadata(local.obj)
		#  local.exp.genes<-subset(local.exp.genes,IsExpress==TRUE)
		#  local.genes<-rownames(local.exp.genes)
		#  EXP.genes[[n]]<-local.genes
		#  setTxtProgressBar(pb, pb$getVal()+1)
		#}
		#UNION.genes<-unique(unlist(EXP.genes))
		###############################################
		print("Combining variable/marker genes, clusters info, and UMIs..")
		###############################################
		file.name.1  <-paste(data.dir,input.files[1],sep="/")
		
		local.obj.1 <-local(get(load(file=file.name.1)))
		local.cells.1<-get_cells_annotation(local.obj.1)
		local.cells.1<-subset(local.cells.1,IsPassed==T)
		#local.obj.1<-local.obj.1[UNION.genes,as.character(local.cells.1$Cell)]
		#local.obj.1<-local.obj.1[,as.character(local.cells.1$Cell)]
		local.vargenes.1<-get_genes_metadata(local.obj.1)
		local.vargenes.1<-subset(local.vargenes.1,IsExpress==T & IsVarGenes==T)
		local.vargenes.1<-rownames(local.vargenes.1)
		
		int.umi<-get_umi_count(local.obj.1)
		int.umi<-int.umi[,as.character(local.cells.1$Cell)]
		set.name<-rep(1,ncol(int.umi))
		
		var.genes<-list()
		var.genes[[1]]<-local.vargenes.1
		########get cluster and marker genes############
		cluster.info<-list()
		cluster.info[[1]]<-get_clusters(local.obj.1)
		
		marker.genes<-list()
		marker.genes [[1]]<-get_marker_genes(local.obj.1)
		
		
		for(i in 2:length(input.files)){
			Sys.sleep(0.5);
			file.name  <-paste(data.dir,input.files[i],sep="/")
			
			local.obj <-local(get(load(file=file.name)))
			local.cells<-get_cells_annotation(local.obj)
			local.cells<-subset(local.cells,IsPassed==T)
			#local.obj<-local.obj[UNION.genes,as.character(local.cells$Cell)]
			#local.obj<-local.obj[,as.character(local.cells$Cell)]
			#######get variable genes####
			local.vargenes<-get_genes_metadata(local.obj)
			local.vargenes<-subset(local.vargenes,IsExpress==T & IsVarGenes==T)
			local.vargenes<-rownames(local.vargenes)
			var.genes[[i]]<-local.vargenes
			#############################
			cluster.info[[i]]<-get_clusters(local.obj)
			marker.genes [[i]]<-get_marker_genes(local.obj)
			#############################
			local.umi<-get_umi_count(local.obj)
			local.umi<-local.umi[,as.character(local.cells$Cell)]
			int.umi<-cbind(int.umi,local.umi)
			set.name.x<-rep(i,ncol(local.umi))
			set.name<-append(set.name,set.name.x)
			setTxtProgressBar(pb, pb$getVal()+1)
		}
	}
	
	sce<-SingleCellExperiment(assays = list(counts = int.umi))
	object <- new("SingCellaR_Int", sce, dir_path_SingCellR_object_files=data.dir,SingCellR_object_files=input.files)
	object@Variable.genes<-var.genes
	object@marker.genes<-marker.genes
	object@SingCellaR.individual.clusters<-cluster.info
	################################
	process_cells_annotation(object,mitochondiral_genes_start_with=mitochondiral_genes_start_with)
	object@meta.data<-cbind(get_cells_annotation(object),data.frame(data_set=set.name))
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("The integrated sparse matrix is created.")
	
}

preprocess_integration_for_ADT <- function(object){
	
	objName <- deparse(substitute(object))
	
	if(!is(object,"SingCellaR_Int")){
		stop("Need to initialize the SingCellaR_Int object")
	}
	data.dir<-object@dir_path_SingCellR_object_files
	
	if (! dir.exists(data.dir)){
		stop("Directory provided does not exist")
	}
	
	input.files<-object@SingCellR_object_files
	if(length(input.files)==1){
		stop("Need more SingCellaR object files as the input!")
	}else{
		pb <- txtProgressBar(max = length(input.files), style = 3)
		###############################################
		print("Checking number of rows..")
		antibody.name<-c()
		set.name<-c()
		for(i in 1:length(input.files)){
			Sys.sleep(0.5);
			file.name  <-paste(data.dir,input.files[i],sep="/")
			local.obj <-local(get(load(file=file.name)))
			#############################
			#############################
			local.umi<-get_umi_count(local.obj)
			antibody.name<-append(antibody.name,rownames(local.umi))
			setTxtProgressBar(pb, pb$getVal()+1)
		}
		antibody.unique<-unique(antibody.name)
		int.umi <- Matrix(nrow = length(antibody.unique), ncol = 1, data = 0, sparse = TRUE)
		int.umi <- as(int.umi, "dgTMatrix") 
		rownames(int.umi)<-antibody.unique
		################################
		print("Combining UMIs..")
		################################
		for(i in 1:length(input.files)){
			Sys.sleep(0.5);
			file.name  <-paste(data.dir,input.files[i],sep="/")
			
			local.obj <-local(get(load(file=file.name)))
			#############################
			#############################
			local.umi<-get_umi_count(local.obj)
			int.umi<-cbindX(as.matrix(int.umi),as.matrix(local.umi))
			set.name.x<-rep(i,ncol(local.umi))
			set.name<-append(set.name,set.name.x)
			setTxtProgressBar(pb, pb$getVal()+1)
		}
	}
	int.umi<-int.umi[,-c(1)]
	int.umi[is.na(int.umi)==T]<-0
	int.umi<- as(int.umi, "dgTMatrix")
	sce<-SingleCellExperiment(assays = list(counts = int.umi))
	object <- new("SingCellaR_Int", sce, dir_path_SingCellR_object_files=data.dir,SingCellR_object_files=input.files)
	################################
	process_cells_annotation(object)
	object@meta.data<-cbind(get_cells_annotation(object),data.frame(data_set=set.name,IsPassed=TRUE))
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("The integrated sparse matrix is created.")	
}
#' Run scanorama integration
#' @param  object The SingCellaR_Int object.
#' @param  useCombinedVarGenesFromIndividualSample is logical. If TRUE, all combined variable genes from individual sample will be used. If FALSE, highly variable genes from the whole combined samples will be used.
#' @export 
#' 

runScanorama <- function(object,useCombinedVarGenesFromIndividualSample=F){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR_Int")){
    stop("Need to initialize the SingCellaR_Int object")
  }
 
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  ##############################
  if(useCombinedVarGenesFromIndividualSample==T){
    VarGenes<-unlist(object@Variable.genes)
    selected.genes<-unique(VarGenes)
  }else{
      genes_info<-get_genes_metadata(object)
      if(c("IsVarGenes") %in% colnames(genes_info)==FALSE){
        stop("No variable genes are found, please identify variable genes!")
      }else{
        selected.genes<-subset(genes_info,IsVarGenes==TRUE)
        selected.genes<-as.character(rownames(selected.genes))
      }
  }
  ##############################
  data.set.ids<-unique(object@meta.data$data_set)
  cells.info<-get_cells_annotation(object)
  cells.info<-subset(cells.info,IsPassed ==TRUE)
  cells.dat<-log1p(get_normalized_umi(object))
  ##############################
  datasets.mat<-list()
  genes_list<-list()
  
  cell.info.update<-data.frame()
  for(i in 1:length(data.set.ids)){
    my.cells<-subset(cells.info,data_set==i)
    each.set.m<-cells.dat[,colnames(cells.dat) %in% as.character(my.cells$Cell)]
    each.set.m<-each.set.m[rownames(each.set.m) %in% selected.genes,]
    datasets.mat[[i]]<-as.matrix(t(as.matrix(each.set.m)))
    genes_list[[i]]<-as.character(selected.genes)
    cell.info.update<-rbind(cell.info.update,my.cells)
  }
  #########import scanorama##########
  if(py_module_available("scanorama")==FALSE){
    stop("The scanorama python module is not installed!. Please install using pip ('pip install scanorama')")
  }
  scanorama <- import('scanorama')
  #########Scanorama integration#####
  print("This process will take time and requires large RAM depending on the number of cells in your integration.")
  integrated.corrected.data <- scanorama$correct(datasets.mat, genes_list,return_dimred=TRUE, return_dense=TRUE)
  ###################################
  data.scanorama<-integrated.corrected.data[[2]][[1]]
  for(j in 2:length(data.set.ids)){
    data.scanorama<-rbind(data.scanorama,integrated.corrected.data[[2]][[j]])
  }
  colnames(data.scanorama)<-integrated.corrected.data[[3]]
  rownames(data.scanorama)<-as.character(cell.info.update$Cell)
  object@Scanorama.integrative.matrix<-t(data.scanorama)
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Scanorama analysis is done!.")
}

#' Run harmony integration
#' @param  object The SingCellaR_Int object.
#' @param  n.dims.use The number of PCA dimensions used for the input for harmony.
#' @param  covariates The unwanted source of variations (e.g. batch, sample_id, etc).
#' @param  harmony.theta The harmony theta parameter. Deafult NULL. This is the diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
#' @param  harmony.sigma The harmony sigma parameter. Default 0.1. This parameter is the width of soft kmeans clusters. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering.
#' @param  harmony.nclust The harmony nclust parameter. Default NULL
#' @param  harmony.tau The harmony tau parameter. Default 0. The parameter is for the protection against overclustering small datasets with large ones. tau is the expected number of cells per cluster.
#' @param  harmony.block.size The harmony block.size parameter. Default 0.05. What proportion of cells to update during clustering. Between 0 to 1. Larger values may be faster but less accurate
#' @param  harmony.max.iter The harmony max iteration parameter. Default 10. Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step.
#' @param  harmony.max.iter.cluster The harmony max iteration cluster. Default 200. Maximum number of rounds to run clustering at each round of Harmony.
#' @param  harmony.epsilon.cluster The harmony epsilon.cluster parameter. Default 1e-05. Convergence tolerance for clustering round of Harmony. Set to -Inf to never stop early.
#' @param  harmony.epsilon.harmony The harmony epsilon parameter. Default 1e-04. Convergence tolerance for Harmony. Set to -Inf to never stop early.
#' @param  n.seed The set seed number.
#' @export 
#' 

runHarmony <- function(object,n.dims.use=30,covariates=c("data_set"),harmony.theta = NULL,harmony.sigma = 0.1,harmony.nclust = NULL,
                       harmony.tau = 0,harmony.block.size = 0.05,harmony.max.iter = 10,harmony.max.iter.cluster = 200,
                       harmony.epsilon.cluster = 1e-05,harmony.epsilon.harmony = 1e-04,n.seed=1){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR_Int")){
    stop("Need to initialize the SingCellaR_Int object")
  }
  ######################################################
  orig.pca<-get_pca.result(object)$x[,1:n.dims.use]
  cell.meta<-get_cells_annotation(object)
  cell.meta<-subset(cell.meta,IsPassed==TRUE)
  ######################################################
  ###Run Harmony########################################
  set.seed(seed=n.seed)
  ######################################################
  harmony_embeddings <- harmony::HarmonyMatrix(orig.pca, 
                                               cell.meta,  
                                               covariates,
                                               theta = harmony.theta,
                                               sigma = harmony.sigma,
                                               nclust = harmony.nclust,
                                               tau = harmony.tau,
                                               block.size = harmony.block.size,
                                               max.iter.harmony = harmony.max.iter,
                                               max.iter.cluster = harmony.max.iter.cluster,
                                               epsilon.cluster = harmony.epsilon.cluster,
                                               epsilon.harmony = harmony.epsilon.harmony,
                                               do_pca = FALSE, verbose=TRUE)
  ######################################################
  object@Harmony.embeddings<-harmony_embeddings
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Harmony analysis is done!.")
}

#' Run Seurat integration
#' @param  object The SingCellaR_Int object.
#' @param  Seurat.split.by The indicated feature for splitting samples for the integration.
#' @param  n.dims.use The number of PCA dimensions used for the input for Seurat.
#' @param  Seurat.variablegenes.method The method for variable genes selection. Default 'vst'
#' @param  Seurat.variablegenes.number The number of highly varible genes. Default 2000
#' @param  use.SingCellaR.varGenes is logical. If TRUE, the highly variable genes identified by SingCellaR will be used.
#' @param  k.anchor Seurat's k anchor parameter. Default 5
#' @param  k.filter Seurat's k filter parameter. Default 200
#' @param  k.score Seurat's k score parameter. Default 30
#' @export 
#' 

runSeuratIntegration <- function(object,Seurat.metadata="",Seurat.split.by="",n.dims.use=30,
                                 Seurat.variablegenes.method = "vst",Seurat.variablegenes.number=2000,
                                 use.SingCellaR.varGenes=F,k.anchor=5,k.filter=200,k.score=30){
	
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR_Int")){
		stop("Need to initialize the SingCellaR_Int object")
	}
	######################################################
	Seuarat_Obj <- Seurat::CreateSeuratObject(get_umi_count(object), meta.data = Seurat.metadata)
	Seurat.list <- Seurat::SplitObject(Seuarat_Obj, split.by = Seurat.split.by)
	
	pb <- txtProgressBar(max = length(Seurat.list), style = 3)
	print("Process each Seurat object..")
	
	if(use.SingCellaR.varGenes==F){
		for (i in 1:length(Seurat.list)) {
			Sys.sleep(0.5);
			Seurat.list[[i]] <- Seurat::NormalizeData(Seurat.list[[i]], verbose = FALSE)
			Seurat.list[[i]] <- Seurat::FindVariableFeatures(Seurat.list[[i]], selection.method = Seurat.variablegenes.method, 
					nfeatures = Seurat.variablegenes.number, verbose = FALSE)
			setTxtProgressBar(pb, pb$getVal()+1)
		}
	}else{
		for (i in 1:length(Seurat.list)) {
			Sys.sleep(0.5);
			gene.info<-get_genes_metadata(object)
			gene.info<-subset(gene.info,IsVarGenes==TRUE)
			var.genes<-as.character(rownames(gene.info))
			Seurat.list[[i]] <- NormalizeData(Seurat.list[[i]], verbose = FALSE)
			Seurat.list[[i]]@assays$RNA@var.features <- var.genes
			setTxtProgressBar(pb, pb$getVal()+1)
		}
	}
	######Seurat Integration###############################
	print("This process will take time and requires large RAM depending on the number of cells in your integration.")
	reference.list <- Seurat.list
	my.anchors     <- Seurat::FindIntegrationAnchors(object.list = reference.list, dims = 1:n.dims.use,k.anchor=k.anchor,k.filter= k.filter,k.score=k.score)
	my.integrated  <- Seurat::IntegrateData(anchorset = my.anchors, dims = 1:n.dims.use)
	my.integrated  <- Seurat::ScaleData(my.integrated, verbose = FALSE)
	my.integrated  <- Seurat::RunPCA(my.integrated, npcs = 100, verbose = FALSE)
	######################################################
	object@Seurat.embeddings <- my.integrated@reductions$pca@cell.embeddings
	######################################################
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("Seurat integrative analysis is done!.")
}

runSeuratIntegration_with_rpca <- function(object,Seurat.metadata="",Seurat.split.by="",n.dims.use=30,
		Seurat.variablegenes.method = "vst",Seurat.variablegenes.number=2000,
		use.SingCellaR.varGenes=F){
	
	objName <- deparse(substitute(object))
	if(!is(object,"SingCellaR_Int")){
		stop("Need to initialize the SingCellaR_Int object")
	}
	
	Seuarat_Obj <- Seurat::CreateSeuratObject(get_umi_count(object), meta.data = Seurat.metadata)
	Seurat.list <- Seurat::SplitObject(Seuarat_Obj, split.by = Seurat.split.by)
	
	pb <- txtProgressBar(max = length(Seurat.list), style = 3)
	print("Process each Seurat object..")
	
	if(use.SingCellaR.varGenes==F){
		for (i in 1:length(Seurat.list)) {
			Sys.sleep(0.5);
			Seurat.list[[i]] <- Seurat::NormalizeData(Seurat.list[[i]], verbose = FALSE)
			Seurat.list[[i]] <- Seurat::FindVariableFeatures(Seurat.list[[i]], selection.method = Seurat.variablegenes.method, 
					nfeatures = Seurat.variablegenes.number, verbose = FALSE)
			setTxtProgressBar(pb, pb$getVal()+1)
		}
	}else{
		for (i in 1:length(Seurat.list)) {
			Sys.sleep(0.5);
			gene.info<-get_genes_metadata(object)
			gene.info<-subset(gene.info,IsVarGenes==TRUE)
			var.genes<-as.character(rownames(gene.info))
			Seurat.list[[i]] <- NormalizeData(Seurat.list[[i]], verbose = FALSE)
			Seurat.list[[i]]@assays$RNA@var.features <- var.genes
			setTxtProgressBar(pb, pb$getVal()+1)
		}
	}
	######Seurat Integration###############################
	print("This process will take time and requires large RAM depending on the number of cells in your integration.")
	
	features <- SelectIntegrationFeatures(object.list = Seurat.list)
	
	Seurat.list <- lapply(X = Seurat.list, FUN = function(x) {
				x <- ScaleData(x, features = features, verbose = FALSE)
				x <- RunPCA(x, features = features, verbose = FALSE)
			})
	
	anchors <- FindIntegrationAnchors(object.list = Seurat.list, 
			reference = 1:length(Seurat.list), 
			reduction = "rpca", 
			dims = 1:n.dims.use)
	integrated_obj <- IntegrateData(anchorset = anchors, dims = 1:n.dims.use)
	integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
	integrated_obj <- RunPCA(integrated_obj, verbose = FALSE)
	######################################################
	object@Seurat.embeddings <- integrated_obj@reductions$pca@cell.embeddings
	######################################################
	assign(objName,object,envir=parent.frame())
	invisible(1)
	print("Seurat integrative analysis is done!.")
}

#' Run combat integration
#' @param  object The SingCellaR_Int object.
#' @param  use.reduced_dim is logical. If TRUE, the reduced dimensions from PCA or nnmf will be used as the input. If FALSE, the normalized UMI will be used. 
#' @param  dim_reduction_method The method name for dimensional reduction.
#' @param  n.dims.use The number of dimensions to be used. Default 30
#' @param  batch_identifier The indicated feature name to separate the batch.
#' @export 
#' 

runCombat<-function(object,use.reduced_dim=T,
                    dim_reduction_method=c("pca","nnmf"),
                    n.dims.use=30,batch_identifier=""){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR")){
    stop("Need to initialize the SingCellaR object")
  }
  dim_reduction_method<-match.arg(dim_reduction_method)
  ######################################
  cells.info<-get_cells_annotation(object)
  cells.used<-subset(cells.info,IsPassed==TRUE)
  rownames(cells.used)<-cells.used$Cell
  drop <- c("Cell")
  pheno<-cells.used[!(colnames(cells.info) %in% drop)]
  ######Run SVA::Combat#################
  modcombat = model.matrix(~1, data=pheno)
  batch<-cells.used[,batch_identifier]
  ######Get embeddings###################
  if(use.reduced_dim==TRUE){
    if(dim_reduction_method=="pca"){
      res.pca<-get_pca.result(object)
      my.embeddings<-as.matrix(res.pca$x[,1:n.dims.use])
    }else if(dim_reduction_method=="nnmf"){
      res.nnmf<-get_nnmf.result(object)
      my.embeddings<-t(res.nnmf$H)[,1:n.dims.use]
    }
    combat_embeddings = sva::ComBat(dat= t(my.embeddings), batch=batch, par.prior=TRUE, prior.plots=FALSE)
    object@Combat.embeddings<-t(combat_embeddings)
    print("Combat analysis is done!")
    print("The batch-free embeddings is now in the Combat.embeddings slot!")
  }else if(use.reduced_dim==FALSE){
    
    genes_info<-get_genes_metadata(object)
    selected.genes<-subset(genes_info,IsVarGenes==TRUE)
    #####################################
    my.select.data<-get_normalized_umi(object)
    SM<-log1p(my.select.data[rownames(my.select.data) %in% rownames(selected.genes),as.character(cells.used$Cell)])
    object@regressout_matrix = as.big.matrix(sva::ComBat(dat= as.matrix(SM), batch=batch, par.prior=TRUE, prior.plots=FALSE))
    print("Combat analysis is done!")
    print("The batch-free matrix is now in the regressout_matrix slot!")
    print("Please continue PCA analysis with 'use.regressout.data=T'.")
    
  }
  #######################################
  assign(objName,object,envir=parent.frame())
  invisible(1)
}

#' Run supervised harmony integration
#' @param  object The SingCellaR_Int object.
#' @param  n.dims.use The number of PCA dimensions used for the input for harmony. Default 30
#' @param  fGSEA.minSize The cutoff minimum number of genes in each gene set. Gene set that contains the number of genes lower than this number will be excluded. Default 10
#' @param  fGSEA.maxSize The cutoff maxiumum number of genes in each gene set. Gene set that contains the number of genes higher than this number will be excluded. Default 500
#' @param  fGSEA.eps The eps paramenter for fgsea, this parameter sets the boundary for calculating the p value. Default 1e-10
#' @param  hcl.height.cutoff The cutree cutoff value for hierarchical clustering. Default 0.25
#' @param  covariates The unwanted source of variations (e.g. batch, sample_id, etc).
#' @param  harmony.sigma The harmony sigma parameter. Default 0.1
#' @param  harmony.tau The harmony tau parameter. Default 0
#' @param  harmony.block.size The harmony block.size parameter. Default 0.05
#' @param  harmony.max.iter The harmony max iteration parameter. Default 10
#' @param  harmony.max.iter.cluster The harmony max iteration cluster. Default 200
#' @param  harmony.epsilon.cluster The harmony epsilon.cluster parameter. Default 1e-05
#' @param  harmony.epsilon.harmony The harmony epsilon parameter. Default 1e-04
#' @param  n.seed The set seed number.
#' @export 
#' 

runSupervised_Harmony <- function(object,n.dims.use=30,fGSEA.minSize=10,fGSEA.maxSize=500,fGSEA.eps = 1e-10,
                                  hcl.height.cutoff=0.25,covariates=c("data_set"),harmony.sigma = 0.1,
                                  harmony.tau = 0,harmony.block.size = 0.05,harmony.max.iter = 10,harmony.max.iter.cluster = 200,
                                  harmony.epsilon.cluster = 1e-05,harmony.epsilon.harmony = 1e-04,n.seed=1){
  
  objName <- deparse(substitute(object))
  if(!is(object,"SingCellaR_Int")){
    stop("Need to initialize the SingCellaR_Int object")
  }
  input.files<-object@SingCellR_object_files
  genes.db<-list()
  for(i in 1:length(input.files)){
    marker.genes<-object@marker.genes[[i]]
    for(j in 1:length(marker.genes[["louvain"]])){
      cl.marker.genes<-marker.genes[["louvain"]][[j]]
      my.genes<-as.character(cl.marker.genes$Gene)
      cl.cluster<-names(marker.genes[["louvain"]][j])
      my.title<-paste(cl.cluster,i,sep="_")
      genes.db[[my.title]]<- my.genes
    }
  }
  set.seed(seed=n.seed)
  #######fGSEA parameters#####
  gsea.minSize=fGSEA.minSize
  gsea.maxSize=fGSEA.maxSize
  #gsea.nperm=fGSEA.nperm
  ############################
  gsea.results<-data.table()
  
  for(i in 1:length(input.files)){
    marker.genes<-object@marker.genes[[i]]
    for(j in 1:length(marker.genes[["louvain"]])){
      
      cl.diff.genes<-marker.genes[["louvain"]][[j]]
      i.index<-which(cl.diff.genes$combined.pval==0)
      
      if(length(i.index) > gsea.minSize){
        cl.diff.genes$combined.pval[i.index]<-min(cl.diff.genes$combined.pval[-c(i.index)])
      }
      cl.diff.genes$prerank.genes<- (-log10(cl.diff.genes$combined.pval))*cl.diff.genes$log2FC
      #########remove -inf genes#########
      inf.index<-which(is.infinite(cl.diff.genes$prerank.genes)==T)
      if(length(inf.index) > 0){
        cl.diff.genes<-cl.diff.genes[-c(inf.index)]
      } 
      ###################################
      prerank.genes<-cl.diff.genes$prerank.genes
      names(prerank.genes)<-cl.diff.genes$Gene
      ###################################
      cl.cluster<-names(marker.genes[["louvain"]][j])
      my.title<-paste(cl.cluster,i,sep="_")
      
      rm.index<-which(names(genes.db) %in% my.title)
      y.genes.db<-genes.db[-rm.index]
      
      print(paste("Processing fGSEA for: ",my.title,sep=""))
      fgseaRes <- fgseaMultilevel(pathways = y.genes.db, 
                        stats = prerank.genes,
                        minSize=gsea.minSize,
                        maxSize=gsea.maxSize,
                        eps=fGSEA.eps)
      
      fgseaRes<-fgseaRes[fgseaRes$ES > 0]
      fgseaRes$cluster<-my.title
      fgseaRes.f<-fgseaRes[,c("pathway","pval","padj","ES","NES","cluster")]
      gsea.results<-rbind(gsea.results,fgseaRes.f)
      
    }
  }
  ####################
  cluster.matching.score<-dcast(gsea.results, cluster ~ pathway, value.var="pval")
  cluster.matching.score<-data.frame(cluster.matching.score)
  rownames(cluster.matching.score)<-cluster.matching.score$cluster
  cluster.matching.score<-cluster.matching.score[,-c(1)]
  cluster.matching.score[is.na(cluster.matching.score)==T]<-1
  cluster.matching.score<--log10(cluster.matching.score+0.0001)
  cluster.matching.score<-cor(cluster.matching.score)
  ####################
  dendrogram.cutoff_height<-hcl.height.cutoff
  hr <- hclust(dist(cluster.matching.score))
  group.clusters<-cutree(hr,h=max(hr$height)*dendrogram.cutoff_height)
  ####################
  group.cluster.table<-data.frame(y_cluster=names(group.clusters),
                                  common_cluster=paste("cl",group.clusters,sep=""))
  ####################
  total_s<-object@SingCellaR.individual.clusters
  
  if(length(total_s)==0){
    stop("Required Louvain's identify clusters for each sample!")
    
  }else{
    individual.cl<-data.frame()
    for(k in 1:length(total_s)){
      ind.cluster<-total_s[[k]]
      ind.cluster$data_set<-k
      ind.cluster$pre_cluster<-paste(ind.cluster$louvain_cluster,ind.cluster$data_set,sep="_")
      individual.cl<-rbind(individual.cl,ind.cluster)
    }
  }
  ######################
  combined_info<-merge(individual.cl,group.cluster.table,by.x="pre_cluster",by.y="y_cluster",all.x=T)
  combined_info<-combined_info[,-which(colnames(combined_info) =="data_set")]
  class(combined_info$common_cluster)<-"character"
  combined_info$common_cluster[is.na(combined_info$common_cluster)==T]<-combined_info$pre_cluster[is.na(combined_info$common_cluster)==T]
  
  cell_anno<-get_cells_annotation(object)
  combined_info<-merge(combined_info,cell_anno,by.x="Cell",by.y="Cell",all.x=TRUE)
  rownames(combined_info)<-combined_info$Cell
  ######################################################
  orig.pca<-get_pca.result(object)$x[,1:n.dims.use]
  cell.meta<-combined_info[rownames(orig.pca),]
  ######################################################
  ccm.clusters<-unique(combined_info$common_cluster)
  n.uniq<-length(grep("cl",ccm.clusters))
  n.common<-length(ccm.clusters)-n.uniq
  show.info.txt<-paste("Identify :", n.common, "matched and", n.uniq, "uniquely un-matched clusters!")
  print(show.info.txt)
  print(paste(length(ccm.clusters)," cluster centers will be used as the pre-assigned clusters for harmony!",sep=""))
  ######################################################  
  R.mat<-matrix(0,nrow(orig.pca), ncol=length(ccm.clusters))
  rownames(R.mat)<-rownames(orig.pca)
  
  for(j in 1:length(ccm.clusters)){
    ccm.j<-ccm.clusters[j]
    sub_info<-subset(cell.meta,common_cluster==ccm.j)
    sub_cells_names<-as.character(sub_info$Cell)
    R.mat[rownames(R.mat) %in% sub_cells_names,j]<-1
  }
  ################Run Harmony###########################
  ######################################################
  harmony_embeddings <- harmony::HarmonyMatrix(orig.pca, 
                                               cell.meta,  
                                               covariates,
                                               sigma = harmony.sigma,
                                               nclust = ncol(R.mat),
                                               tau = harmony.tau,
                                               block.size = harmony.block.size,
                                               max.iter.harmony = harmony.max.iter,
                                               max.iter.cluster = harmony.max.iter.cluster,
                                               epsilon.cluster = harmony.epsilon.cluster,
                                               epsilon.harmony = harmony.epsilon.harmony,
                                               cluster_prior=t(R.mat),do_pca = FALSE, verbose=TRUE)
  ######################################################
  object@SupervisedHarmony.embeddings<-harmony_embeddings
  assign(objName,object,envir=parent.frame())
  invisible(1)
  print("Supervised harmony analysis is done!.")
}