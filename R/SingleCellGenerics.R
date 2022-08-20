#' Get path directory of 10x data
#' @param  object The SingCellaR object.
#' @export 

setGeneric(
		name="get_dir_path_10x_matrix",
		def=function(object){
			standardGeneric("get_dir_path_10x_matrix")
		})
setMethod("get_dir_path_10x_matrix",
		signature(object = "SingCellaR"),
		function (object){
			object@dir_path_10x_matrix
		})
########
#' Get unique sample ID
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_sample_uniq_id",
		def=function(object){
			standardGeneric("get_sample_uniq_id")
		})
setMethod("get_sample_uniq_id",
		signature(object = "SingCellaR"),
		function (object){
			object@sample_uniq_id
		})
########
#' Get UMI count matrix
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_umi_count",
		def=function(object){
			standardGeneric("get_umi_count")
		})
setMethod("get_umi_count",
		signature(object = "SingCellaR"),
		function (object){
			assay(object,"counts")
		})

########
#' Get normalized UMI count matrix
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_normalized_umi",
		def=function(object){
			standardGeneric("get_normalized_umi")
		})
setMethod("get_normalized_umi",
		signature(object = "SingCellaR"),
		function (object){
		  assay(object,"normalized.umi")
		})
########
#' Get regressed out log-scaled matrix
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
  name="get_regressout_log_data",
  def=function(object){
    standardGeneric("get_regressout_log_data")
  })
setMethod("get_regressout_log_data",
          signature(object = "SingCellaR"),
          function (object){
            object@regressout_matrix
          })
########
#' Get cell annotation
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_cells_annotation",
		def=function(object){
			standardGeneric("get_cells_annotation")
		})
setGeneric(
		name="get_cells_annotation<-",
		def=function(object,value){
			standardGeneric("get_cells_annotation<-")
		})
setMethod("get_cells_annotation",
		signature(object = "SingCellaR"),
		function (object){
			object@meta.data
		})
setReplaceMethod(
		f="get_cells_annotation",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,meta.data=value)
		})
########
#' Get gene's metadata
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_genes_metadata",
		def=function(object){
			standardGeneric("get_genes_metadata")
		})
setGeneric(
		name="get_genes_metadata<-",
		def=function(object,value){
			standardGeneric("get_genes_metadata<-")
		})
setMethod("get_genes_metadata",
		signature(object = "SingCellaR"),
		function (object){
			object@genes.info
		})
setReplaceMethod(
		f="get_genes_metadata",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,genes.info=value)
		})
########
#' Get PCA analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_pca.result",
		def=function(object){
			standardGeneric("get_pca.result")
		})
setGeneric(
		name="get_pca.result<-",
		def=function(object,value){
			standardGeneric("get_pca.result<-")
		})
setMethod("get_pca.result",
		signature(object = "SingCellaR"),
		function (object){
			object@pca.result
		})
setReplaceMethod(
		f="get_pca.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,pca.result=value)
		})
########
#' Get NNMF analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
  name="get_nnmf.result",
  def=function(object){
    standardGeneric("get_nnmf.result")
  })
setGeneric(
  name="get_nnmf.result<-",
  def=function(object,value){
    standardGeneric("get_nnmf.result<-")
  })
setMethod("get_nnmf.result",
          signature(object = "SingCellaR"),
          function (object){
            object@nnmf.result
          })
setReplaceMethod(
  f="get_nnmf.result",
  signature="SingCellaR",
  definition=function(object,value){
    initialize(object,nnmf.result=value)
  })
########
#' Get LSI analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_lsi.result",
		def=function(object){
			standardGeneric("get_lsi.result")
		})
setGeneric(
		name="get_lsi.result<-",
		def=function(object,value){
			standardGeneric("get_lsi.result<-")
		})
setMethod("get_lsi.result",
		signature(object = "SingCellaR"),
		function (object){
			object@lsi.result
		})
setReplaceMethod(
		f="get_lsi.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,lsi.result=value)
		})
########
#' Get tSNE analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_tsne.result",
		def=function(object){
			standardGeneric("get_tsne.result")
		})
setGeneric(
		name="get_tsne.result<-",
		def=function(object,value){
			standardGeneric("get_tsne.result<-")
		})
setMethod("get_tsne.result",
		signature(object = "SingCellaR"),
		function (object){
			object@tsne.result
		})
setReplaceMethod(
		f="get_tsne.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,tsne.result=value)
		})

########
#' Get diffusion map analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_dfm.result",
		def=function(object){
			standardGeneric("get_dfm.result")
		})
setGeneric(
		name="get_dfm.result<-",
		def=function(object,value){
			standardGeneric("get_dfm.result<-")
		})
setMethod("get_dfm.result",
		signature(object = "SingCellaR"),
		function (object){
			object@diffusionmap.result
		})
setReplaceMethod(
		f="get_dfm.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,diffusionmap.result=value)
		})

########
#' Get UMAP analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_umap.result",
		def=function(object){
			standardGeneric("get_umap.result")
		})
setGeneric(
		name="get_umap.result<-",
		def=function(object,value){
			standardGeneric("get_umap.result<-")
		})
setMethod("get_umap.result",
		signature(object = "SingCellaR"),
		function (object){
			object@umap.result
		})
setReplaceMethod(
		f="get_umap.result",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,umap.result=value)
		})
########
########
#setGeneric(
#		name="get_tsne.result.high.contour",
#		def=function(object){
#			standardGeneric("get_tsne.result.high.contour")
#		})
#setGeneric(
#		name="get_tsne.result.high.contour<-",
#		def=function(object,value){
#			standardGeneric("get_tsne.result.high.contour<-")
#		})
#setMethod("get_tsne.result.high.contour",
#		signature(object = "SingCellaR"),
#		function (object){
#			object@tsne.result.high.contour
#		})
#setReplaceMethod(
#		f="get_tsne.result.high.contour",
#		signature="SingCellaR",
#		definition=function(object,value){
#			initialize(object,tsne.result.high.contour=value)
#		})
########
#' Get graph from KNN-graph analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_knn_graph.graph",
		def=function(object){
			standardGeneric("get_knn_graph.graph")
		})
setGeneric(
		name="get_knn_graph.graph<-",
		def=function(object,value){
			standardGeneric("get_knn_graph.graph<-")
		})
setMethod("get_knn_graph.graph",
		signature(object = "SingCellaR"),
		function (object){
			object@knn_graph.graph
		})
setReplaceMethod(
		f="get_knn_graph.graph",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,knn_graph.graph=value)
		})
########
#' Get graph-layout from KNN-graph analysis result
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_knn_graph.layout",
		def=function(object){
			standardGeneric("get_knn_graph.layout")
		})
setGeneric(
		name="get_knn_graph.layout<-",
		def=function(object,value){
			standardGeneric("get_knn_graph.layout<-")
		})
setMethod("get_knn_graph.layout",
		signature(object = "SingCellaR"),
		function (object){
			object@knn_graph.layout
		})
setReplaceMethod(
		f="get_knn_graph.layout",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,knn_graph.layout=value)
		})
########
########
setGeneric(
		name="get_knn_graph.kmeans.cluster",
		def=function(object){
			standardGeneric("get_knn_graph.kmeans.cluster")
		})
setGeneric(
		name="get_knn_graph.kmeans.cluster<-",
		def=function(object,value){
			standardGeneric("get_knn_graph.kmeans.cluster<-")
		})
setMethod("get_knn_graph.kmeans.cluster",
		signature(object = "SingCellaR"),
		function (object){
			object@knn_graph.kmeans.cluster
		})
setReplaceMethod(
		f="get_knn_graph.kmeans.cluster",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,knn_graph.kmeans.cluster=value)
		})
########
#' Get cluster information
#' @param  object The SingCellaR object.
#' @export 
########
setGeneric(
		name="get_clusters",
		def=function(object){
			standardGeneric("get_clusters")
		})
setGeneric(
		name="get_clusters<-",
		def=function(object,value){
			standardGeneric("get_clusters<-")
		})
setMethod("get_clusters",
		signature(object = "SingCellaR"),
		function (object){
			object@sc.clusters
		})
setReplaceMethod(
		f="get_clusters",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,sc.clusters=value)
		})

########
#' Get marker genes information
#' @param  object The SingCellaR object.
#' @export 
#' 
setGeneric(
		name="get_marker_genes",
		def=function(object){
			standardGeneric("get_marker_genes")
		})
setGeneric(
		name="get_marker_genes<-",
		def=function(object,value){
			standardGeneric("get_marker_genes<-")
		})
setMethod("get_marker_genes",
		signature(object = "SingCellaR"),
		function (object){
			object@marker.genes
		})
setReplaceMethod(
		f="get_marker_genes",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,marker.genes=value)
		})

########
setGeneric(
		name="get_igraph.graph",
		def=function(object){
			standardGeneric("get_igraph.graph")
		})
setGeneric(
		name="get_igraph.graph<-",
		def=function(object,value){
			standardGeneric("get_igraph.graph<-")
		})
setMethod("get_igraph.graph",
		signature(object = "SingCellaR"),
		function (object){
			object@igraph.graph
		})
setReplaceMethod(
		f="get_igraph.graph",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,igraph.graph=value)
		})
########
#' Get FA2 graph-layout
#' @param  object The SingCellaR object.
#' @export 
#' 
########
setGeneric(
		name="get_fa2_graph.layout",
		def=function(object){
			standardGeneric("get_fa2_graph.layout")
		})
setGeneric(
		name="get_fa2_graph.layout<-",
		def=function(object,value){
			standardGeneric("get_fa2_graph.layout<-")
		})
setMethod("get_fa2_graph.layout",
		signature(object = "SingCellaR"),
		function (object){
			object@fa2_graph.layout
		})
setReplaceMethod(
		f="get_fa2_graph.layout",
		signature="SingCellaR",
		definition=function(object,value){
			initialize(object,fa2_graph.layout=value)
		})
########
