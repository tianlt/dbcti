#define class ###############################
#' Define cell_trajectory class
#'
#' @slot raw_data data.frame.
#' @slot normalized_data data.frame.
#' @slot cell_annotation data.frame.
#' @slot nmf_result list.
#' @slot selected_feature list.
#' @slot specified_gene character.
#' @slot tsne_data list.
#' @slot contour_plot list.
#' @slot distribution_estimation list.
#' @slot point_possibility list.
#' @slot connect_cluster list.
#' @slot trajectory list.
#' @slot pseudotime list.
#' @slot trajectory_plot list.
#'
#' @return define a cell trajectory class
#' @export
#'
cell_trajectory <- setClass(
  Class = 'cell_trajectory',
  slots = c(
    raw_data = 'data.frame',
    normalized_data = 'data.frame',
    cell_annotation = 'data.frame',
    nmf_result = 'list',
    selected_feature = 'list', #nmf result and so on
    specified_gene = 'character',
    tsne_data = 'list',
    contour_plot = 'list',
    distribution_estimation = 'list',
    point_possibility = 'list',
    connect_cluster = 'list',
    trajectory = 'list',
    pseudotime = 'list',
    trajectory_plot = 'list'
  )
)


#create_object##############################
#' create cell_trajectory object
#'
#' @param data single-cell RNA-seq data frame with genes in row and cells in column
#' @param normalized indication of if data being normalized
#'
#' @return a cell_trajectory object
#' @export
#'
#' @examples
#' create_object(sample_data_df,  normalized = FALSE)
create_object <- function(data, normalized = FALSE){
  if (!is.data.frame(data)) stop('data must be in form of data frame')
  object <- new(Class = 'cell_trajectory')
  if (normalized == FALSE) object@raw_data <- data else object@normalized_data <- data
  return(object)
}


#normalize ######################
#' Normalize the raw counts data
#'
#' @param object a cell_trajectory object
#' @param gene_cri criteria for the sum of counts of any gene in all cells to be greater than
#' @param cell_cri criteria for cells containing genes counts greater than
#' @param scale_factor a scale factor for normalization
#'
#' @return values in the normalized_data slot of cell_trajectory object
#' @export
#'
#' @examples
#' normalize(sample_data, gene_cri = 0, cell_cri = 0, scale_factor = 10000)
normalize <- function(object, gene_cri = 0, cell_cri = 0, scale_factor = 10000){
  data <- object@raw_data
  data <- data[apply(data, 1, sum) >= gene_cri, apply(data, 2, sum) >= cell_cri]
  normalized_data <- as.data.frame(t(t(data)/apply(data, 2, sum)))*scale_factor

  object@normalized_data <- normalized_data
  return(object)
}



#filter_data ######################
#' filter data but not perform normalization
#'
#' @param object a cell_trajectory object
#' @param gene_cri gene_cri criteria for the sum of counts of any gene in all cells to be greater than
#' @param cell_cri criteria for cells containing genes counts greater than
#' @param use_normalized_data if or not use normalized data as input
#'
#' @return values in the normalized_data slot of cell_trajectory object
#' @export
#'
#' @examples
#' filter_data(sample_data, gene_cri = 0, cell_cri = 0, use_normalized_data = TRUE)
filter_data <- function(object, gene_cri = 0, cell_cri = 0, use_normalized_data = FALSE){
  if (use_normalized_data == FALSE) data <- object@raw_data else data <- object@normalized_data
  data <- data[apply(data, 1, sum) >= gene_cri, apply(data, 2, sum) >= cell_cri]
  object@normalized_data <- data
  return(object)
}


#feature selection function###################

##feature selection  # check if provided marker gene is duplicated
#NMF step 1 ##################################
#gene matrix as input and output w and h matrix
#' Select features via NMF step 1
#'
#' @param object a cell_trajectory object
#' @param use_normalized_data if or not use normalized data as input
#' @param rank rank being remained for the decomposed matrix
#' @param thread thread number used for parallel computation
#'
#' @return values in the nmf_result slot of cell_trajectory object
#' @export
#'
#' @importFrom NMF nmf
#' @examples
#' data_nmf(sample_data, use_normalized_data = TRUE, thread = 1)
data_nmf<-function(object, use_normalized_data = TRUE, rank=10, thread=1){
  if (use_normalized_data == FALSE) data_truncted<-object@raw_data[apply(object@raw_data, 1, function(x) sum(x)!=0),] else data_truncted<-object@normalized_data[apply(object@normalized_data, 1, function(x) sum(x)!=0),]
  data_nmf<-NMF::nmf(data_truncted,rank=rank,.opt=paste0('vP',thread), seed=12345)
  object@nmf_result<-list(W=data_nmf@fit@W, H=data_nmf@fit@H, rank=rank)
  return(object)
}




#nmf step 2 #if featrue has dumplicated? ################################
#KNN NMF W matrix as input output as selected gene name
#' Select features via NMF step 2
#'
#' @param object a cell_trajectory object
#' @param k number of nearest neighbour
#' @param feature the specifeid features
#'
#' @return values in the selected_feature slot of cell_trajectory object
#' @export
#' @importFrom FNN get.knnx
#' @examples
#' feature_selection_knn(sample_data,feature = c('1','2'), k = 5)
feature_selection_knn<-function(object, k=10, feature){
  #feature ratio
  feature_remained<-feature[feature %in% rownames(object@nmf_result$W)]
  ratio<-length(feature_remained)/length(feature)


  #knn
  data_knn<-FNN::get.knnx(object@nmf_result$W,object@nmf_result$W[feature_remained,],k=k)
  knn_vec_matrix<-as.vector(data_knn[["nn.index"]])
  selected_gene<-rownames(object@nmf_result$W[knn_vec_matrix[!duplicated(knn_vec_matrix)],])
  gene_length<-length(selected_gene)
  object@selected_feature<-list(selected_gene=selected_gene,gene_length=gene_length, mapped_gene_ratio=ratio, method='knn', k=k, gene_for_select = feature)
  return(object)
}




#select_cor_feature########################
#' Select features based on correlation
#'
#' @param object a cell_trajectory object
#' @param feature the specifeid features
#' @param use_normalized_data if or not use normalized data as input
#' @param k number of highest correlation genes for each target
#'
#' @return values in the selected_feature slot of cell_trajectory object
#' @export
#'
#' @examples
#' select_cor_feature(sample_data, c('1','2'), use_normalized_data = FALSE, k = 5)
select_cor_feature <- function(object, feature, use_normalized_data = FALSE, k){
  if (use_normalized_data == FALSE) data = object@raw_data else data = object@normalized_data
  feature_remained<-feature[feature %in% rownames(data)]
  ratio<-length(feature_remained)/length(feature)

  cor <- as.data.frame(cor(t(data[feature_remained,]), t(data)))
  cor_res <- as.vector(t(apply(cor, 1, function(x) names(cor)[sort(head(order(x, decreasing = TRUE), k))])))
  cor_res <- cor_res[!duplicated(cor_res)]
  gene_length<-length(cor_res)
  object@selected_feature <- list(selected_gene = cor_res, gene_length=gene_length, mapped_gene_ratio=ratio, method = 'cor', k=k, gene_for_select = feature)
  return(object)
}


#select_var_feature########################
#' Select features based on variation
#'
#' @param object a cell_trajectory object
#' @param use_normalized_data if or not use normalized data as input
#' @param n number of features with the highest variation
#'
#' @return values in the selected_feature slot of cell_trajectory object
#' @export
#'
#' @examples
#' select_var_feature(sample_data, use_normalized_data = FALSE, n = 10)
select_var_feature <- function(object, use_normalized_data = FALSE, n = 2000){
  if (use_normalized_data == FALSE) data = object@raw_data else data = object@normalized_data
  var <- apply(data, 1, var)
  object@selected_feature <- list(selected_gene = names(sort(var, decreasing = TRUE)[1:n]), method = 'var')
  return(object)
}


#tsneplot #################################
#input object, perplexity, max_iter, use_normalized_data...
#if specified gene = true, specified gene will be used for tsne plot
#' Tsne plot for cell_trajectory object
#'
#' @param object a cell_trajectory object
#' @param perplexity perplexity value for tsne
#' @param max_iter max iteration number
#' @param specified_gene feature names if specified
#' @param use_normalized_data if or not use normalized data as input
#' @param pca if pca will be performed
#' @param check_duplicates if or not duplicates are checked at tsne step
#' @param title title of plot, default as none
#' @param file file the plot to be written to, defalut as empity
#'
#' @return values in the tsne_data slot of cell_trajectory object
#' @export
#' @importFrom Rtsne Rtsne
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @examples
#' tsneplot(sample_data, use_normalized_data = TRUE, perplexity = 5)
tsneplot<-function(object, perplexity=30, max_iter = 500, specified_gene = FALSE, use_normalized_data = TRUE, pca=TRUE, check_duplicates=FALSE, title='', file='', initial_dims = 50){
  set.seed(12345)
  if (specified_gene == TRUE) {
    if (use_normalized_data == FALSE) data_feature <- scale(object@raw_data)[object@specified_gene,] else data_feature <- scale(object@normalized_data)[object@specified_gene,]
  }
  else {
    if (use_normalized_data == FALSE) data_feature <- scale(object@raw_data)[object@selected_feature$selected_gene,] else data_feature <- scale(object@normalized_data)[object@selected_feature$selected_gene,]
  }
  tsne <- Rtsne::Rtsne(t(data_feature), dims = 2, perplexity=perplexity, verbose=FALSE, max_iter = max_iter, pca=pca, check_duplicates = check_duplicates, initial_dims = initial_dims)
  if (nchar(file)>0) {
    pdf(file = file)
    plot(tsne$Y, main=title)
    dev.off()
  }
  plot(tsne$Y, main=title)
  object@tsne_data <- tsne
  return(object)
}


#trajectory inference function ####################
#contour plot input matrix cells in 2d space(after tsne) output gorup estimation
#' Title
#'
#' @param object a cell_trajectory object
#' @param n number of grid points in each direction
#'
#' @return values in the contour_plot slot of cell_trajectory object
#' @export
#' @importFrom MASS kde2d
#' @importFrom mclust densityMclust
#' @examples
#' contour_plot(sample_data)
contour_plot<-function(object, n=1000){
  data<-as.data.frame(object@tsne_data$Y)
  data_kde<-MASS::kde2d(data[,1],data[,2], n=1000)
  dense_group<-mclust::densityMclust(data = data)
  contour(data_kde)
  text(data,labels=dense_group[["classification"]], col=dense_group[["classification"]])
  object@contour_plot <- list(dense_group)
  return(object)
}



#distribution_estimation
#density estimation and cluster input marix of cells in 2d dimension, classification index by mclust, group number to cell state list(vertex), ndraw draw from distribution, expansion, expand the estimated distribution
#output fitted distribution parameter list and drawed distribution list
#i=state j=cluster in i
#fit the distribution and then build distribution based on fitted parameter

#' Estimate distribution
#'
#' @param object a cell_trajectory object
#' @param ndraw number of draws
#' @param expansion eapansion parameter for fitted distribution
#' @param ... cluster index
#'
#' @return values in the distribution_estimation slot of cell_trajectory object
#' @export
#' @importFrom MASS mvrnorm
#' @examples
#' distribution_estimation(sample_data, ndraw = 50, expansion = 1.5, ... = 1,2,3)
distribution_estimation<-function(object, ndraw = 1000, expansion = 1.5, ...){
  distribution_fit_list<-distribution_draw_list<-list()
  data = as.data.frame(object@tsne_data$Y)

  #build new cluster index

  cluster_index_new <- rep(0, length(object@contour_plot[[1]]$classification))
  names(cluster_index_new) <- 1:length(object@contour_plot[[1]]$classification)
  for (i in 1:length(list(...))) {
    n<-length(list(...)[[i]])
    #situation where cluster in state i >=2
    if (n>=2) {
      cell_selector<-object@contour_plot[[1]][["classification"]]==list(...)[[i]][1]
      for (j in 2:n) {
        cell_selector<-cell_selector | object@contour_plot[[1]][["classification"]]==list(...)[[i]][j]
      }
      assign(paste('d',i, sep = '_'), fr(as.matrix(data[cell_selector,])))
      #for new cluster index
      cluster_index_new[cell_selector] <- rep(i, sum(cell_selector))
    }
    #situation where cluster in state i >=2
    if (n==1) {
      assign(paste('d',i, sep = '_'), fr(as.matrix(data[object@contour_plot[[1]][["classification"]]==list(...)[[i]],])))
      cluster_index_new[object@contour_plot[[1]][["classification"]]==list(...)[[i]]] <- rep(i, sum(object@contour_plot[[1]][["classification"]]==list(...)[[i]]))
    }
    #cluster_index_new[cluster_index_new==2] <- rep(i, sum(cluster_index_new==2))


  }
  #examine the mean of mixture gaussian distribution

  for (k in 1:length(list(...))) {
    d_name<-paste('d', k, sep = '_')
    distribution_fit_list[[d_name]]<-eval(parse(text=d_name))

    #draw from distribution
    ds_name<-paste('ds', k, sep = '_')
    mean_k<-eval(parse(text=d_name))$mu
    var_k<-eval(parse(text=d_name))$S*expansion

    distribution_draw_list[[ds_name]]<-MASS::mvrnorm(n = ndraw, mu = mean_k, Sigma = var_k)
  }
  #new cluster index


  result<-list(distribution_fit_list=distribution_fit_list, distribution_draw_list=distribution_draw_list, cluster_index=cluster_index_new)
  object@distribution_estimation <- result
  return(object)
}


#test_distribution_estimation###############
#test_hescmt <- distribution_estimation(test_hescmt, ndraw = 1000, expansion = 1.5, ... = 1,2,3)
#test_distribution_estimation<-distribution_estimation(data = dat, group = test_contour_plot, ndraw = 1000, ... = c(1,2), c(3,4), 5)

#generate vague point possibility
#input r as radius, ds_list, 2d cell data
#ouput counted_distribution and counted_possibility, na value in counted_possibility means no count for all distribution. and distant point that not belong to none distribution
#' Calculate point possibility
#'
#' @param object a cell_trajectory object
#' @param r distance from the data point for the counting area
#'
#' @return values in the point_possibility slot of cell_trajectory object
#' @export
#'
#' @examples
#' point_possibility(sample_data , r = 2)
point_possibility<-function(object, r){
  ds_list = object@distribution_estimation$distribution_draw_list
  data = as.data.frame(object@tsne_data$Y)

  data_range<-data.frame(cbind(data-r, data+r))
  colnames(data_range)<-c('x-', 'y-', 'x+', 'y+')
  count_1=0 #all points counted
  counted_point<-data.frame()
  counted_distribution<-data.frame(matrix(nrow = nrow(data), ncol = 1))
  for (k in 1:length(ds_list)) {
    #distribution<-eval(parse(text = paste('ds', k, sep = '_')))
    distribution<-ds_list[[paste('ds', k, sep = '_')]]

    for (i in 1:nrow(data_range)) {
      count_2=0 #point j counted
      for (j in 1:nrow(distribution)) {
        if (data_range[i,1]<=distribution[j,1] && data_range[i,3]>=distribution[j,1] && data_range[i,2] <= distribution[j,2]&& data_range[i,4] >= distribution[j,2]) {

          count_1<-count_1+1
          count_2<-count_2+1

        }

      }
      counted_point<-rbind(counted_point, count_2)
    }
    counted_distribution<-cbind(counted_distribution,counted_point)
    counted_point<-data.frame()
  }
  counted_distribution<-counted_distribution[,-1]
  colnames(counted_distribution)<-paste('ds', 1:length(ds_list), sep = '_')
  #counted_possibility matrix
  counted_possibility<-counted_distribution/rowSums(counted_distribution)

  #distant point ratio, ratio of point belongs to neither distribution
  distant_point_ratio<-sum(is.na(counted_possibility[,1]))/nrow(counted_possibility)

  result<-list(counted_possibility=counted_possibility, counted_distribution=counted_distribution,
               total_counts=count_1, distant_point_ratio=distant_point_ratio)
  object@point_possibility <- result
  return(object)
}


#find vague point input counted_possibility, sum criteria and difference criteria, vageu point criteria
#sum_cri the sum of both distribution exceed some value, diff_cri the difference of both distribution exceed some value, vague_cri
#output ratio list, cluster_connection matrix

#' Connect cluster
#'
#' @param object a cell_trajectory object
#' @param sum_cri criteria set for sum of a data point belonging to two distribution
#' @param diff_cri criteria set for difference between a data point belonging to two distribution
#' @param vague_cri criteria set for the ratio of vague point in two distribution
#'
#' @return values in the connect_cluster slot of cell_trajectory object
#' @export
#'
#' @examples
#' connect_cluster(sample_data)
connect_cluster<-function(object, sum_cri = 0.8, diff_cri = 0.5, vague_cri = 0.01){

  counted_possibility = object@point_possibility$counted_possibility
  n=ncol(counted_possibility)
  vague_ratio<-list()
  cluster_connection<-matrix(0, ncol = n, nrow = n)

  counted_possibility<-counted_possibility[!is.na(counted_possibility)[,1], ]
  #i=first distribution for comparison and j= the second
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {

      ratio<-sum((counted_possibility[,i]+counted_possibility[,j])>=sum_cri &
                   abs(counted_possibility[,i]-counted_possibility[,j])<=diff_cri)/
        sum((counted_possibility[,i]+counted_possibility[,j])>=sum_cri)

      if (is.na(ratio)) ratio = 0

      list_name<-paste('ratio',i,j, sep = '_')
      vague_ratio[[list_name]]<-ratio

      #determine cluster_trajectory matrix
      if (ratio>=vague_cri) cluster_connection[i,j]<-cluster_connection[j,i]<-1

    }
  }
  rownames(cluster_connection)<-colnames(cluster_connection)<-1:nrow(cluster_connection)
  result<-list(cluster_connection=cluster_connection, vague_ratio=vague_ratio)
  object@connect_cluster <- result
  return(object)
}


#fit principle curve iteratively
#input connection matrix(n*n) from above step, cluster_df collumn as x, y and cluster index, and corresponding cell name, iter_n
#i.e. test<-infer_trajectory(connection_matrix = connection_matrix, tsne_df = dat, cluster_index = cluster_index)
#' Title
#'
#' @param object a cell_trajectory object
#' @param iter_n iteration number
#'
#' @return values in the trajectory slot of cell_trajectory object
#' @export
#' @importFrom princurve principal_curve
#' @examples
#' infer_trajectory(sample_data, iter_n = 50)
infer_trajectory<-function(object, iter_n = 50){
  connection_matrix = object@connect_cluster$cluster_connection
  tsne_df = as.data.frame(object@tsne_data$Y)
  cluster_index = object@distribution_estimation$cluster_index


  #start####################
  if (nrow(connection_matrix)<2) stop('not enough clusters')
  n_cluster=length(unique(cluster_index))
  lines_list<-result_list<-list()
  for (i in 1:n_cluster){
    cluster_name_i<-rownames(connection_matrix)[i]
    #build result list for iteration
    result_list[[cluster_name_i]]<-tsne_df[cluster_index==cluster_name_i, ]

  }
  ###########for cluster wiht no other connection with other clusters################
  isolated_cluster<-c()
  for (i in 1:n_cluster) {
    cluster_name_i<-rownames(connection_matrix)[i]
    if (sum(connection_matrix[, i]) == 0) {
      fit_data<-as.matrix(result_list[[cluster_name_i]])
      #fit_data<-as.matrix(tsne_df[(cluster_index==cluster_name_i | cluster_index==cluster_name_j),])
      fit_i<-princurve::principal_curve(fit_data)
      resulted_i<-as.data.frame(fit_i$s)

      colnames(resulted_i)<-c('x', 'y')
      result_list[[cluster_name_i]]<-resulted_i
      lines_list_name<-paste('fit', cluster_name_i, sep = '_')
      lines_list[[lines_list_name]] <- fit_i
      isolated_cluster<-append(isolated_cluster, cluster_name_i)
    }
  }

  for (k in 1:iter_n) {
    for (i in 1:n_cluster) {
      cluster_name_i<-rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }

      #dataframe for summarise all regression projected point for each cluster

      assign(paste('fc', cluster_name_i, 'x', sep = '_'), data.frame(matrix(0, nrow = sum(cluster_index==cluster_name_i), ncol = 1)))
      assign(paste('fc', cluster_name_i, 'y', sep = '_'), data.frame(matrix(0, nrow = sum(cluster_index==cluster_name_i), ncol = 1)))

      #index for each cluster
      assign(paste('cluster_index', cluster_name_i, sep = '_'), names(cluster_index)[cluster_index==cluster_name_i])
    }

    #fit principle curve based on connection matrix###################
    for (i in 1:(nrow(connection_matrix)-1)) {
      for (j in (i+1):nrow(connection_matrix)) {
        if (connection_matrix[i,j]==1) {
          cluster_name_i<-rownames(connection_matrix)[i]
          cluster_name_j<-rownames(connection_matrix)[j]
          ################if cluster is the same as above (2 cluster may combine in to1)
          fit_data<-as.matrix(rbind(result_list[[cluster_name_i]], result_list[[cluster_name_j]]))
          #fit_data<-as.matrix(tsne_df[(cluster_index==cluster_name_i | cluster_index==cluster_name_j),])
          assign(paste('fit', cluster_name_i, cluster_name_j, sep = '_'), princurve::principal_curve(fit_data))
          fitted<-eval(parse(text = paste('fit', cluster_name_i, cluster_name_j, sep = '_')))$s

          fitted_i_x<-fitted[eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_'))),][,1]
          fitted_i_y<-fitted[eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_'))),][,2]

          fitted_j_x<-fitted[eval(parse(text = paste('cluster_index', cluster_name_j, sep = '_'))),][,1]
          fitted_j_y<-fitted[eval(parse(text = paste('cluster_index', cluster_name_j, sep = '_'))),][,2]

          #x and y for each principle curve for each cluster
          assign(paste('fc', cluster_name_i, 'x', sep = '_'), cbind(eval(parse(text = paste('fc', cluster_name_i, 'x', sep = '_'))), fitted_i_x))
          assign(paste('fc', cluster_name_i, 'y', sep = '_'), cbind(eval(parse(text = paste('fc', cluster_name_i, 'y', sep = '_'))), fitted_i_y))

          assign(paste('fc', cluster_name_j, 'x', sep = '_'), cbind(eval(parse(text = paste('fc', cluster_name_j, 'x', sep = '_'))), fitted_j_x))
          assign(paste('fc', cluster_name_j, 'y', sep = '_'), cbind(eval(parse(text = paste('fc', cluster_name_j, 'y', sep = '_'))), fitted_j_y))
        }
      }
    }

    #average result###################
    for (i in 1:n_cluster) {
      cluster_name_i<-rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      x_df<-eval(parse(text = paste('fc', cluster_name_i, 'x', sep = '_')))
      y_df<-eval(parse(text = paste('fc', cluster_name_i, 'y', sep = '_')))
      assign(paste('fc', cluster_name_i, 'x', sep = '_'), as.data.frame(x_df[,-1], row.names = rownames(x_df)))
      assign(paste('fc', cluster_name_i, 'y', sep = '_'), as.data.frame(y_df[,-1], row.names = rownames(y_df)))
      data_x<-eval(parse(text = paste('fc', cluster_name_i, 'x', sep = '_')))
      data_y<-eval(parse(text = paste('fc', cluster_name_i, 'y', sep = '_')))

      assign(paste('loc', cluster_name_i, 'x', sep = '_'), apply(as.data.frame(data_x), 1, mean))
      assign(paste('loc', cluster_name_i, 'y', sep = '_'), apply(as.data.frame(data_y), 1, mean))
      assign(paste('loc', cluster_name_i, sep = '_'), data.frame(x=eval(parse(text = paste('loc', cluster_name_i, 'x', sep = '_'))), y=eval(parse(text = paste('loc', cluster_name_i, 'y', sep = '_')))))
      result_list[[cluster_name_i]]<-eval(parse(text = paste('loc', cluster_name_i, sep = '_')))
    }
  }

  #for lines data
  for (i in 1:(nrow(connection_matrix)-1)) {
    for (j in (i+1):nrow(connection_matrix)) {
      cluster_name_i<-rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      cluster_name_j<-rownames(connection_matrix)[j]
      if (connection_matrix[i,j]==1) {
        fitted<-eval(parse(text = paste('fit', cluster_name_i, cluster_name_j, sep = '_')))
        lines_list_name<-paste('fit', cluster_name_i, cluster_name_j, sep = '_')
        lines_list[[lines_list_name]]<-fitted
      }
    }
  }
  result<-list(result_list=result_list, lines_list=lines_list)
  object@trajectory <- result
  return(object)
}



#connection_matrix, trajectory, start_state_name(character), separate==FALSE, cluster_index
#' calculate pseudotime for each cell
#'
#' @param object a cell_trajectory object
#' @param start_state_name index for the state that starts
#'
#' @return values in the pseudotime slot of cell_trajectory object
#' @export
#' @import igraph
#' @importFrom stringr str_detect
#' @importFrom stringr str_split
#' @examples
#' calculate_pseudotime(sample_data, start_state_name = c('1','2'))
calculate_pseudotime<-function(object, start_state_name){
  connection_matrix = object@connect_cluster$cluster_connection
  trajectory = object@trajectory
  cluster_index = object@distribution_estimation$cluster_index

  #build graph
  g<-igraph::graph_from_adjacency_matrix(connection_matrix, mode = 'undirected')
  #if separated
  if (igraph::components(g)$no == 1) separated = FALSE else separated = TRUE
  #build cluster index for each cluster
  for(i in 1:nrow(connection_matrix)){
    cluster_name_i<-rownames(connection_matrix)[i]
    assign(paste('cluster_index', cluster_name_i, sep = '_'), names(cluster_index)[cluster_index==cluster_name_i])
  }
  #separated == false ###################
  if (separated==FALSE) {
    #calculate distance of each cluster i
    for(i in 1:nrow(connection_matrix)){
      d_list<-list()
      cluster_name_i<-rownames(connection_matrix)[i]
      fit_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
      part<-names(trajectory$lines_list)[fit_index]
      #for each part j of cluster i calculate the mean of j
      for (j in 1:length(part)) {
        part_index_j<-part[j]
        tra_lam<-trajectory$lines_list[[part_index_j]]$lambda
        d_part_j<-tra_lam[names(tra_lam) %in% eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))]
        d_range_j<-max(d_part_j)-min(d_part_j)
        d_list<-append(d_list, d_range_j)
      }
      assign(paste('d', cluster_name_i, sep = '_'), mean(unlist(d_list)))
    }

    #assign psdeotime
    psdeotime_list<-list()
    for(i in 1:nrow(connection_matrix)){
      cluster_name_i<-rownames(connection_matrix)[i]
      path<-names(unlist(igraph::shortest_paths(g, cluster_name_i, start_state_name)[[1]]))


      cl_point<-eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))
      #find the start point of the cluster
      #situation where point is not in the cluster of start
      if (length(path)>=2) {
        if (as.numeric(path[1]) < as.numeric(path[2])) curve_name<-paste('fit', path[1], path[2], sep = '_') else curve_name<-paste('fit', path[2], path[1], sep = '_')

        ordered<-trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ][rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ]) %in% cl_point, ]
        name_first<-rownames(ordered)[1]
        name_last<-rownames(ordered)[nrow(ordered)]

        match_first<-match(name_first, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        match_last<-match(name_last, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', path[2], sep = '_'))), rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,])))
        if (avg_other_cl<=match_first) start_point<-name_first else start_point<-name_last

        #distance
        d_within_cluster<-abs(trajectory$lines_list[[curve_name]]$lambda[start_point]-trajectory$lines_list[[curve_name]]$lambda[cl_point])
        d_total<-d_within_cluster
        for (j in path[-1]) d_total <- d_total + eval(parse(text = paste('d', j, sep = '_')))

      }

      #situation where cluster is state 0
      if (length(path)==1){
        part_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
        part<-names(trajectory$lines_list)[part_index]
        path_1_df<-data.frame(matrix(nrow = length(cl_point)), row.names = cl_point)
        n=0
        for (j in part) {
          n=n+1
          ordered<-trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ][rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ]) %in% cl_point, ]
          name_first<-rownames(ordered)[1]
          name_last<-rownames(ordered)[nrow(ordered)]

          match_first<-match(name_first, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))
          match_last<-match(name_last, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))

          #other cluster
          if (cluster_name_i == unlist(stringr::str_split(j, '_'))[2]) other_cluster<-unlist(stringr::str_split(j, '_'))[3] else other_cluster<-unlist(stringr::str_split(j, '_'))[2]


          avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', other_cluster, sep = '_'))), rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,])))
          if (avg_other_cl<=match_first) {
            start_point<-name_first

            #all distance
            d_within_cluster<-eval(parse(text = paste('d', cluster_name_i, sep = '_')))-abs(trajectory$lines_list[[j]]$lambda[start_point]-trajectory$lines_list[[j]]$lambda[cl_point])

            path_1_df[,n]<-d_within_cluster
          } else {
            start_point<-name_first

            #all distance
            d_within_cluster<-trajectory$lines_list[[j]]$lambda[cl_point]

            path_1_df[,n]<-d_within_cluster
          }
        }


        #average distance
        d_total<-apply(path_1_df, 1, mean)
      }
      psdeotime_list[[cluster_name_i]]<-d_total
    }
  }
  #separated == true ####################
  if (separated == TRUE) {
    #if the start_state_name meet requirement
    g_comp<-igraph::components(g)
    print(g_comp)
    if (length(start_state_name)<=1) stop('start_state_name for separated graph should be a at least length of 2 vertex')


    start_cluster_list<-c()
    for (i in start_state_name) {
      separated_comp<-g_comp$membership[match(i, names(igraph::V(g)))]
      start_cluster_list<-append(start_cluster_list, separated_comp)
    }
    if (sum(duplicated(start_cluster_list)) >= 1) stop('multiple start state name from same graph part')


    #calculate distance of each cluster i
    for(i in 1:nrow(connection_matrix)){
      d_list<-list()
      cluster_name_i<-rownames(connection_matrix)[i]
      fit_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
      part<-names(trajectory$lines_list)[fit_index]
      #for each part j of cluster i calculate the mean of j
      for (j in 1:length(part)) {
        part_index_j<-part[j]
        tra_lam<-trajectory$lines_list[[part_index_j]]$lambda
        d_part_j<-tra_lam[names(tra_lam) %in% eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))]
        d_range_j<-max(d_part_j)-min(d_part_j)
        d_list<-append(d_list, d_range_j)
      }
      assign(paste('d', cluster_name_i, sep = '_'), mean(unlist(d_list)))
    }

    #isolated part situation##########################
    psdeotime_list<-list()
    isolated_cluster<-c()
    for (i in 1:nrow(connection_matrix)) {
      cluster_name_i<-rownames(connection_matrix)[i]
      if (sum(connection_matrix[, i]) == 0) {
        psdeotime_list[[cluster_name_i]]<-trajectory$lines_list[[paste('fit', cluster_name_i, sep = '_')]]$lambda
        isolated_cluster<-append(isolated_cluster, cluster_name_i)
      }
    }
    #
    #assign psdeotime

    for(i in 1:nrow(connection_matrix)){
      cluster_name_i<-rownames(connection_matrix)[i]
      if (cluster_name_i %in% isolated_cluster) {
        next
      }
      #build cluster and their start state relationship
      separated_comp_i<-g_comp$membership[match(cluster_name_i, names(igraph::V(g)))]
      for (k in start_state_name) {
        separated_comp_k_start<-g_comp$membership[match(k, names(igraph::V(g)))]
        if (separated_comp_k_start == separated_comp_i) break
      }

      path<-names(unlist(igraph::shortest_paths(g, cluster_name_i, names(separated_comp_k_start))[[1]]))
      cl_point<-eval(parse(text = paste('cluster_index', cluster_name_i, sep = '_')))
      #line################
      #find the start point of the cluster
      #situation where point is not in the cluster of start
      if (length(path)>=2) {
        if (as.numeric(path[1]) < as.numeric(path[2])) curve_name<-paste('fit', path[1], path[2], sep = '_') else curve_name<-paste('fit', path[2], path[1], sep = '_')

        ordered<-trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ][rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord, ]) %in% cl_point, ]
        name_first<-rownames(ordered)[1]
        name_last<-rownames(ordered)[nrow(ordered)]

        match_first<-match(name_first, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        match_last<-match(name_last, rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,]))
        avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', path[2], sep = '_'))), rownames(trajectory$lines_list[[curve_name]]$s[trajectory$lines_list[[curve_name]]$ord,])))
        if (avg_other_cl<=match_first) start_point<-name_first else start_point<-name_last

        #distance
        d_within_cluster<-abs(trajectory$lines_list[[curve_name]]$lambda[start_point]-trajectory$lines_list[[curve_name]]$lambda[cl_point])
        d_total<-d_within_cluster
        for (j in path[-1]) d_total <- d_total + eval(parse(text = paste('d', j, sep = '_')))

      }

      #situation where cluster is state 0
      if (length(path)==1){
        part_index<-stringr::str_detect(names(trajectory$lines_list), cluster_name_i)
        part<-names(trajectory$lines_list)[part_index]
        path_1_df<-data.frame(matrix(nrow = length(cl_point)), row.names = cl_point)
        n=0
        for (j in part) {
          n=n+1
          ordered<-trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ][rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord, ]) %in% cl_point, ]
          name_first<-rownames(ordered)[1]
          name_last<-rownames(ordered)[nrow(ordered)]

          match_first<-match(name_first, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))
          match_last<-match(name_last, rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,]))

          #other cluster
          if (cluster_name_i == unlist(stringr::str_split(j, '_'))[2]) other_cluster<-unlist(stringr::str_split(j, '_'))[3] else other_cluster<-unlist(stringr::str_split(j, '_'))[2]


          avg_other_cl<-mean(match(eval(parse(text = paste('cluster_index', other_cluster, sep = '_'))), rownames(trajectory$lines_list[[j]]$s[trajectory$lines_list[[j]]$ord,])))
          if (avg_other_cl<=match_first) {
            start_point<-name_first

            #all distance
            d_within_cluster<-eval(parse(text = paste('d', cluster_name_i, sep = '_')))-abs(trajectory$lines_list[[j]]$lambda[start_point]-trajectory$lines_list[[j]]$lambda[cl_point])

            path_1_df[,n]<-d_within_cluster
          } else {
            start_point<-name_first

            #all distance
            d_within_cluster<-trajectory$lines_list[[j]]$lambda[cl_point]

            path_1_df[,n]<-d_within_cluster
          }
        }
        #average distance
        d_total<-apply(path_1_df, 1, mean)
      }
      psdeotime_list[[cluster_name_i]]<-d_total
    }
  }
  object@pseudotime <- psdeotime_list
  return(object)
}

#fr##############
#' Title
#'
#' @param x refer to prada
#' @param y refer to prada
#' @param scalefac refer to prada
#' @param method refer to prada
#' @param noise refer to prada
#' @param gateName refer to prada
#'
#' @return fitted distribution parameters
#'
#' @import prada
#' @importFrom rrcov CovMcd
#' @importFrom MASS cov.rob

fr<-function (x, y = NA, scalefac = 1, method = "covMcd", noise,
              gateName = "fitNorm")
{
  if (is(x, "cytoFrame"))
    x <- exprs(x)[, 1:2]
  if (!(is.matrix(x) && ncol(x) == 2)) {
    if (!length(x) == length(y) || !is.numeric(x) || !is.numeric(y))
      stop("'x' and 'y' must be numeric vectors of equal length")
    x <- cbind(x, y)
  }
  xorig <- x
  if (!missing(noise)) {
    if (is.logical(noise))
      noise <- which(noise)
    if (!is.numeric(noise) || length(noise) > nrow(x) ||
        length(noise) == 0)
      stop("'noise' should be an index or logical vector not longer than x")
    x <- x[-noise, , drop = FALSE]
  }
  if (!is.numeric(scalefac))
    stop("'scalefac' must be numeric")
  cov <- switch(method, covMcd = {
    nmax <- 50000
    if (nrow(x) > nmax) tmp <- CovMcd(x[sample(nrow(x), nmax),
                                        ]) else tmp <- CovMcd(x)
                                        list(center = tmp@center, cov = tmp@cov)
  }, cov.rob = {
    cov.rob(x)
  }, stop("'method' must be one of 'covMcd' or 'cov.rob'"))
  mu <- cov$center
  S <- cov$cov
  Sinv <- solve(S)
  w <- rbind(xorig[, 1], xorig[, 2]) - mu
  z <- Sinv %*% w
  p <- exp(-0.5 * (z[1, ] * w[1, ] + z[2, ] * w[2, ]))
  sel <- p > exp(-0.5 * scalefac^2)
  gfun <- function(x = x, cov, scalefac) {
    mu <- cov$center
    S <- cov$cov
    Sinv <- solve(S)
    w <- rbind(x[, 1], x[, 2]) - mu
    z <- Sinv %*% w
    p <- exp(-0.5 * (z[1, ] * w[1, ] + z[2, ] * w[2, ]))
    return(p > exp(-0.5 * scalefac^2))
  }
  cn <- colnames(x)
  if (is.null(cn))
    colnames(xorig) <- c("x", "y")
  gate <- new("gate", name = gateName, gateFun = function(x) gfun(x = x,
                                                                  cov = cov, scalefac = scalefac), colnames = colnames(xorig),
              logic = "&", type = "fitNorm")
  return(invisible(list(mu = mu, S = S, p = p, sel = sel, scalefac = scalefac,
                        data = xorig, gate = gate)))
}

#plot_trejactory##########################
######################
#' construct trajectory plot for cell_trajectory object
#'
#' @param object a cell_trajectory object
#' @param plot_title title of the plot
#' @param width jitter width
#' @param height jitter height
#'
#' @return values in the trajectory_plot slot of cell_trajectory object
#' @export
#' @import igraph
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @examples
#' plot_trajectory(sample_data)
plot_trajectory<-function(object, plot_title = '', width = 0.2, height = 0.2){

  pseudotime = object@pseudotime
  trajectory = object@trajectory
  connection_matrix = object@connect_cluster$cluster_connection

  g<-igraph::graph_from_adjacency_matrix(connection_matrix, mode = 'undirected')
  g_comp<-igraph::components(g)
  n <- tail(sort(table(g_comp$membership)), 1) + 1
  ordered_pseudotime<-pseudotime[[1]]
  ordered_trajectory<-trajectory$result_list[[1]]
  for (i in pseudotime[-1]) ordered_pseudotime<-rbind(as.matrix(i), as.matrix(ordered_pseudotime))
  for (i in trajectory$result_list[-1]) ordered_trajectory<-rbind(as.matrix(i), as.matrix(ordered_trajectory))

  ordered_pseudotime<-ordered_pseudotime[order(as.numeric(rownames(ordered_pseudotime))),]
  ordered_trajectory<-ordered_trajectory[order(as.numeric(rownames(ordered_trajectory))),]

  gg_data<-as.data.frame(ordered_trajectory)

  g<-ggplot()

  for (i in trajectory[["lines_list"]]) {
    line <- as.data.frame(i$s[i$ord,])
    colnames(line) <- c('x', 'y')
    g <- g + geom_path(data = line, aes_string('x','y'), size = 1.5)
  }

  g <- g + geom_point(data = gg_data, aes_string('x','y', fill = ordered_pseudotime), color = 'black', size=8, pch =21, stroke = 0.5, position = position_jitter(width = width, height = height)) +
    scale_fill_gradientn(name = 'Time', colors = colorRampPalette(c("#1b98e0", "red"))(5)) + theme_classic()+xlab('x')+ylab('y') + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                                                                                                                         axis.text.y=element_blank(),axis.ticks=element_blank(),axis.title.x=element_blank(),
                                                                                                                                         axis.title.y=element_blank(), legend.title=element_text(size=18, face="bold"),legend.text=element_text(size=16))


  if (length(plot_title) >= 1) g <- g + ggtitle(label = plot_title)
  result = list(ordered_pseudotime = ordered_pseudotime, ordered_trajectory = ordered_trajectory, plot = g)
  object@trajectory_plot<- result
  return(object)
}
