checkTree <- function(t, cs) {
  if (is.list(t)) {
    if (length(t) != 2) {
      stop("Tree must be binary")
    }
    cs <- checkTree(t[[1]], cs)
    checkTree(t[[2]], cs)
  } else {
    if (is.null(t) || is.na(t)) {
      stop("Node must not be null nor NA")
    }
    if (!(t %in% cs)) {
      stop(paste("Cluster in a tree must be mentioned exactly once:"), t)
    }
    cs[cs != t]
  }
}

deletionName <- "###____toDelete____$$$"

cats <- function(file, depth, ...) {
  cat(spaces(depth), ..., file = file, sep = "")
}

writePopulationStructure <- function(tree, classes, fd) {
  tree$Do(function(node) {
    if (node$level == 1) {
      cats(fd, 0, "hierarchy:\n")
      return()
    }
    margin <- node$level - 2
    cats(fd, margin, "", ifelse(node$level == 2, "  ", "- "))
    if (data.tree::isLeaf(node)) {
      cats(fd, 0, "cluster:\n")
      cats(fd, margin + 2 , "id: ", node$id, "\n")
      cats(fd, margin + 2, "name: ", node$name, "\n")
    } else {
      cats(fd, 0, "split:\n")
    }
  })
}

recTree <- function(tree, hier, classes) {
  if (length(hier) == 1) {
    leaf <- tree$AddChild(hier[[1]])
    leaf$id <- which(classes == hier[[1]])
  } else {
    node <- tree$AddChild("-#")
    node$id <- ""
    recTree(node, hier[[1]], classes)
    recTree(node, hier[[2]], classes)
  }
}

#' Constructor for clustering object
#' @param classification named with sample names character vector of 
#' classes.
#' @param hierarchy tree-like structure, a list containing cluster name(leafs) or
#' two lists with the same structure. See examples
#' @examples
#' samples <- c("C1", "C1", "C2", "C2", "C3", "C4")
#' names(samples) <- paste0("sample", 1:6)
#' hier <- list("C1", list(list("C2", "C4"), "C3"))
#' clustering(samples, hier)
#' @export
clustering <- function(classification, hierarchy) {
  classes <- as.character(unique(classification))
  map <- setNames(1:length(classes), classes)
  classification <- setNames(map[classification], names(classification))
  left <- checkTree(hierarchy, classes)
  if(length(left) != 0) {
    stop(paste0("Following classes are not mentioned in herarchy: ", left))
  }
  
  tree <- data.tree::Node$new("Population structure")
  recTree(tree, hierarchy, classes)
  structure(list(classes = classes, samples = classification,
                 hier = tree), class = "clustering")
}

normalizeClustering <- function(clustering) {
  ids <- clustering$hier$Get(function(x) x$id, filterFun = data.tree::isLeaf)
  ids <- ids[order(ids)]
  map <- list()
  for (i in 1:length(ids)) {
    map[ids[i]] <- i
  }
  clustering$classes <- names(ids)
  samples <- names(clustering$samples)
  clustering$samples <- unlist(map[clustering$samples])
  names(clustering$samples) <- samples
  clustering$hier$Do(function(node) node$id = map[[node$id]], 
                     filterFun = data.tree::isLeaf)
  clustering
}

clusterID <- function(clustering, cluster) {
  ret <- which(clustering$classes == cluster) 
  if (length(ret) == 0) {
    stop("No such cluster")
  }
  ret[1]
}

findNode <- function(tree, id) {
  ret <- tree$Get(function(node) {
    node
  }, filterFun = function(node) { node$isLeaf && node$id == id })
  ret[[1]]
}

#' Remove cluster from clusutering object with all samples in it
#' @param clustering clustering object 
#' @param a name of the cluster to delete or its id
#' @export
removeCluster <- function(clustering, cluster) {
  clustering$hier <- data.tree::Clone(clustering$hier)
  id <- clusterID(clustering, cluster)  
  tree <- clustering$hier
  node <- findNode(tree, id)
  
  if (node$parent$isRoot) {
    stop("Can't delete the root of hierarchy tree")
  }
  
  samples <- node$Get(function(node) node$id, filterFun = data.tree::isLeaf)
  clustering$samples <- clustering$samples[!(clustering$samples %in% samples)]
  
  sibling <- node$siblings[[1]]
  
  p <- node$parent
  p$name <- deletionName
  gp <- p$parent
  p$AddSiblingNode(sibling)
  gp$RemoveChild(deletionName)
  gp$id <- ""
  
  clustering$samples <- clustering$samples[clustering$samples != id]
  normalizeClustering(clustering)
}

#' Merge a cluster with its sibling in a tree
#' @param clsutering clustering object
#' @param cluster a name of the cluster or its id
#' @export
mergeCluster <- function(clustering, cluster) {
  clustering$hier <- data.tree::Clone(clustering$hier)
  id <- clusterID(clustering, cluster)
  node <- findNode(clustering$hier, id)
  
  if (node$parent$isRoot) {
    stop("Can't merge the root of hierarchy tree")
  }
  
  p <- node$parent
  clusters <- p$Get(function(node) node$id, filterFun = data.tree::isLeaf)
  clusterNames <- p$Get(function(node) node$name, filterFun = data.tree::isLeaf)
  id <- clusters[1]
  clustering$samples[clustering$samples %in% clusters] <- id
  
  gp <- p$parent
  p$name <- deletionName
  newNode <- data.tree::Node$new(paste(clusterNames, collapse = " | "))
  newNode$id <- id
  p$AddSiblingNode(newNode)
  gp$RemoveChild(deletionName)
  gp$id <- ""
  
  normalizeClustering(clustering)
}

#' @export
print.clustering <- function(x, ...) {
  cat("Classes: ", paste(x$classes, collapse = ", "), "\n\n")
  print(x$hier, "id")
}

#' @export
`[<-.clustering` <- function(x, i = NULL, j = NULL, value) {
  if (is.null(i)) {
    if (length(value) != length(x$classes)) {
      stop("number of items to replace is not a multiple of replacement length")
    }
    for (j in 1:length(value)) {
      x[j] <- value[j]
    }
    return(x)
  }
  if (i %in% x$classes) {
    x$classes[x$classes == i] <- value
    x$hier$Do(function(node) { node$name <- value }, 
              filterFun = function(node) { node$name == i })
    return(x)
  } 
  if (i >= 1 && i <= length(x$classes)) {
    x$classes[i] <- value
    x$hier$Do(function(node) { node$name <- value }, 
              filterFun = function(node) { data.tree::isLeaf(node) && node$id == i })
    return(x)
  }
  stop(paste0("No such cluster: ", i))
}