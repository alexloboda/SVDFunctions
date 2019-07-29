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

cats <- function(file, depth, ...) {
  cat(spaces(depth), ..., file = file, sep = "")
}

writePopulationStructure <- function(tree, classes, fd, depth) {
  if (length(tree) == 1) {
    cats(fd, 0, "cluster:\n")
    id <- tree[[1]]
    cats(fd, depth + 1, "id: ", id, "\n")
    cats(fd, depth + 1, "name: ", classes[id], "\n")
  } else {
    cats(fd, 0, "split:\n")
    for (i in 1:length(tree)) {
      cats(fd, depth + 1, "- ")
      writePopulationStructure(tree[[i]], classes, fd, depth + 2)
    }
  }
}

plainTree <- function(tree, classes) {
  if (length(tree) == 1) {
    which(classes == tree[[1]])
  } else {
    list(plainTree(tree[[1]], classes), plainTree(tree[[2]], classes))
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
  hierarchy <- plainTree(hierarchy, classes)
  if(length(left) != 0) {
    stop(paste0("Following classes are not mentioned in herarchy: ", left))
  }
  structure(list(classes = classes, samples = classification,
                 hier = hierarchy), class = "clustering")
}

#' @export
print.clustering <- function(x, ...) {
  cat("Classes: ", paste(x$classes, collapse = ", "), "\n\n")
  writePopulationStructure(x$hier, x$classes, stdout(), 0)
}

#' @export
`[<-.clustering` <- function(x, i = NULL, j = NULL, value) {
  if (i %in% x$classes) {
    x$classes[x$classes == i] <- value
    return(x)
  } 
  if (i >= 1 && i <= length(x$classes)) {
    x$classes[i] <- value
    return(x)
  }
  stop(paste0("No such cluster: ", i))
}