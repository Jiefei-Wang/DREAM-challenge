function(object, quantiles) {
  rsme.cor <- function(x, y) {
    return (sqrt(mean((x[lower.tri(x)] - y[lower.tri(y)])^2)))
  }
  
  insitu.genes <- colnames(object@insitu.matrix)
  insitu.cor <- cor(object@insitu.matrix)
  rsme.scores = data.frame("quantile" = quantiles,
                           "score" = 0)
  
  for (quantile in quantiles) {
    gene.thresholds <- sapply(insitu.genes, function(gene) quantile(as.numeric(object@data[gene, object@data[gene, ] > 0]), quantile))
    names(gene.thresholds) <- gsub("\\..*", "", names(gene.thresholds))
    binarized.data <- object@data[insitu.genes, ]
    binarized.data[binarized.data >= 0] <- 0
    binarized.data[apply(object@data[insitu.genes, ], 2, function(cell) cell > gene.thresholds)] <- 1
    rsme.scores$score[which(rsme.scores$quantile == quantile)] <- rsme.cor(insitu.cor, cor(t(binarized.data)))
  }
  best.quantile <- rsme.scores$quantile[which.min(rsme.scores$score)]
  best.gene.thresholds <- sapply(insitu.genes, function(gene) quantile(as.numeric(object@data[gene, object@data[gene, ] > 0]),
                                                                       best.quantile))
  names(best.gene.thresholds) <- gsub("\\..*", "", names(best.gene.thresholds))
  
  binarized.data <- object@data[insitu.genes, ]
  binarized.data[binarized.data >= 0] <- 0
  binarized.data[apply(object@data[insitu.genes, ], 2, function(cell) cell > best.gene.thresholds)] <- 1
  object@binarized.data <- binarized.data
  
  return (object)
}