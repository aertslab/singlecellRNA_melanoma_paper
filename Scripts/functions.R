# 401 DC alignment functions
library(car)

plot_genes_along_trajectory<-function(dgem.inter.expr, dgem.inter.error, pseudo.time, expr.mat, gene, line) {
  selected.inter<-dgem.inter.expr[gene,]
  
  df.i = data.frame(traj = pseudo.time, value=(selected.inter), error=dgem.inter.error[gene,])
  # df = data.frame(traj = pseudo.time, expr = expr.mat[gene,])
  # df.M = melt(df, id.vars = 'traj')
  # Plot of an example gene and its interpolation with error bars
  ggplot(df.i, aes(x=traj,y=value)) + 
    geom_errorbar(aes(ymin=value-error/2, ymax=value+error/2)) + 
    geom_line(size = 2) + 
    # geom_point(data=df.M, aes(x=traj,y=value)) +
    ggtitle(paste0(line, " - ", gene))
}

save_inter_gene_expr_heatmap<-function(x, genes = NULL, metric, file.name, width = 600, height = 2400) {
  x.no.na<-na.omit(x)
  genes.sorted<-names(sort(metric, decreasing = T))
  genes.sorted.filtered<-genes.sorted[genes.sorted%in%row.names(x.no.na)]
  if(!is.null(genes)) {
    genes.sorted.filtered<-genes.sorted.filtered[genes.sorted.filtered%in%genes]
  }
  x.no.na.filtered.sorted<-x.no.na[genes.sorted.filtered,]
  png(file.name, width = width, height = height)
  print(Heatmap(x.no.na.filtered.sorted, name = "ISGE", use_raster = TRUE, cluster_rows = F, cluster_columns = F, show_row_names = F, show_column_names = F))
  dev.off()
  invisible(x.no.na.filtered.sorted)
}

plot_gene_expr_lambda<-function(line, lowess, lowess.fit, dgem, gene) {
  lowess.fit.with.meta.data<-merge(x = lowess.fit, y = meta_data, by = "row.names", all.x = T, all.y = F)
  row.names(lowess.fit.with.meta.data)<-lowess.fit.with.meta.data$Row.names
  gene.expr<-as.data.frame(dgem[gene, ])
  lowess.fit.with.meta.data[row.names(gene.expr), "gene"]<-gene.expr[,1]
  lowess.fit.with.meta.data<-lowess.fit.with.meta.data[with(lowess.fit.with.meta.data, order(gene)), ]
  lowess.fit.with.meta.data[names(lowess$lambda), "PrinCurveLambda"]<-lowess$lambda
  
  p1<-ggplot(data = lowess.fit.with.meta.data, aes(x = DM1, y = DM2, color = gene)) + 
    geom_point(size = 3) +
    ggtitle(paste0(line, " - ", gene, " expression"))
  print(p1)
  p2<-ggplot(data = lowess.fit.with.meta.data, aes(x = DM1, y = DM2, color = PrinCurveLambda)) + 
    geom_point(size = 3) +
    ggtitle(paste0(line, " - lambda"))
  print(p2)
}

# 500 clustering functions

get_first_closest<-function(vec, val) {
  down<-if(sign(vec[length(vec)]-vec[1]) < 0) T else F
  for(i in seq_along(along.with = vec)) {
    if(down) {
      if(vec[i] < val)
        return (i)
    } else {
      if(vec[i] > val)
        return (i)
    }
  }
}

get_half_life<-function(x) {
  return (apply(X = na.omit(x), MARGIN = 1, FUN = function(r) {
    return (get_first_closest(vec = r, val = 0.5)*(abs(r[1]-r[length(r)]))*sign(r[1]-r[length(r)]))
  }))
}

get_half_life_sorted<-function(x) {
  return (sort(x = get_half_life(x = x), decreasing = T))
}

fit_logit<-function(y,x) {
  suppressWarnings({
    gene<-y[1]
    # print(gene)
    y<-as.numeric(y[2:length(y)])
    out<-tryCatch({
      start.vals<-coef(lm(logit(y)~x))
      nls(y~phi1/(1+exp(-(phi2+phi3*x))), start=list(phi1=1,phi2=start.vals[1],phi3=start.vals[2]))
    }, error=function(cond) {
      return(NA)
    })
    diff.start.end<-y[1]-y[length(y)]
    abs.diff.start.end<-abs(diff.start.end)
    sign.diff.start.end<-sign(diff.start.end)
    if(is.na(out))
      return (cbind(data.frame("gene"=NA, "diff.start.end"=abs.diff.start.end, "sign.diff.start.end"=sign.diff.start.end, "phi1"=NA,"phi2"=NA,"phi3"=NA,"std.res.error"=NA, "rank.metric"=NA), t(x = data.frame(rep(NA, length(y))))))
    std.res.error<-sqrt(deviance(out)/df.residual(out))
    rank.metric<-1/std.res.error*abs.diff.start.end*sign.diff.start.end
    pred.y<-coef(out)[1]/(1+exp(-(coef(out)[2]+coef(out)[3]*x)))
    return (cbind(data.frame("gene"=gene, "diff.start.end"=abs.diff.start.end, "sign.diff.start.end"=sign.diff.start.end, "phi1"=coef(out)[1],"phi2"=coef(out)[2],"phi3"=coef(out)[3],"std.res.error"=std.res.error, "rank.metric"=rank.metric), t(x = data.frame(pred.y))))
  })
}

get_fitted_logit_half_life<-function(x, traj) {
  library(tictoc)
  tic()
  x.no.na<-na.omit(x)
  print("Fitting logit for each gene...")
  x.no.na.logit.model<-do.call(what = "rbind", args = apply(X = cbind(data.frame("gene"=row.names(x.no.na)), x.no.na), MARGIN = 1, FUN = fit_logit, x = traj))
  print(dim(x.no.na.logit.model))
  x.no.na.logit.model.no.na<-na.omit(x.no.na.logit.model)
  print(dim(x.no.na.logit.model.no.na))
  print("Computing half-life of fitted logit...")
  x.no.na.logit.model.no.na.half.life<-setNames(object = as.numeric(get_half_life(x = x.no.na.logit.model.no.na[,9:208])), nm = row.names(x.no.na.logit.model.no.na[,9:208]))
  warning(paste0(sum(is.na(x = x.no.na.logit.model.no.na.half.life))), " genes have no numeric (NA) half-life.")
  x.no.na.logit.model.no.na.half.life.sorted<-sort(x = x.no.na.logit.model.no.na.half.life, decreasing = T)
  # By default sort is removing NAs -> subset the logit data.frame
  toc()
  return (list("logit.model"=x.no.na.logit.model.no.na[names(x.no.na.logit.model.no.na.half.life.sorted),], "logit.model.half.life"=x.no.na.logit.model.no.na.half.life.sorted))
}

aligner<-function(query.interpolator, ref.interpolator, sig.calc = F, num.perm = 200) {
  return (list("alignment"=globalAlign(x = query.interpolator$scaledData
                                       , y = ref.interpolator$scaledData
                                       , scores = list(query = query.interpolator$traj, 
                                                       ref = ref.interpolator$traj)
                                       , sigCalc = sig.calc, numPerm = num.perm)
               , "query.interpolator"=query.interpolator
               , "ref.interpolator"=ref.interpolator))
}

scaler<-function(aligner, title) {
  alignment.mapping<-data.frame("query"=aligner$alignment$align[[1]]$index1, "ref"=aligner$alignment$align[[1]]$index2)
  # Apply the transformation of pseudotime values of the query trajectory based on the alignment, i.e.: 
  # each element of the query pseudotime will be transformed into the mean pseudotime values of the reference trajectory which were linked to it by the alignment.
  query.pt.scaled.to.ref.pt<-do.call(what = "c", args = lapply(X = seq_along(aligner$query.interpolator$traj), function(pt) {
    mean(aligner$ref.interpolator$traj[alignment.mapping$ref[alignment.mapping$query == pt]])
  }))
  spl <- smooth.spline(aligner$query.interpolator$traj, query.pt.scaled.to.ref.pt)
  plot(aligner$query.interpolator$traj, query.pt.scaled.to.ref.pt,  xlab = "Raw Query Pseudotime", ylab = "Scaled Query Pseudotime", main = title)
  lines(spl, col = "blue")
  abline(a=0,b=1,col="red",lwd=3)
  return (list("scale"=function(x) { return (predict(spl$fit, x = x)) }))
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

plotCellAlignment<-function(alignment, title = NULL) {
  library(pheatmap)
  costMat = alignment$localCostMatrix
  costMat = t(apply(costMat, 1, function(x) {
    return(as.numeric(x))
  }))
  linearInd = sub2ind(nrow(costMat), alignment$align[[1]]$index1, 
                      alignment$align[[1]]$index2)
  costMat[linearInd] = NA
  costMat = data.frame(costMat, row.names = 1:nrow(costMat))
  colnames(costMat) = 1:ncol(costMat)
  
  if(is.null(title)) {
    title<-"alignment plot"
  }
  
  if (!is.null(alignment$ptShift)) {
    annotCols = data.frame(ptShift = abs(alignment$ptShift), 
                           sign = factor(sign(alignment$ptShift)), row.names = colnames(costMat))
    return(pheatmap(costMat, cluster_cols = F, cluster_rows = F, 
             border_color = NA, main = title, show_rownames = F, 
             show_colnames = F, annotation_col = annotCols))
  }
  else {
    return (pheatmap(costMat, cluster_cols = F, cluster_rows = F, 
             border_color = NA, main = title, show_rownames = F, 
             show_colnames = F))
  }
  return(NA)
}


