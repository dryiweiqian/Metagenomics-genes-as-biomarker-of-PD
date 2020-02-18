## load packages
library(pROC)
library(mRMRe)
library(e1071)
library(data.table)


## mRMR
mRMR <- function(category.abund, 
                 sample.factor, 
                 gene.number,
                 thread.number){
  
  ## 数据准备和参数设置
  library(mRMRe)
  set.thread.count(thread.number)
  data.mRMR <- data.frame(group = factor(sample.factor, ordered = T), t(category.abund))
  data.mRMR <- mRMR.data(data = data.mRMR)
  
  ## mRMR 算法
  mRMR.result <- mRMR.classic(data = data.mRMR, target_indices = c(1), 
                              feature_count = gene.number)
  mRMR.selected.index <- solutions(mRMR.result)[[1]]
  
  ## 返回结果
  ## 没想到，mRMR 增加了 group 一列，导致 index 整体后移，不知结果怎样
  data.frame(selected.index = mRMR.selected.index - 1,
             selected.scores = scores(mRMR.result),
             selected.names = rownames(category.abund)[mRMR.selected.index - 1])
}


## ROC plot
plot.ROC <- function(response.real,
                     response.decision.value,
                     output.prefix,
                     boot.number = 1000){
  
  library(pROC)
  ## 计算 roc 曲线信息
  response.real <- as.numeric(factor(response.real)) - 1
  result.roc <- roc(response.real, response.decision.value, ci=TRUE)
  AUC <- paste("AUC = ", as.character(as.numeric(result.roc$auc), sep = ""))
  CI <- paste("Confidence Interval: ", 
              as.character(round(100 * as.numeric(result.roc$ci)[1], 1)),
              "-", 
              as.character(round(100*as.numeric(result.roc$ci)[3], 1)), 
              "%", sep = "")
  
  ## 绘制 roc 曲线
  pdf(paste(output.prefix, ".roc.pdf", sep = ""))
  plot.roc(result.roc, add = FALSE, axes = TRUE, col="red", 
           legacy.axes = TRUE, xlab = "1-Specificity")
  # 添加 CI 信息
  legend.temp <- legend("bottomright", legend = c(" ", " "), 
                        text.width = abs(strwidth(CI)), xjust = 1, bty = "n")
  text(legend.temp$rect$left + legend.temp$rect$w, legend.temp$text$y, 
       c(AUC, CI), pos = 2)
  # 添加误差限
  ci.sp.obj<-ci.se(result.roc, c(seq(0,1,0.01)), boot.n = boot.number)
  plot(ci.sp.obj, type="bar", col="black")
  dev.off()
}

## set parameter
file_mgs <- "../data/mgs.genus.combined.results.txt"
file_gene2mgs <- "../data/mgs.gene.txt"
file_gene_abund <- "../data/gene.abundance.txt"
sample.factor <- c(rep("Healthy",40),rep("PD",40))
number.marker <- 25
threads <- 2

## read gene abund
gene.abund <- fread(file_gene_abund, header = T, sep = "\t", data.table = F)
rownames(gene.abund) <- gene.abund[[1]]
gene.abund <- gene.abund[, 2:ncol(gene.abund)]

## filter mgs
mgs.result <- fread(file = file_mgs, header = T, sep = "\t", data.table = F)
selected.mgs.names <- mgs.result$mgs[with(mgs.result, 
                                          gene.number.stage1 >= 50 & 
                                          gene.number.stage1 < 3000 & 
                                          species.ratio.stage1 >= 0.9 &
                                          anno.species.stage1 != "Ruthenibacterium")]

selected.mgs.names <- mgs.result$mgs[with(mgs.result, 
                                          gene.number.stage1 >= 50 & 
                                          anno.species.number.stage1 == 1 & 
                                          species.ratio.stage1 >= 0.9)]

## filter genes based on mgs
gene2mgs <- fread(file_gene2mgs, header = T, sep = "\t", data.table = F)
colnames(gene2mgs) <- c("mgs", "gene")
selected.mgs.genes <- gene2mgs$gene[gene2mgs$mgs %in% selected.mgs.names]


## mRMR
selected.markers <- mRMR(gene.abund[rownames(gene.abund) %in% selected.mgs.genes, ], sample.factor, number.marker, threads)


## SVM model
train.data <- t(gene.abund[rownames(gene.abund) %in% selected.markers$selected.names, ])
train.y <- factor(sample.factor, levels = unique(sample.factor), ordered = T)

library(e1071)
train.svm.model <- svm(x = train.data, y = train.y, kernel = "linear")

data.frame(table(predict = predict(train.svm.model, train.data), observe = train.y))

## ROC
train.decision.value <- attr(predict(train.svm.model, train.data, decision.values = T), "decision.values")
train.decision.value <- as.vector(train.decision.value[, 1])
plot.ROC(train.y, train.decision.value, 
         paste("marker-based-on-upmgs", ".trainset", sep = ""))


