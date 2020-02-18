## 加载包
## import packages
library(ape)
library(cluster)
library(ggplot2)
library(ggrepel)
library(clusterSim)
library(data.table)

## 配置参数
## set parameters
path.gene.abund <- "../data/gene.abundance.txt"
path.gene.genus <- "../data/gene2genus.txt"
path.sample.group <- "../sample2group.txt"
output.prefix <- "enterotype"
cluster.number <- 2
cluster.labels <- c("Bacteriodes", "Prevotella")
## 这里选择 cluster.number = 2 而不是 3，因为两者 CH 指标差不多，而当其为2时 silhouette 指标更大
## we set cluster.number to 2 not 3, as it has similar CH index, and greater silhouette index
plot.width <- 8
plot.height <- 8

## 导入基因丰度数据
## import gene abundance data
gene.abund <- fread(file = path.gene.abund, header = T, sep = "\t", data.table = F)
rownames(gene.abund) <- gene.abund[[1]]
gene.abund <- gene.abund[, 2:ncol(gene.abund)]

## 导入样本信息表
## import sample:group information
sample.group <- fread(file = path.sample.group, header = T, sep = "\t", data.table = F)

## 导入基因种属对应表
## import gene:genus information
gene.genus <- fread(file = path.gene.genus, header = F, data.table = F)
colnames(gene.genus) <- c("gene", "genus")
gene.genus.dict <- setNames(gene.genus$genus, gene.genus$gene)

## 取交集
## get the intersects of gene abundances and gene annotations
genus.gene.abund <- gene.abund[rownames(gene.abund) %in% names(gene.genus.dict), ]
genus.gene.abund_genus <- gene.genus.dict[rownames(genus.gene.abund)]

## 计算种属丰度表
## get genus abundance by summing up the according genes
genus.abund <- aggregate(x = genus.gene.abund, by = data.frame(genus.gene.abund_genus), FUN = sum)
rownames(genus.abund) <- genus.abund[[1]]
genus.abund <- genus.abund[, -1]

## 按样本归一化
## normalize genus in each sample
genus.abund <- apply(genus.abund, 2, function(i) i/sum(i))
 


## reference: http://enterotype.embl.de/enterotypes.html
## Jensen 距离
## Jensen distance
dist.JSD <- function(inMatrix, 
                     pseudocount=0.000001, ...) {
  
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

## PAM 聚类函数
## PAM clustering
pam.clustering=function(x,k) { 
  # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

## 作用：PCoA 可视化
## mymds：cmdscale(dist) 的输出结果
## PCoA 
pcoa.plot <- function(sample.dist,
                      sample.group1,
                      sample.group2,
                      name.sample.group1 = deparse(substitute(sample.group1)),
                      name.sample.group2 = deparse(substitute(sample.group2)),
                      cluster.labels = NULL){
  
  ## PCoA
  library(ape)
  sample.pcoa <- ape::pcoa(sample.dist)
  sample.mds <- sample.pcoa$vectors[, 1:2]
  sample.relative.eig <- sample.pcoa$values$Relative_eig[1:2]
  sample.mds <- data.frame(sample.mds)
  sample.order <- rownames(sample.mds)
  colnames(sample.mds) <- c("PC1", "PC2")
  sample.mds <- data.frame(Sample = rownames(sample.mds), sample.mds)
  
  ## 计算中心点
  ## calculate cluster center
  sample.mds$sample.group1 <- sample.group1
  colnames(sample.mds)[ncol(sample.mds)] <- name.sample.group1
  centroids <- aggregate(formula = cbind(PC1, PC2) ~ get(name.sample.group1), data = sample.mds, FUN = mean)
  colnames(centroids)[1] <- name.sample.group1
  sample.mds <- merge(sample.mds, centroids, by = name.sample.group1, suffixes=c("",".centroid"))
  sample.mds <- sample.mds[match(sample.order, sample.mds$Sample), ]
  
  ## 绘制坐标点及分组信息
  ## plot pcoa
  library(ggplot2)
  cbbPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7",  
                  "#56B4E9", "#F0E442", "#009E73", "#000000")
  sample.mds$sample.group2 <- sample.group2
  colnames(sample.mds)[ncol(sample.mds)] <- name.sample.group2

  pcoa.graph <- ggplot(data = sample.mds, 
                       aes(x = PC1, y = PC2, color = as.character(sample.group1))) +
                geom_point(aes(shape = as.character(sample.group2)), size = 3) +
                scale_shape_manual(values = c(1, 17)) +
                scale_color_manual(values = cbbPalette) +
                labs(color = name.sample.group1, shape = name.sample.group2) + 
                theme_bw() +
                theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                      plot.margin=unit(c(1.1,1.1,1.1,1.1),"cm"),
                      legend.background = element_rect(size=.5, color = "black"),
                      axis.title=element_text(size=14)) +
                xlab(paste0("PC1(", format(as.numeric(sample.relative.eig[1])*100, digits = 4), "%)")) +
                ylab(paste0("PC2(", format(as.numeric(sample.relative.eig[2])*100, digits = 4), "%)"))
  
  ## 添加中心圆圈
  ## add clustering circles
  pcoa.graph <- pcoa.graph +
    geom_segment(aes(x = PC1.centroid, y = PC2.centroid, 
                     xend = PC1, yend = PC2, 
                     color = as.character(sample.group1))) +
    geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
    geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    stat_ellipse(type = "norm", linetype=2, level = 0.70)
  
  if(!is.null(cluster.labels)){
    library(ggrepel)
    pcoa.graph <- pcoa.graph + 
      geom_label_repel(data = centroids, aes(x = PC1, y = PC2, fill = cluster.labels, label = cluster.labels), 
                       color = 'white', box.padding = unit(0.25, "lines"), show.legend=FALSE) +
      scale_fill_manual(values = cbbPalette)
  }

  pcoa.graph
}


## 计算 JS 距离 
## calculate JSD
genus.JSD.dist <- dist.JSD(genus.abund)
  
## 计算 CH 指标
## calculate CH index
nclusters <- NULL
sclusters <- NULL
choose.cluster.number <- min(20, ncol(gene.abund)-1)
for (k in 1:choose.cluster.number) {
    if (k==1) {
      nclusters[k] <- NA
      sclusters[k] <- NA
    } else {
    genus.cluster_temp <- pam.clustering(genus.JSD.dist, k)
    nclusters[k] <- index.G1(t(genus.abund), genus.cluster_temp,
                             d = genus.JSD.dist, centrotypes = "medoids")
    genus.cluster_temp <- pam(genus.JSD.dist, k, diss = T)
    sclusters[k] <- mean(silhouette(genus.cluster_temp$clustering,
                                            genus.JSD.dist)[,"sil_width"])
    }
}

## 如果没有设置 cluster.number, 则取CH指标最大的聚类个数
## if cluster.number not set, then choose that maximize the CH index
nclusters[is.na(nclusters)] <- 0
if((! exists("cluster.number")) | is.null(cluster.number)){
  cluster.number <- which(nclusters == max(nclusters))
}
genus.cluster <- pam(as.dist(genus.JSD.dist), cluster.number, diss=TRUE)

## silhouette验证指标（越接近1越好）
genus.cluster.silhouette <- mean(silhouette(genus.cluster$clustering,
                                            genus.JSD.dist)[,"sil_width"])


## pcoa图
Enterotype <- genus.cluster$clustering[colnames(gene.abund)]
Phenotype <- sample.group$Group
enterotype.plot <- pcoa.plot(genus.JSD.dist, Enterotype, Phenotype, cluster.labels = cluster.labels)
my.x <- enterotype.plot$data$PC1
my.y <- enterotype.plot$data$PC2
  
enterotype.plot <- enterotype.plot +
                   labs(title = paste("Average Silhouette Width: ", 
                                as.character(round(genus.cluster.silhouette, 3), 
                                sep = ""))) +
                   theme(plot.title = element_text(hjust = 0.5, size = 20),
                         legend.title = element_text(size=15),
                         legend.text = element_text(size=12
),
                         text = element_text(size=20))
  
## 输出图表
### enterotype 聚类指标图
pdf(paste(output.prefix, ".CHindex.pdf", sep = ""))
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
dev.off()

pdf(paste(output.prefix, ".MeanSilindex.pdf", sep = ""))
plot(sclusters, type="h", xlab="k clusters", ylab="Mean Silhouette index")
dev.off()

### enterotype JSD 距离的 PCOA 图
ggsave(enterotype.plot, filename = paste(output.prefix, ".pcoa.pdf", sep = ""), width = plot.width, height = plot.height)
ggsave(enterotype.plot, filename = paste(output.prefix, ".pcoa.png", sep = ""), width = plot.width, height = plot.height)

### enterotype PCOA 向量
write.table(enterotype.plot$data,
            file = paste(output.prefix, ".result.xls", sep = ""),
            row.names = F, col.names = T, sep = "\t", quote = F)
  
### genus 丰度
write.table(genus.abund,
            file = paste(output.prefix, ".genus.abund.xls", sep = ""),
            row.names = T, col.names = NA, sep = "\t", quote = F)