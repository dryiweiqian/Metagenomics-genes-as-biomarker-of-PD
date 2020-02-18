## 加载包
## import packages
library(ape)
library(vegan)
library(ggplot2)

## set parameter
file_dist <- "../data/sample.bray.curtis.txt"
file_sample2group <- "../data/sample2group.txt"
output.prefix <- "stage1"
color.palette <- c("orangered3","royalblue")
width <- 8
height <- 8

## read distance file
sample.dist <- read.table(file = file_dist, header = T, sep = "\t", row.names = 1)
sample.dist <- dist(sample.dist)

## read sample2group; with SampleName in the 1st column, SampleGroup in the 2nd column
sample2group <- read.table(file = file_sample2group, header = T, sep = "\t")
dict.sample2group <- setNames(sample2group[[2]], sample2group[[1]])

## PCoA
## 作用：beta 多样性比较，绘制散点图
## sample.bray.curtis.dist：样本间的bray curtis距离矩阵，其中行名为样本，列名为样本
## sample.factor：各样本对应的组别，只能是两组，样本顺序需与上述矩阵相同
## sample.pair: 样本对，先一类样本，再后续另一类样本，样本对在各类样本中保持对应顺序

## 根据 sample.bray.curtis.dist 计算各样本的 mds 坐标
## calculate mds
sample.pcoa <- ape::pcoa(sample.dist)
sample.relative.eig <- sample.pcoa$values$Relative_eig[1:2]
sample.mds <- cmdscale(sample.dist)
sample.mds <- data.frame(sample.mds)
colnames(sample.mds) <- c("PCOA1", "PCOA2")
sample.mds <- data.frame(Sample = rownames(sample.mds), sample.mds)
sample.mds$Phenotype <- dict.sample2group[sample.mds$Sample]
 
## Anosim
library(vegan)
sample.anosim <- anosim(sample.dist, as.factor(sample.mds$Phenotype), 
                        permutations = 1000)
beta.Rsquare <- sample.anosim$statistic
beta.pvalue <- sample.anosim$signif

## 绘制坐标点及分组信息
## plot PCoA
library(ggplot2)
pcoa.graph <- ggplot(data = sample.mds) + 
    geom_point(aes(x = PCOA1, y = PCOA2, color = Phenotype), size = 2) +
    scale_color_manual(values = color.palette) 

## plot paired samples' lines  
if(!is.null(sample.pair)){
  pair.data <- merge(sample.mds, data.frame(Sample=sample.pair, Pair=factor(rep(1:(length(sample.pair)/2),times=2))), by = "Sample", all.x = T)
  
  pcoa.graph <- pcoa.graph +
    geom_line(data = pair.data, 
              aes(x = PCOA1, y = PCOA2, group = Pair), color = "grey")
}
  
## set theme
pcoa.graph <- pcoa.graph +
    theme_bw() +
    theme(legend.justification = c(1,1), 
          legend.position = c(1,1),
          legend.background = element_rect(fill="white", color = "black", linetype="solid"),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.title=element_text(size=15),
          legend.text=element_text(size=14),
          text = element_text(size=15))
  
## add PCoA Relative Eig
pcoa.graph <- pcoa.graph +
    labs(title = paste("ANOSIM R = ", as.character(round(beta.Rsquare, 6)),
                       " , P = ", as.character(round(beta.pvalue, 6)), sep = "")) +
                xlab(paste0("PCoA1(", format(as.numeric(sample.relative.eig[1])*100, digits = 4), "%)")) +
                ylab(paste0("PCoA2(", format(as.numeric(sample.relative.eig[2])*100, digits = 4), "%)"))
  
## output
ggsave(pcoa.graph, filename = paste(output.prefix, ".beta.diversity.pdf", sep = ""),
       width = width, height =  height)

############################################################

