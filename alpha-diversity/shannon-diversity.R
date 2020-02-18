## 加载包
## load packages
library(coin)
library(reshape)
library(ggplot2)

## 配色板
## color
color.palette <- c("orangered3","royalblue",
                     "#0072B2", "#D55E00", "#E69F00", "#CC79A7",  
                     "#56B4E9", "#F0E442", "#009E73", "#000000")


## 配置参数
## set parameter
file_diversity <- "../data/sample.alpha.diversity.txt"
sample.factor <- c(rep("Healthy", 40), rep("PD", 40))
names(sample.factor) <- colnames(test)
sample.group.uniq <- c("Healthy", "PD")
mycategory <- "" 
outPrefix <- "shannon_index"
mysigvalue <- 0.008375 ## wilcoxon test
height <- 6
width <- 6

data.plot <- read.table(file = file_diversity, header = T, sep = "\t")
colnames(data.plot) <- c("Sample", "Value")

data.plot$Phenotype <- sapply(data.plot$Sample, FUN = function(sample) {
    sample.factor[colnames(category.abund) == sample]
})
data.plot$Value <- as.numeric(data.plot$Value)

main <- ggplot(data = data.plot) +
    geom_boxplot(aes(x=Phenotype, y = Value, fill=Phenotype), lwd=0.5, width=0.6) +
    geom_jitter(aes(x = Phenotype, y = Value, group = Phenotype)) +
    scale_fill_manual(values = color.palette, labels = sample.group.uniq) +
    theme_bw() + ggtitle(mycategory) +
    theme(axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.margin=unit(c(1.5,1.5,1.5,1.5),"cm"),
          plot.title = element_text(hjust = 0.5)) +
    #coord_cartesian(ylim=c(0, max.ylim)) +
    ## 相对丰度 (<1) 取对数，变负数
    #coord_cartesian(ylim=c(max.ylim, 0)) +
    ylab(label = "Shannon Index")


#mysigvalue <- pvalue(wilcox_test(as.numeric(category.abund[mycategory, ]) ~ as.factor(sample.factor), distribution = "exact"))
if(TRUE){
   
   #### 计算差异线所在位置
   x <- which(levels(data.plot$Category) == mycategory)
   x1 <-  1 
   x2 <-  2
   y1 <- 12.5
   y2 <- 12.6
      
   #### 差异指标
   if(mysigvalue < 0.05){
     add.text <- paste("p = ", as.character(format(mysigvalue, digits = 2)), sep = "")
   } 
      
   #### 绘制差异线
   main <- main + 
        geom_path(data = data.frame(x=c(x1,x1,x2,x2), y=c(y1,y2,y2,y1)),
                  aes(x = x, y = y), col = "gray12", size = 0.6) +
        annotate("text", label = add.text, x = (x1+x2)/2, y=12.7, 
                 size = 4, colour = "gray12")
    
}
  

#plot(main)


ggsave(main, filename=paste0(outPrefix, ".boxplot.pdf"), height=height, width=width)
ggsave(main, filename=paste0(outPrefix, ".boxplot.png"), height=height, width=width)
