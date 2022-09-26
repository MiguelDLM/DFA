#Install packages for first time

require(devtools)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/BaylorEdPsych/BaylorEdPsych_0.5.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/mvnmle/mvnmle_0.1-11.1.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
packageurl <- "https://cran.r-project.org/src/contrib/Archive/MissMech/MissMech_1.0.2.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

#Load libraries

packages <- c("MissMech","kableExtra","parallel","doParallel","stringr","DAAG","scales","grid","gridExtra","cowplot","GGally","Amelia","BaylorEdPsych","caret","data.table","dplyr","finalfit","ggforce","ggfortify","ggplot2","ggpubr","heplots","klaR","MASS","MVN","reshape2")
lapply(packages, require, character.only = TRUE)

#Multicore
no_cores <- detectCores() - 1
cl <- makePSOCKcluster(no_cores)
registerDoParallel(cl)

#Load csv data and subsets on R

setwd("C:/Articulo FDA")
original.data <- read.csv(file = 'DeGusta.csv', na.strings=c("","NA"))
original.data <- original.data [, 2:10]
original.data$Habitat.Category <- as.factor(original.data$Habitat.Category)
#Statistical assuptions

#Normality

split_data <- original.data %>%
  group_split(Habitat.Category)%>%
  setNames(unique(original.data$Habitat.Category))

MVN <- lapply(split_data, function(split_data){
  MVN <- mvn(split_data[,2:9], mvnTest = "hz", covariance = TRUE, desc = TRUE, transform = "none",univariateTest = "SW") 
})  

MVN$L$multivariateNormality$Group <- names(MVN[1])
MVN$O$multivariateNormality$Group <- names(MVN[2]) 
MVN$F$multivariateNormality$Group <- names(MVN[3]) 
MVN$H$multivariateNormality$Group <- names(MVN[4])

multi_normality <- rbind.data.frame(MVN[[1]][[1]],MVN[[2]][[1]],MVN[[3]][[1]],MVN[[4]][[1]])


normality <- cbind.data.frame(lapply(MVN, function(MVN){
  univariate <- MVN$univariateNormality
  
}))

normality<- normality[,c(1,2,3,4,5,8,9,10,13,14,15,18,19,20)]
normality <- kbl(normality) %>%
  kable_classic() %>%
  add_header_above(c(" " = 2, "L" = 3, "O" = 3, "F" = 3, "H"=3))



#Distribution graph

melt.data <- melt(original.data)

distribution <- ggplot(data = melt.data, aes(x = value, fill=Habitat.Category)) + 
  geom_density(alpha = 0.7) + geom_line(stat = "density") +
  scale_fill_manual(values=c("#fde725","#5ec962","#21918c","#3b528b","#440154"),
                    labels = c("Forest", "Heavy Cover", "Light Cover", "Open"))+
  facet_wrap(~variable,scales = "free")+
  labs(x = 'Size in mm', y = 'Density') + theme_bw() + theme(legend.position=c(.85,.15))+
  guides(fill=guide_legend(title = "Habitat Category", ncol=2))

#Homoscedasticity

BoxM_test <- boxM(original.data[, 2:9],original.data[,"Habitat.Category"])
BoxM_test
#Following Jamshidian et al., 2014, an artificial incomplete dataset was created with a
#from a complete dataset, creating different missing pattern for each group
#This incomplete dataset is used in the non-parametric test for
#Homoscedasticity, while the complete dataset is used as "imputed" data
split_data$L$LM <- NA
split_data$O$TD <- NA
split_data$O$WI <- NA
split_data$H$TI <- NA
split_data$H$LM <- NA
split_data$F$TP <- NA
split_data$F$LM <- NA

paramet_dataset <- rbindlist(split_data)

paramet_dataset <- as.matrix(paramet_dataset[,2:9])
original.data2 <- as.matrix(original.data[,2:9])
non_para <- TestMCARNormality(data=paramet_dataset, imputed.data = original.data2)
print(non_para)

#Transform data


log.data <- original.data %>%
  mutate_at(2:9, list(log= ~ log(.)))

log.data<- log.data[, c(1,10:17)]

#PCA

dfs <- list(original.data=original.data,log.data=log.data)

PCAs <- lapply(dfs,function(dfs) {
  df <- prcomp(dfs[, 2:9], center = TRUE, scale = TRUE)
  pca.scores <- cbind(Habitat.Category = dfs$Habitat.Category, as.data.table(df$x))
})

#Plots

#Plot PCA's
PCAs_plot <- lapply(dfs,function(dfs) {
  
  df <- prcomp(dfs[, 2:9], center = TRUE, scale = TRUE)
  pca.scores <- cbind(Habitat.Category = dfs$Habitat.Category, as.data.table(df$x))
  
  percentage <- round(df$sdev^2 / sum(df$sdev^2) * 100, 2)
  percentage <- paste(as.character(percentage), "%", sep="")
  
  plot.pca <- ggplot()+geom_point(data=pca.scores, aes(x=PC1, y=PC2, shape=Habitat.Category, colour=Habitat.Category))+
    theme_classic()+xlab(percentage[1]) + ylab(percentage[2]) + guides(colour = "none", shape="none") +scale_colour_manual(values=c("#fde725","#35b779","#31688e","#440154"))
  
  xhist <- 
    axis_canvas(plot.pca, axis = "x") + 
    stat_density(data = pca.scores,aes(x = PC1, fill=Habitat.Category, alpha=0.7)) +
    scale_fill_manual(values=c("#fde725","#35b779","#31688e","#440154"))
  
  yhist <- 
    axis_canvas(plot.pca, axis = "y", coord_flip = TRUE) + 
    stat_density(data = pca.scores,aes(x = PC2, fill=Habitat.Category, alpha=0.7)) + coord_flip()+
    scale_fill_manual(values=c("#fde725","#35b779","#31688e","#440154"))
  
  final.plot <- plot.pca %>%
    insert_xaxis_grob(xhist, grid::unit(0.2, "null"), position = "top") %>%
    insert_yaxis_grob(yhist, grid::unit(0.2, "null"), position = "right") %>%
    ggdraw()
})

#Get legend

legend.pca <- get_legend(ggplot()+geom_point(data=PCAs$original.data, aes(x=PC1, y=PC2, shape=Habitat.Category, colour=Habitat.Category))+
                           theme_classic()+ theme(legend.position = "bottom")+scale_colour_manual(values=c("#fde725","#35b779","#31688e","#440154"))+
                           guides(colour=guide_legend(title = "Habitat Category"),shape=guide_legend(title = "Habitat Category")))

#Merge plots

Compare_PCA <- plot_grid(PCAs_plot$original.data,
                         PCAs_plot$log.data,
                         align = 'vh',
                         labels = c("A","B"),
                         label_size = 10,
                         hjust = -1,
                         nrow = 1)


Compare_PCA <- plot_grid(legend.pca,Compare_PCA, ncol = 1, rel_heights = c(0.1,1))

y.grob <- textGrob("Principal Component 2", 
                   gp=gpar(fontface="bold", fontsize=10), rot=90)

x.grob <- textGrob("Principal Component 1", 
                   gp=gpar(fontface="bold", fontsize=10))

grid.arrange(arrangeGrob(Compare_PCA, left = y.grob, bottom = x.grob))

#Multiple PCAs

subsets.pca <- lapply(PCAs, function(PCAs){
  subset <- list(eight_vars= PCAs, seven_vars= PCAs[,1:8] ,
                 six_vars= PCAs[, 1:7],five_vars= PCAs[, 1:6],
                 four_vars= PCAs[, 1:5], three_vars= PCAs[, 1:4],
                 two_vars= PCAs[, 1:3])
})

#LDA accuracy for multiple PCAs

LDAs <- lapply(subsets.pca,lapply, function(subsets.pca){
  lda <- lda(Habitat.Category ~ ., data = subsets.pca,prior=c(1,1,1,1)/4, CV=TRUE)
  predictions <- confusion(original.data$Habitat.Category,lda$class)
  percent <- predictions$overall

})

#QDA accuracy for multiple PCAs

QDAs <- lapply(subsets.pca,lapply, function(subsets.pca){
  tryCatch({qda <- qda(Habitat.Category ~ ., data = subsets.pca, prior=c(1,1,1,1)/4, CV=TRUE)
  predictions <- confusion(original.data$Habitat.Category,qda$class)
  percent <- predictions$overall
  },
  error=function(e) NA)
})

#Plots
dfa_list <- list(LDAs=LDAs,QDAs=QDAs)

plot_multi_dfa <-  lapply(dfa_list, function(dfa_list){
  
  compare.tables <- cbind.data.frame(treatment=c("Original data","log"),rbind.data.frame(dfa_list$original.data,dfa_list$log.data))
  melt.data <- melt(compare.tables,id.vars="treatment") 
  plot <- ggplot(melt.data, aes(x= variable, y=value, group= treatment, colour=treatment))+ geom_line(aes(size=0.3))+
    geom_point(aes(x= variable, y=value, colour=treatment,shape=treatment, size=0.4))+ guides(size = "none") + theme_classic()+
    labs(y= "Classification success", x = "Number of components") +scale_y_continuous(labels = scales::percent, limits=c(0.54,0.68)) + 
    scale_x_discrete(labels=c("eight_vars"="8","seven_vars"="7","six_vars"="6", "five_vars"="5","four_vars"="4","three_vars"="3","two_vars"="2"), limits=rev)+ 
    theme(legend.position = "bottom") +scale_colour_manual(values=c("#fde725","#440154"))
  
})

legend.PCAs <-get_legend(plot_multi_dfa$LDAs)
compare_accuracy_dfa <- plot_grid(plot_multi_dfa$LDAs + theme(legend.position = "none"),
                                  plot_multi_dfa$QDAs + theme(legend.position = "none"),
                                  align = 'vh',
                                  labels = c("A","B"),
                                  label_size = 10,
                                  hjust = -1)

compare_accuracy_dfa<- plot_grid(legend.PCAs,compare_accuracy_dfa,ncol = 1, rel_heights = c(0.1,1))

#RDAs


cv_5_rand <- trainControl(method = "LOOCV", search = "random", allowParallel=TRUE)
dfs2 <- list(original_vars = dfs, PCAs = list(original.data= subsets.pca$original.data$four_vars, log.data= subsets.pca$log.data$four_vars))

RDAs <- lapply(dfs2,lapply,function(RDAs){
  
  set.seed(1337)
  fit <- train(Habitat.Category ~ ., data = RDAs, method = "rda",
               trControl = cv_5_rand, tuneLength = 1000)				
  Hp <- as.data.table(fit$results)
  Hp <- Hp %>% arrange (desc(Accuracy))
})

plot_rdas <- lapply(RDAs,lapply, function(RDAs){
  
  RDAs <- RDAs%>% arrange (Accuracy)
  
  plot <-ggplot(RDAs, aes(gamma, lambda)) + geom_point(aes(fill = Accuracy, size= Accuracy),shape=21, alpha=0.7)+
    scale_fill_viridis_c(guide = "legend", limits=c(0.35,0.66)) +
    scale_size_continuous(limits=c(0.35,0.66))+ guides(fill= guide_legend(nrow = 1))+
    theme_classic()+
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none")
  
  ymin <- min(RDAs$Accuracy)
  
  xhist <- 
    axis_canvas(plot, axis = "x") + 
    geom_ribbon(data = RDAs,
              aes(x = gamma, y= Accuracy, ymin=ymin, ymax=Accuracy, alpha=0.5)) +
    geom_area()
  yhist <- 
    axis_canvas(plot, axis = "y", coord_flip = TRUE) + 
    geom_ribbon(data = RDAs,
              aes(x = lambda, y= Accuracy, ymin=ymin, ymax=Accuracy, alpha=0.5)) + coord_flip()
  
  
  plot<- plot %>%
    insert_xaxis_grob(xhist, grid::unit(0.2, "in"), position = "top") %>%
    insert_yaxis_grob(yhist, grid::unit(0.2, "in"), position = "right") %>%
    ggdraw()
})

legend.RDAs <-get_legend(ggplot(RDAs$original_vars$original.data, aes(gamma, lambda)) + geom_point(aes(fill = Accuracy, size= Accuracy),shape=21)+
                           scale_fill_viridis_c(guide = "legend", limits=c(0.35,0.66)) +
                           scale_size_continuous(limits=c(0.35,0.66))+ guides(fill= guide_legend(title= "Classification success",nrow = 1),size=guide_legend(title= "Classification success"))+
                           theme_classic() + theme (legend.position = "bottom"))

compare_RDAs<- plot_grid(plot_rdas$original_vars$original.data + theme(legend.position = "none"),
                         plot_rdas$original_vars$log.data + theme(legend.position = "none"),
                         plot_rdas$PCAs$original.data + theme(legend.position = "none"),
                         plot_rdas$PCAs$log.data + theme(legend.position = "none"),
                         align = 'vh',
                         labels = c("A","B","C","D"),
                         label_size = 10,
                         ncol = 2,
                         hjust = -1)

compare_RDAs<- plot_grid(legend.RDAs,compare_RDAs,ncol = 1, rel_heights = c(0.1,1))
y.grob <- textGrob("Lambda", 
                   gp=gpar(fontface="bold", fontsize=10), rot=90)

x.grob <- textGrob("Gamma", 
                   gp=gpar(fontface="bold", fontsize=10))

grid.arrange(arrangeGrob(compare_RDAs, left = y.grob, bottom = x.grob))

#Table RDAs

top_accuracy <- lapply(RDAs, lapply, function(RDAs){
  
  top <- RDAs[1,]
  
})

compare.table.RDA <- cbind.data.frame(Treatment=c("Original data","Log","PCA-original.data","PCA-log"),rbind.data.frame(
  top_accuracy$original_vars$original.data,top_accuracy$original_vars$log.data,
  top_accuracy$PCAs$original.data,top_accuracy$PCAs$log.data))
compare.table.RDA$Accuracy <- percent(compare.table.RDA$Accuracy, accuracy = 0.01)

#LDA

LDA <- lapply(dfs2,lapply,function(dfs2){

  lda <-  lda(Habitat.Category ~ ., data = dfs2, prior=c(1,1,1,1)/4, crossval = TRUE)
  predictions<- lda %>% predict(dfs2)
  mean <- mean(predictions$class == dfs2$Habitat.Category)
  percent <- percent(mean,accuracy = 0.01)
})

#QDA
QDA <- lapply(dfs2,lapply,function(dfs2){
  
  qda <-  qda(Habitat.Category ~ ., data = dfs2, prior=c(1,1,1,1)/4, crossval = TRUE)
  
  predictions<- qda %>% predict(dfs2)
  mean <- mean(predictions$class == dfs2$Habitat.Category)
  percent <- percent(mean,accuracy = 0.01)
})

#Compare tables

dfas <- cbind.data.frame(compare.table.RDA[1:4],Accuracy=t(cbind.data.frame(QDA)),Accuracy=t(cbind.data.frame(LDA)))

rownames(dfas) <- NULL

compare_dfas <- kbl(dfas) %>%
  kable_classic() %>%
  add_header_above(c(" " = 1, "RDA" = 3, "QDA" = 1, "LDA" = 1))

#QDA
qda.res <- qda(Habitat.Category ~ ., data = original.data, prior=c(1,1,1,1)/4, crossval = TRUE)
predictions.qda<- qda.res %>% predict(original.data)
ct <- table(original.data$Habitat.Category, predictions.qda$class)
mean(predictions.qda$class == original.data$Habitat.Category)


#Summarize data
sum <- lapply (as.list(original.data[,2:9]), Summarize)


lda.plot <- partimat(Habitat.Category ~ ., data = original.data, method = "lda", main="lda")
qda.plot <- partimat(Habitat.Category ~ ., data = original.data, method = "qda", main="qda")