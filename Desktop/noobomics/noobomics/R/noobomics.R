



require(RColorBrewer)
require(graphics)
require(NMF)
require(ggbiplot)
require(stats)

noobomics <- function(df1, df2){
  #create color spectrum for heatmap
  neg_spectrum <- colorRampPalette(brewer.pal(11,'Spectral'))(50)
  
  #environment being set is due to a bug in scope of ggplot2 calls inside of other functions
  environment = environment()
  
  #pre-process data into useable melted together dataframe
  matched_ids <- intersect(df1[,1], df2[,1])
  matched_df1_data <- subset(df1, df1[,1] %in% matched_ids == TRUE)
  matched_df2_data <- subset(df2, df2[,1] %in% matched_ids == TRUE)
  combined_dfs <- cbind(matched_df1_data, matched_df2_data)

  #perform log2() on columns using [,2] and [,5] as reference(WT) as those should be 1st column of real data in each df
  #TODO - function that makes this user controlled as to how they transform data, but for now defaults to log2(sample/ref) relative to [,2] & [,6]
  log2_transformed_combined_df <- combined_dfs
  log2_transformed_combined_df[,3] <- log2(combined_dfs[,3])-log2(combined_dfs[,2])
  log2_transformed_combined_df[,4] <- log2(combined_dfs[,4])-log2(combined_dfs[,2])
  log2_transformed_combined_df[,2] <- log2(combined_dfs[,2])-log2(combined_dfs[,2])
  log2_transformed_combined_df[,7] <- log2(combined_dfs[,7])-log2(combined_dfs[,6])
  log2_transformed_combined_df[,8] <- log2(combined_dfs[,8])-log2(combined_dfs[,6])
  log2_transformed_combined_df[,6] <- log2(combined_dfs[,6])-log2(combined_dfs[,6])

  #Histograms with density of distributions
  columns_to_use <- c(3,4,7,8)
  par(mfrow = c(2,2))
  sapply(columns_to_use, function(i){ 
    title_i <- paste(c("Distribution of Log2(", paste(names(log2_transformed_combined_df[i])), ")"), collapse = "")
    h <- hist(log2_transformed_combined_df[,i],  main = title_i, xlab="Log2(Counts/Reference)", col=i, border="black")
    xfit<-seq(min(log2_transformed_combined_df[,i]),max(log2_transformed_combined_df[,i]),length=40)
    yfit<-dnorm(xfit,mean=mean(log2_transformed_combined_df[,i]),sd=sd(log2_transformed_combined_df[,i]))
    yfit <- yfit*diff(h$mids[1:2])*length(log2_transformed_combined_df[,i])
    lines(xfit, yfit, col="blue", lwd=2)
  })

  #scatterplots of Pearson CCs
  par(mfrow = c(2,3))
  sapply(seq(1,3,1), function(i){ 
    xrange <- find_plot_range(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+1]])
    yrange <- find_plot_range(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+1]])
    mainlabel <- paste(c(" Pearson \nR=",signif(cor(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+1]], use = "everything", method = "pearson"),4)),collapse="")
    plot(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+1]], col=i,main=mainlabel,xlim=xrange, ylim=yrange,xlab=paste(c("Log2(", paste(names(log2_transformed_combined_df[columns_to_use[i]])), ")"), collapse = ""),ylab=paste(c("Log2(", paste(names(log2_transformed_combined_df[columns_to_use[i+1]])), ")"), collapse = ""))
    #lines(coef(lm(log2_transformed_combined_df[,columns_to_use[i+1]] ~ log2_transformed_combined_df[,columns_to_use[i]])),lwd=2,lty=1,xlim=xrange, ylim=yrange)
    #legend("topright",legend=paste(c("Pearson CC\nR=",signif(cor(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+1]], use = "everything", method = "pearson"),4)),collapse=""))
    abline(c(0,1), lty=2)
  })
  sapply(seq(1,2,1), function(i){ 
    xrange <- find_plot_range(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+2]])
    yrange <- find_plot_range(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+2]])
    mainlabel <- paste(c(" Pearson \nR=",signif(cor(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+2]], use = "everything", method = "pearson"),4)),collapse="")
    plot(log2_transformed_combined_df[,columns_to_use[i]],log2_transformed_combined_df[,columns_to_use[i+2]],col=i+3,main=mainlabel,xlim=xrange, ylim=yrange,xlab=paste(c("Log2(", paste(names(log2_transformed_combined_df[columns_to_use[i]])), ")"), collapse = ""),ylab=paste(c("Log2(", paste(names(log2_transformed_combined_df[columns_to_use[i+2]])), ")"), collapse = ""))
    abline(c(0,1), lty=2)
  })
  xrange <- find_plot_range(log2_transformed_combined_df[,columns_to_use[1]],log2_transformed_combined_df[,columns_to_use[4]])
  yrange <- find_plot_range(log2_transformed_combined_df[,columns_to_use[1]],log2_transformed_combined_df[,columns_to_use[4]])
  mainlabel <- paste(c(" Pearson \nR=",signif(cor(log2_transformed_combined_df[,columns_to_use[1]],log2_transformed_combined_df[,columns_to_use[4]], use = "everything", method = "pearson"),4)),collapse="")
  plot(log2_transformed_combined_df[,columns_to_use[1]],log2_transformed_combined_df[,columns_to_use[4]],col=6,main=mainlabel,xlim=xrange, ylim=yrange,xlab=paste(c("Log2(", paste(names(log2_transformed_combined_df[columns_to_use[1]])), ")"), collapse = ""),ylab=paste(c("Log2(", paste(names(log2_transformed_combined_df[columns_to_use[4]])), ")"), collapse = ""))
  abline(c(0,1), lty=2)
  
  #display only 1 plot at a time again, and plot heatmap of results  
  par(mfrow = c(1,1))
  aheatmap(log2_transformed_combined_df[,c(3,4,7,8)],main="Heatmap of Unscaled Log2() Data",scale = 'none',labRow = log2_transformed_combined_df[,1], color = rev(neg_spectrum), cexRow = 1.5,cellwidth = 20, Colv = FALSE)
  
  #calculate principal components and significant labels, eventually significance expression will be user controlled
  log2_transformed_combined_df.pca <- princomp(log2_transformed_combined_df[,c(3,4,7,8)])
  significant_hit_labels <- subset(log2_transformed_combined_df, abs(log2_transformed_combined_df.pca$scores[,1]) > 0.5 & abs(log2_transformed_combined_df.pca$scores[,2]) > 0.5)
  significant_hit_labels_complete <- c(sapply(seq(1,nrow(log2_transformed_combined_df),1), function(i) ifelse(log2_transformed_combined_df[i,1] %in% significant_hit_labels[,1], yes=as.vector(log2_transformed_combined_df[i,1]), no="*")))
  
  #plot the biplot in a nicely formatted way
  ggbiplot(pcobj = log2_transformed_combined_df.pca, environment = environment(), scale=1, obs.scale=1, circle=T, labels=significant_hit_labels_complete,labels.size = 3, varname.size = 5)+
    theme(panel.background = element_rect(fill="white",colour = "black"))+
    scale_x_continuous(limits = c(range(log2_transformed_combined_df.pca$scores[,1])[1]-0.5,range(log2_transformed_combined_df.pca$scores[,1])[2]+0.5))+
    scale_y_continuous(limits = c(range(log2_transformed_combined_df.pca$scores[,2])[1]-0.5,range(log2_transformed_combined_df.pca$scores[,2])[2]+0.5))+
    labs(title="Principal Component Analysis of Log2 Transformed Data")
  
}



















