setwd("E:/ADI BUANA/seminar")

dataset <- read.csv("miskin.csv")
prov <- read.csv("provinsi.csv")
rownames(dataset) <- prov$wilayah

library(factoextra)
library(cluster)
library(psych)
library(corrplot)

round(cor(dataset),3) 

corrplot(cor(dataset),type="lower",method="circle")
pairs.panels(dataset)

summary(dataset)

kmo <- function(x)
{
  x <- subset(x, complete.cases(x)) # Omit missing values
  r <- cor(x) # Correlation matrix
  r2 <- r^2 # Squared correlation coefficients
  i <- solve(r) # Inverse matrix of correlation matrix
  d <- diag(i) # Diagonal elements of inverse matrix
  p2 <- (-i/sqrt(outer(d, d)))^2 # Squared partial correlation coefficients
  diag(r2) <- diag(p2) <- 0 # Delete diagonal elements
  KMO <- sum(r2)/(sum(r2)+sum(p2))
  MSA <- colSums(r2)/(colSums(r2)+colSums(p2))
  return(list( KMO = KMO, MSA = MSA))
}

kmo(dataset)

Bartlett.spher.test <- function(x)
{
  method <- "Bartlett's test of sphericity"
  data.name <- deparse(substitute(x))
  x <- subset(x, complete.cases(x)) # Omit missing values
  n <- nrow(x)
  p <- ncol(x)
  chisq <- (1-n+(2*p+5)/6)*log(det(cor(x)))
  df <- p*(p-1)/2
  p.value <- pchisq(chisq, df, lower.tail=FALSE)
  names(chisq) <- "X-squared"
  names(df) <- "df"
  return(structure(list(statistic=chisq, 
                        parameter=df, p.value=p.value,
                        method=method, data.name=data.name), 
                   class="htest"))
}

Bartlett.spher.test(dataset)


fviz_nbclust(dataset, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)


# 1. Elbow Method
fviz_nbclust(dataset, kmeans, method = "wss") 

fviz_nbclust(dataset, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")

## 2. Metode silhouette
fviz_nbclust(dataset, kmeans, method = "silhouette") 

fviz_nbclust(dataset, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

## 3. Gap Statistic
set.seed(123)
gap_stat <- clusGap(dataset, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 38) 
fviz_gap_stat(gap_stat)

fviz_nbclust(dataset, kmeans, nstart = 25, method = "gap_stat", nboot = 38)+
  labs(subtitle = "Gap statistic method")

set.seed(222)
km.res <- kmeans(dataset, 4, nstart = 25)
print(km.res)


fviz_cluster(km.res, data = dataset,
             palette = c("green", "blue", "red", "brown"),
             ellipse.type = "euclid", # Concentration ellipse
             star.plot = T, # Add segments from centroids to items
             repel = T, # Avoid label overplotting (slow)
             ggtheme = theme_minimal()
)


res.dist <- dist(dataset, method = "euclidean")
res.hc <- hclust(d = res.dist, method = "ward.D2")

fviz_dend(res.hc, k = 4, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("green", "purple", "blue", "red"),
          color_labels_by_k = TRUE, # color labels by groups
          ggtheme = theme_gray() # Change theme
)


fviz_dend(res.hc, cex = 0.6, k = 4,
          k_colors = "ucscgb", type = "circular")

library(igraph)

fviz_dend(res.hc,cex = 0.7, k = 4, k_colors = "jco",
          type = "phylogenic", repel = TRUE)


fviz_dend(res.hc, cex = 0.7,
          k = 4, k_colors = "npg",
          type = "phylogenic",
          repel = T)



