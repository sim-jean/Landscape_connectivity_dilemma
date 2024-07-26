---
title: "Final extension"
author: "Simon Jean"
date: "2024-04-30"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)

library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(landscapeR)
library(raster)
```

## Analysis of real looking landscapes in light of our results
Set landscape size : 

```{r, parameter and functions}
size = 200

viz_landscape = function(landscape){
 melt(landscape)%>%
    mutate(value = as.factor(value))%>%
    ggplot(aes(x= Var1, y = Var2, fill = value))+
    geom_tile()+
    scale_fill_manual(values= c('white', 'orange', 'red'))+
    theme_minimal()
}

next_landscape = function(landscape, treatment){
  vector_version = c(landscape)
  next_landscape = pmin((1+vector_version), rep(2,length(vector_version)))
  next_landscape = next_landscape * (1 - c(treatment))
  return(next_landscape)
}

```

Now originate landscape : 

```{r, }
landscapes = list()
nlandscape = 15
for(x in 1:nlandscape){
    m  = matrix(0, size, size)
    r  = raster(m, xmn = 0, xmx = size, ymn = 0, ymx = size)
    num1 = round(runif(1, 0, size))
    size1 = round(runif(num1, 0, size))
    rr = makeClass(r, num1, size1, val= 1)

    num2 = round(runif(1, 0, size))
    size2 = round(runif(num2, 0, size))
    rr = makeClass(rr, num2, size2, val = 2)
    plot(rr)
    landscapes = append(landscapes, list(rr@data@values))
}


```



```{r, own}
#fictitious_landscape = matrix(round(runif(size^2,0,2 )), nrow = size, ncol = size)

#viz_landscape(landscapes[1])

viz_landscape1 = matrix(landscapes[1], size, size)
viz_landscape(viz_landscape1)

treatment = matrix(round(runif(size^2, 0, 1)), nrow = size, ncol = size)

viz_landscape(matrix(next_landscape(fictitious_landscape, treatment),100,100))
```