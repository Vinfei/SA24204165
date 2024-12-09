## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
library(ggplot2)

## ----example1, echo=FALSE-----------------------------------------------------
data(mtcars)
p <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point() +
  labs(title = "汽车重量与每加仑英里数的关系",
       x = "重量 (1000磅)",
       y = "每加仑英里数")
print(p)

## ----example2, echo=FALSE-----------------------------------------------------
kable(head(mtcars), caption = "汽车数据集的前六行")

