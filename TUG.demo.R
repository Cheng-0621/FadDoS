if (!require("plotly")) install.packages("plotly")
if (!require("paletteer")) install.packages("paletteer")
library(plotly)
library(paletteer)

dat <- readRDS("../FadDoS/GitHub/TUG.demo.rds")

fig <- dat %>%
  plot_ly(
    x = ~x,
    y = ~y,
    z = ~z,
    frame = ~time,
    color = ~index,
    colors = c(paletteer_c("ggthemes::Red-Green-Gold Diverging", 25)),
    type = 'scatter3d',
    mode = 'markers',
    marker = list(size = 8)
  )  %>% 
  layout(font=list(size=15), scene=list(xaxis = list(tickfont =list(size=15)), yaxis = list(tickfont =list(size=15)), zaxis = list(tickfont =list(size=15))))
fig
