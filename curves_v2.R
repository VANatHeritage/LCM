library(RColorBrewer)
library(here)

proj.name <- "LCM_data_full_10m_20190829_v2"
setwd(here("outputs", proj.name))
rules <- read.csv("rules.csv")
rules <- rules[order(rules$category, rules$siteindex, decreasing = T),]
rules$curve <- NULL

# dev.new(width=24,height=16,noRStudioGD = TRUE)
val <- function(d, ds, si) {
  # modified power for site intensity; affects y-intecept
  y <- (-0.25*log(si))
  # distance decay function
  x <- (1/(1+exp(-2*((d-(ds*0.5))/(ds*0.25)))))
  return(x^y)
}

# assign line colors, styles
rules$col <- NA
cats <- table(rules$category)
rules$col[rules$category=="development"] <- brewer.pal(cats["development"], "Dark2")
rules$col[rules$category=="transport"] <- brewer.pal(cats["transport"], "Set2")
rules$col[rules$category=="landcover"] <- brewer.pal(cats["landcover"], "Set3")
rules$lty <- NA
rules$lty[rules$category=="development"] <- 3 # dotted
rules$lty[rules$category=="transport"] <- 2 # dashed
rules$lty[rules$category=="landcover"] <- 1 # solid

# create plot
png("curves.png", width = 8, height = 5, units = "in", res = 600)

nx <- max(rules$distthresh)*1.1
plot(1, type = "n", xlim = c(1, nx), ylim = c(0, 1), 
     xlab = "Distance (m)", ylab = "Site Impact", 
     bty = "n", main = "Distance Decay x Site Intensity")
for (s in 1:nrow(rules)) 
  lines(x = seq(1, nx, by = 5), y = val(seq(1, nx, by = 5), ds = rules$distthresh[s], si = rules$siteindex[s]), type = "l", col = rules$col[s], lty = rules$lty[s], lwd = 2.5)
legend(x = nx*0.8, y = .9, legend = rules$name, col = rules$col, lty = rules$lty, lwd = 2, cex = 0.7, bty = "n")

dev.off()
