library(RColorBrewer)
library(here)

proj.name <- "LCM_data_1mlrg_1m_20190729"
setwd(here("outputs", proj.name))
rules <- read.csv("rules.csv")
rules <- rules[order(rules$weight),]
rules$curve <- NULL

# add curves
rules$att <- paste(rules$a,rules$b, rules$scalar, sep = "-")
curves <- data.frame(att = unique(rules$att))
curves$curve <- paste("y",1:nrow(curves), sep = "")
rules <- merge(rules, curves, by = "att")

rules <- rules[!duplicated(rules$curve), ]
rules <- rules[order(rules$curve), ]
png("curves.png", width = 6, height = 4, units = "in", res = 600)

# dev.new(width=24,height=16,noRStudioGD = TRUE)
val <- function(x, a, b, scalar, weight) {
    (1/(1 + (exp(((x/scalar) - a) * b)))) * weight
}
plotcol <- brewer.pal(6, "Accent")
plot(1, type = "n", xlim = c(0, 2000), ylim = c(0, max(rules$weight)), xlab = "Distance", ylab = "Weight", bty = "n")
for (s in 1:nrow(rules)) {
    lines(val(1:2000, a = rules$a[s], b = rules$b[s], scalar = rules$scalar[s], weight = rules$weight[s]), type = "l", col = plotcol[s], lwd = 2)
}
legend("topright", legend = rules$curve, col = plotcol, lty = c(1, 1), cex = 0.8, bty = "n")

dev.off()
