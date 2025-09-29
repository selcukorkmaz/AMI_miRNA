png("troponin_ckmb_plot.png", width = 1800, height = 1200, res = 300)

# Plot layout: 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))  # Set margins

#### A. Troponin I ####
grup <- c(rep("Kontrol", 8), rep("AMI", 7))
konsantrasyon <- c(
  2.284910965, 1.380506092, 1.352389878, 1.272727273,
  0.724461106, 1.183692596, 1.389878163, 1.774133083,
  2.997188379, 2.844111215, 3.072164948, 2.90659169,
  3.072164948, 2.75663855, 2.880037488
)
df <- data.frame(Grup = grup, Konsantrasyon = konsantrasyon)
means <- tapply(df$Konsantrasyon, df$Grup, mean)
ses <- tapply(df$Konsantrasyon, df$Grup, function(x) sd(x)/sqrt(length(x)))
pval <- t.test(Konsantrasyon ~ Grup, data = df)$p.value

bp <- barplot(
  means,
  beside = TRUE,
  ylim = c(0, 4.5),
  col = c(rgb(1, 0, 0, 0.2), rgb(0, 0, 1, 0.2)),
  border = c("red", "blue"),
  ylab = "Troponin-I (ng/ml)",
  names.arg = names(means),
  font.lab = 2,
  font.axis = 2
)
arrows(x0 = bp, y0 = means - ses, x1 = bp, y1 = means + ses,
       angle = 90, code = 3, length = 0.1, lwd = 2, col = c("red", "blue"))
y_max <- max(means + ses) + 0.4
segments(bp[1], y_max, bp[2], y_max)
segments(bp[1], y_max, bp[1], y_max - 0.2)
segments(bp[2], y_max, bp[2], y_max - 0.2)
text(x = mean(bp), y = y_max + 0.2, labels = "p < 0,001", cex = 1.2, font = 2)

# Add 'A' to upper-left
usr <- par("usr")
text(x = usr[1] + 0.2, y = usr[4] - 0.2, labels = "A", font = 1, cex = 1.5)

#### B. CK/MB ####
grup <- c(rep("Kontrol", 8), rep("AMI", 7))
konsantrasyon <- c(
  4.616780045, 5.603174603, 4.934240363, 5.693877551,
  5.149659864, 5.274376417, 4.854875283, 4.616780045,
  9.397581255, 8.354497354, 8.0143613, 8.490551776,
  10.98488284, 8.16553288, 11.41950113
)
df <- data.frame(Grup = grup, Konsantrasyon = konsantrasyon)
means <- tapply(df$Konsantrasyon, df$Grup, mean)
ses <- tapply(df$Konsantrasyon, df$Grup, function(x) sd(x)/sqrt(length(x)))
pval <- t.test(Konsantrasyon ~ Grup, data = df)$p.value

bp <- barplot(
  means,
  beside = TRUE,
  ylim = c(0, 12.5),
  col = c(rgb(1, 0, 0, 0.2), rgb(0, 0, 1, 0.2)),
  border = c("red", "blue"),
  ylab = "CK/MB (ng/ml)",
  names.arg = names(means),
  font.lab = 2,
  font.axis = 2
)
arrows(x0 = bp, y0 = means - ses, x1 = bp, y1 = means + ses,
       angle = 90, code = 3, length = 0.1, lwd = 2, col = c("red", "blue"))
y_max <- max(means + ses) + 0.4
segments(bp[1], y_max, bp[2], y_max)
segments(bp[1], y_max, bp[1], y_max - 0.2)
segments(bp[2], y_max, bp[2], y_max - 0.2)
text(x = mean(bp), y = y_max + 0.4, labels = "p < 0,001", cex = 1.2, font = 2)

# Add 'B' to upper-left
usr <- par("usr")
text(x = usr[1] + 0.2, y = usr[4] - 0.48, labels = "B", font = 1, cex = 1.5)

dev.off()

