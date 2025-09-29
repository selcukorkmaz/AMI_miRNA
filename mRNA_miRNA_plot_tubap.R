# PNG cihazı (4 gen için 2x2 düzen)
png("mrna_barplots.png", width = 2000, height = 1600, res = 300)

par(mfrow = c(2, 2), mar = c(3, 6, 2, 1))

plot_gene <- function(name, mean_k, sd_k, mean_a, sd_a, pval, logfc=NULL) {
  means <- c(mean_k, mean_a)
  sds <- c(sd_k, sd_a)
  ylim_max <- max(means + sds) * 1.3
  bp <- barplot(
    means,
    ylim = c(0, ylim_max+0.5),
    names.arg = c("Control", "AMI"),
    ylab = name,
    col = c("gray80", "gray40"),
    font.lab = 2,
    font.axis = 2
  )
  arrows(bp, means - sds, bp, means + sds, angle = 90, code = 3, length = 0.1, lwd = 2)

  # p etiketi
  y_pos <- max(means + sds) * 1.15

  if(!is.null(logfc)){
  p_label <- if (pval < 0.001) {
    paste0("p < 0.001", " logFC = ", logfc)
  } else {
    paste0("p = ", formatC(pval, digits = 3, format = "f"), ", logFC = ", logfc)
  }
  }else{
    p_label <- if (pval < 0.001) {
      "p < 0.001"
    } else {
      paste0("p = ", formatC(pval, digits = 3, format = "f"))
    }

  }


  segments(bp[1], y_pos, bp[2], y_pos)
  segments(bp[1], y_pos, bp[1], y_pos - 0.1)
  segments(bp[2], y_pos, bp[2], y_pos - 0.1)
  text(x = mean(bp), y = y_pos + 0.3, labels = p_label, cex = 0.8, font = 1)
}

# 4 gen
plot_gene("Fbn1 mRNA Expression\n(Fold Change vs. Control)",     1.014, 0.097, 1.489, 0.477, 0.366, 0.554)
plot_gene("Hamp mRNA Expression\n(Fold Change vs. Control)",     1.247, 0.500, 2.915, 1.047, 0.201, 1.225)
plot_gene("Col3a1 mRNA Expression\n(Fold Change vs. Control)",   1.012, 0.084, 1.804, 0.234, 0.019, 0.834)
plot_gene("Cyp2e1 mRNA Expression\n(Fold Change vs. Control)",   1.663, 0.849, 0.114, 0.069, 0.043,-3.867)

dev.off()

# ---- miRNA barplotları tek panelde ----


png("miRNA_barplots.png", width = 1750, height = 800, res = 300)

par(mfrow = c(1, 3), mar = c(3, 6, 2, 1))

plot_gene("miR-1b Expression\n(Fold Change vs. Control)",     1.011, 0.081, 0.501, 0.043, 0.001, -1.013)
plot_gene("miR-214-3p Expression\n(Fold Change vs. Control)", 1.046, 0.186, 2.388, 0.114, 0.0000001, 1.191)
plot_gene("miR-133b-3p Expression\n(Fold Change vs. Control)", 1.028, 0.131, 0.352, 0.074, 0.004, -1.546)

dev.off()


png("troponin.png", width = 1750, height = 1100, res = 300)

par(mfrow = c(1, 2), mar = c(3, 6, 2, 2))

plot_gene("Creatine Kinase-MB (ng/ml)",     5.093, 1.403, 9.261, 0.413, 0.00001)
plot_gene("Cardiac Troponin-I (ng/ml)", 1.420, 0.454, 2.933, 0.119, 0.0000001)

dev.off()
