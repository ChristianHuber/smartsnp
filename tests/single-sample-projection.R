library(smartsnp)

set.seed(17)
genotypes <- matrix(sample(0:2, 80 * 20, replace = TRUE), nrow = 80)
genotypes[1, ] <- 0
genotypes[2, 3] <- 9
input <- tempfile()
write.table(genotypes, input, row.names = FALSE, col.names = FALSE, quote = FALSE)
groups <- rep(c("A", "B"), each = 10)

cases <- list(
  list(sample_project = 20),
  list(sample_project = 20, snp_remove = 4),
  list(sample_project = 20, missing_impute = "remove"),
  list(sample_project = c(19, 20))
)

for (name in c("smart_pca", "smart_mva")) {
  for (options in cases) {
    set.seed(42)
    arguments <- c(list(snp_data = input, sample_group = groups), options)
    if (name == "smart_mva") {
      arguments <- c(arguments, list(permanova = FALSE, permdisp = FALSE))
    }
    result <- do.call(getExportedValue("smartsnp", name), arguments)
    if (name == "smart_mva") result <- result$pca
    coordinates <- result$pca.sample_coordinates
    stopifnot(nrow(coordinates) == 20,
              all(is.finite(as.matrix(coordinates[, c("PC1", "PC2")]))))
  }
}
unlink(input)
