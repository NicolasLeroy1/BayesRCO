# Generate complex overlapping annotations
set.seed(42)
n_snps <- 2854
n_cats <- 5

# Create matrix
mat <- matrix(0, nrow=n_snps, ncol=n_cats)

for(i in 1:n_snps) {
    # Ensure at least one category, allow up to 3
    n_assign <- sample(1:3, 1, prob=c(0.6, 0.3, 0.1))
    cats <- sample(1:n_cats, n_assign)
    mat[i, cats] <- 1
}

write.table(mat, "toy_dataset/complex_annot.txt", row.names=FALSE, col.names=FALSE, sep=" ")
cat("Generated toy_dataset/complex_annot.txt with", n_snps, "SNPs and", n_cats, "categories.\n")
