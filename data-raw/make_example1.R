# Generate the example1 dataset for the package
# Run this script from the package root directory

# Source the sim_func1 function
source("R/data.R")

example1 <- sim_func1(n = 1000, seed = 123)

save(example1, file = "data/example1.rda", compress = "xz")

cat("Created data/example1.rda with", nrow(example1), "rows\n")
