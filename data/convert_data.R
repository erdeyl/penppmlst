# Script to convert R penppml data to Stata format
# Run this in R with the penppml and haven packages installed

library(penppml)
library(haven)

# Load trade data
data(trade)
data(countries)

# Convert to Stata format
write_dta(trade, "trade.dta", version = 14)
write_dta(countries, "countries.dta", version = 14)

# Print summary for verification
cat("Trade data:\n")
cat("  Rows:", nrow(trade), "\n")
cat("  Cols:", ncol(trade), "\n")
cat("  Variables:", paste(names(trade), collapse = ", "), "\n\n")

cat("Countries data:\n")
cat("  Rows:", nrow(countries), "\n")
cat("  Cols:", ncol(countries), "\n")
cat("  Variables:", paste(names(countries), collapse = ", "), "\n")
