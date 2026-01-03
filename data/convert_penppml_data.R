# =============================================================================
# Convert R penppml test data to Stata format
# =============================================================================
# Run this script in R to create the Stata test datasets
# Requires: haven package (install.packages("haven") if needed)
#
# Usage:
#   1. Open R or RStudio
#   2. Set working directory to this folder
#   3. Run: source("convert_penppml_data.R")
# =============================================================================

library(haven)

# Path to the penppml R package data
penppml_data_path <- "../../../penppml_r/data"

# Load the RDA files
load(file.path(penppml_data_path, "trade.rda"))
load(file.path(penppml_data_path, "countries.rda"))

cat("Loaded trade data:", nrow(trade), "rows,", ncol(trade), "columns\n")
cat("Loaded countries data:", nrow(countries), "rows,", ncol(countries), "columns\n")

# Convert factors to character/numeric for Stata compatibility
trade$exp <- as.character(trade$exp)
trade$imp <- as.character(trade$imp)
trade$time <- as.numeric(as.character(trade$time))
trade$agreement <- as.character(trade$agreement)

countries$iso <- as.character(countries$iso)
countries$name <- as.character(countries$name)
countries$region <- as.character(countries$region)
countries$subregion <- as.character(countries$subregion)

# Save full datasets
write_dta(trade, "trade.dta", version = 14)
write_dta(countries, "countries.dta", version = 14)

cat("\nSaved: trade.dta\n")
cat("Saved: countries.dta\n")

# Create Americas subset (smaller, for quick tests)
americas_iso <- countries$iso[countries$region == "Americas"]
trade_americas <- trade[(trade$exp %in% americas_iso) & (trade$imp %in% americas_iso), ]
write_dta(trade_americas, "trade_americas.dta", version = 14)

cat("Saved: trade_americas.dta (", nrow(trade_americas), " rows)\n")

# Print variable names for reference
cat("\n=== Trade dataset variables ===\n")
cat(paste(names(trade), collapse = "\n"))

cat("\n\n=== Countries dataset variables ===\n")
cat(paste(names(countries), collapse = "\n"))

cat("\n\nData conversion complete!\n")
