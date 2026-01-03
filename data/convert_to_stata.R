# Convert penppml R data to Stata format
library(haven)

# Load the RDA files
load("C:/Users/erdey/Claude/penppml_r/data/trade.rda")
load("C:/Users/erdey/Claude/penppml_r/data/countries.rda")

# Convert trade data - handle factors carefully
trade_df <- as.data.frame(trade)
trade_df$exp <- as.character(trade_df$exp)
trade_df$imp <- as.character(trade_df$imp)
trade_df$time <- as.numeric(as.character(trade_df$time))
trade_df$agreement <- as.character(trade_df$agreement)

# Clean variable names for Stata (replace dots and special chars with underscores)
names(trade_df) <- gsub("[^a-zA-Z0-9_]", "_", names(trade_df))

# Convert countries data
countries_df <- as.data.frame(countries)
for(col in names(countries_df)) {
  if(is.factor(countries_df[[col]])) {
    countries_df[[col]] <- as.character(countries_df[[col]])
  }
}
# Clean variable names for Stata
names(countries_df) <- gsub("[^a-zA-Z0-9_]", "_", names(countries_df))

# Write to Stata format
write_dta(trade_df, "C:/Users/erdey/Claude/penppmlst/data/trade.dta", version = 14)
write_dta(countries_df, "C:/Users/erdey/Claude/penppmlst/data/countries.dta", version = 14)

cat("Trade data:", nrow(trade_df), "rows,", ncol(trade_df), "columns\n")
cat("Countries data:", nrow(countries_df), "rows,", ncol(countries_df), "columns\n")
cat("Trade variables:", paste(names(trade_df), collapse = ", "), "\n")
