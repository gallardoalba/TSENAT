# Quick grammar/syntax scanner
files <- c("TSENAT.Rmd", "TSENAT_appendix_A.Rmd", "TSENAT_appendix_B.Rmd")

possible_issues <- list()

for (file in files) {
  content <- readLines(file)
  issues <- data.frame(File = character(), Line = integer(), Issue = character(), Context = character(), stringsAsFactors = FALSE)
  
  for (i in seq_along(content)) {
    line <- content[i]
    
    # Check for common issues
    # 1. Inconsistent spacing around punctuation
    if (grepl("\\s+\\.", line)) {
      issues <- rbind(issues, data.frame(File = file, Line = i, Issue = "Space before period", Context = substr(line, 1, 60)))
    }
    
    # 2. Multiple spaces (except in code)
    if (!grepl("```", line) && grepl("  {2,}", line) && !grepl("^    ", line)) {
      issues <- rbind(issues, data.frame(File = file, Line = i, Issue = "Multiple spaces", Context = substr(line, 1, 60)))
    }
    
    # 3. Check for common verb agreement issues (simplified)
    if (grepl("\\bData\\s+has\\b|\\bData\\s+have\\b", line, ignore.case = TRUE)) {
      if (grepl("\\bData\\s+has\\b", line, ignore.case = TRUE)) {
        issues <- rbind(issues, data.frame(File = file, Line = i, Issue = "Data should use 'have' not 'has'", Context = substr(line, 1, 60)))
      }
    }
    
    # 4. Check for "utilize" vs "use" (style preference)
    if (grepl("\\butilize\\b", line, ignore.case = TRUE)) {
      issues <- rbind(issues, data.frame(File = file, Line = i, Issue = "'Use' preferred over 'utilize'", Context = substr(line, 1, 60)))
    }
  }
  
  if (nrow(issues) > 0) {
    possible_issues[[file]] <- issues
  }
}

# Print results
for (file in names(possible_issues)) {
  cat("\n=== ", file, " ===\n", sep = "")
  print(possible_issues[[file]])
}

if (length(possible_issues) == 0) {
  cat("\nNo major grammar/syntax issues detected in automated scan.\n")
}
