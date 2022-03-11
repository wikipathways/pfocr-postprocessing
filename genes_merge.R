## Merges new PFOCR gene results with prior results

## UPDATE SCRIPT TO NEW BATCH DATE
# Find/replace: 20210515

## SET WORKING DIR to source file location after copying into new batch folder
## e.g., Dropbox (Gladstone)/Pathway Figure OCR/20210515


## Read in new genes rds
genes_new <- readRDS("../pfocr_genes_20210515.rds")

## Read in old genes rds
genes_old <- readRDS("~/Dropbox (Gladstone)/PFOCR_25Years/pfocr_genes.rds")

## Restore original column names (if needed)
genes_new_final <- genes_new %>%
  select(-pfocr_year) 

names(genes_new_final) <- names(genes_old)
saveRDS(genes_new_final, "../pfocr_genes_draft.rds")
