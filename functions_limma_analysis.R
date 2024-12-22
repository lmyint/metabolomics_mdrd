#####################################################################
## Functions needed to conduct the main analysis
#####################################################################

#' Convert a numeric vector into percentiles
percentile <- function(x) {
    rank(x, ties.method = "min")*100/length(x)
}

#' Add QC column and its percentile version to row_data
#' @param x the QC variable
#' @param col_name the name of the QC variable
#' @return the same row_data supplemented with 2 additional columns
.add_qc_metrics <- function(row_data, x, col_name) {
    perc_col_name <- paste0("perc_", col_name)
    row_data[[col_name]] <- x
    row_data[[perc_col_name]] <- percentile(x)
    row_data
}

#' Add necessary variables to col_data and make some checks
#' @param col_data The colData object
#' @param primary_obj Either a `SummarizedExperiment` or an abundance
#'                    matrix. `col_data` will be subset to have the
#'                    same samples as in `primary_obj`.
#' @return A new version of `col_data`
prepare_col_data <- function(col_data, primary_obj) {
    ## Add a binary smoke_status variable
    col_data <- col_data %>%
        mutate(smoke_status = smoke_years > 0)
    
    if (is(primary_obj, "SummarizedExperiment")) {
        subject_ids_primary_obj <- colData(primary_obj)$SUBJECT_ID
        num_samples_primary_obj <- ncol(primary_obj)
    } else if (is(primary_obj, "matrix")) {
        if (nrow(primary_obj) < ncol(primary_obj)) {
            warning("More columns than rows. (Assuming columns=samples.)")
        }
        subject_ids_primary_obj <- colnames(primary_obj)
        num_samples_primary_obj <- ncol(primary_obj)
    }
    
    ## Make sure that all subject IDs from primary_obj
    ## are contained in supplied col_data
    stopifnot(all(subject_ids_primary_obj %in% col_data$SUBJECT_ID))
    
    ## Make sure that supplied col_data uses exactly 
    ## the subject IDs in primary_obj
    col_data_temp <- tibble(SUBJECT_ID = subject_ids_primary_obj) %>%
        left_join(col_data)
    
    stopifnot(nrow(col_data_temp)==num_samples_primary_obj)
    
    col_data_temp
}

run_analysis <- function(col_data, se, symptom_var_type = c("quant", "cat"), metab_type = c("keep_NA", "QRILC"), adj_formula = c("default", "no_uun"), only_xeno = FALSE) {
    symptom_var_type <- match.arg(symptom_var_type)
    metab_type <- match.arg(metab_type)
    adj_formula <- match.arg(adj_formula)
    
    ## Create binary smoke status variable and align
    ## supplied `col_data` with subject IDs in colData(se)
    col_data <- prepare_col_data(col_data, se)
    
    ## Add metabolite IDs to row data to facilitate joining with topTable
    row_data <- rowData(se)
    row_data$metab_id <- rownames(row_data)
    rowData(se) <- row_data
    
    ## Regardless of metab_type, QC info is for non-imputed data
    ## Add QC info about each metabolite to give context for top hits
    ## % missing, mean, SD, CV AND percentile versions
    log_abund <- assay(se, "log_abund") ## Non-imputed data
    mean_abund <- rowMeans(log_abund, na.rm = TRUE)
    sd_abund <- rowSds(log_abund, na.rm = TRUE)
    cv_abund <- sd_abund/mean_abund
    perc_missing <- rowMeans(is.na(log_abund))
    
    row_data <- .add_qc_metrics(row_data, mean_abund, "mean_abund")
    row_data <- .add_qc_metrics(row_data, sd_abund, "sd_abund")
    row_data <- .add_qc_metrics(row_data, cv_abund, "cv_abund")
    row_data <- .add_qc_metrics(row_data, perc_missing, "perc_missing")
    rowData(se) <- row_data
    
    ## Extract desired metabolomics data for modeling
    if (metab_type=="keep_NA") {
        log_abund <- assay(se, "log_abund")
    } else if (metab_type=="QRILC") {
        log_abund <- assay(se, "log_abund_qrilc")
    }
    
    old_opt <- options()$na.action
    options(na.action = "na.pass")
    
    ## Loop over symptom numbers
    symptom_df <- get_symptom_definitions()
    
    results <- lapply(seq_len(nrow(symptom_df)), function(i) {
        ## Get symptom score
        col_data$symptom_score <- get_symptom_score(symptom_df$num[i], col_data, symptom_var_type = symptom_var_type)
        
        ## Create design matrix
        if (adj_formula=="default") {
            model_matrix_formula <- ~symptom_score + age + female + racecat + diab + cause_ckd_ckdbc + smoke_status + sys + upro + uun + gfr + bmi + study + bp + diet
        } else if (adj_formula=="no_uun") {
            ## Remove uun from adjustment set
            model_matrix_formula <- ~symptom_score + age + female + racecat + diab + cause_ckd_ckdbc + smoke_status + sys + upro + gfr + bmi + study + bp + diet
        }
        
        design <- model.matrix(model_matrix_formula, data = col_data)
        colnames(design) <- colnames(design) %>%
            str_replace_all(" ", "_")
        colnames(design)[str_detect(colnames(design), "Intercept")] <- "intercept"
        
        ## Remove samples with NAs in `design`
        samples_keep <- complete.cases(design)
        design <- design[samples_keep,,drop=FALSE]
        log_abund <- log_abund[, samples_keep, drop=FALSE]
        
        if (sum(samples_keep)==0) {
            stop("No samples remain")
        }
        
        ## Run limma
        fit <- lmFit(log_abund, design = design)
        fit <- eBayes(fit)
        
        if (symptom_var_type=="quant") {
            tt <- topTable(fit, coef = "symptom_score", number = Inf)
        } else if (symptom_var_type=="cat") {
            select_coeffs <- c(
                "symptom_scorecat1_less",
                "symptom_scorecat1_most",
                "symptom_scorecat23_less",
                "symptom_scorecat23_most"
            )
            tt <- topTable(fit, coef = select_coeffs, number = Inf)
        }
        res <- summarize_results(tt, se, only_xeno = only_xeno)
        attr(res, "samples_removed") <- sum(!samples_keep)
        res
        
    })
    names(results) <- symptom_df$descrip
    
    ## Reset NA options
    options(na.action = old_opt)
    
    results
}


#####################################################################
## Functions for summarizing and displaying results
#####################################################################

#' Summarize results of differential analysis
#' ASSUMES THAT METABS AND SAMPLES IN tt AND se MATCH UP
#'
#' @param tt Result of topTable from `limma`
#' @param se SummarizedExperiment object
#' @param only_xeno Only keep Xenobiotics in the top table? (Default: FALSE)
#' @return A topTable augmented with BH-adjusted p-values and row data,
#'         with metabolites that are not of interest filtered out
summarize_results <- function(tt, se, only_xeno = FALSE) {
    stopifnot(nrow(tt)==nrow(se))
    row_data <- rowData(se) %>% as.data.frame()
    
    ## Rownames of the topTable correspond to the row names of row data
    ## Explicitly add row names as a variable and join with row data
    tt$metab_id <- rownames(tt)
    tt <- tt %>%
        left_join(row_data)
    
    ## Filter out NAs for SUPER_PATHWAY: indicates unidentified peaks
    tt <- tt %>% filter(!is.na(SUPER_PATHWAY))
    if (!only_xeno) {
        ## Filter out Xenobiotics
        tt <- tt %>% filter(SUPER_PATHWAY != "Xenobiotics")
    } else {
        ## Filter to only include Xenobiotics
        tt <- tt %>% filter(SUPER_PATHWAY == "Xenobiotics")
    }
    
    ## Obtain FDR corrected p-values
    tt$adj_pval <- p.adjust(tt$P.Value, method = "BH")
    tt <- tt %>% arrange(adj_pval)
    
    ## Add a column giving the percentile of logFC (absolute value)
    if ("logFC" %in% colnames(tt)) {
        tt$perc_logFC <- percentile(abs(tt$logFC))
    }
    
    tt
}

kable_scroll_print <- function(dat, caption, row_height_px = 75) {
    height_px <- 100+row_height_px*nrow(dat)
    knitr::kable(dat, caption = caption) %>%
        kable_styling() %>%
        scroll_box(width = "100%", height = paste0(height_px, "px")) %>%
        print()
    cat("\n")
}

#' Display results and write results to file
#' @param results A list of `topTable`'s from `runAnalysis()`
#' @param file_name_descrip String. A descriptor prepended to 
#'                          individual symptom file names.
display_write_results <- function(results, file_name_descrip) {
    symptom_df <- get_symptom_definitions()
    for (i in seq_along(results)) {
        tt <- results[[i]]
        sub_tt <- tt %>%
            filter(adj_pval < 0.1) %>%
            select(P.Value, adj_pval, HMDB_ID, BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY, PLATFORM, logFC, perc_logFC, mean_abund:perc_perc_missing)
        write_csv(sub_tt, path = paste0("../results/", file_name_descrip, "_", symptom_df$alt_descrip[i], ".csv"))
        kable_scroll_print(sub_tt, caption = names(results)[i])
        cat("\n")
    }
}
