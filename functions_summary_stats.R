#####################################################################
## Functions for providing summaries of quantitative
## and categorical variables
#####################################################################

#' Display summary statistics for a quantitative variable
#' @param x The quantitative variable
#' @param include_miss If `TRUE`, include the number and fraction
#'                     of missing values.
#' @param plot If `TRUE`, show a density plot.
summarize_quant_var <- function(x, include_miss=TRUE, plot=FALSE) {
    cat("Range:", min(x, na.rm = TRUE), "-", max(x, na.rm = TRUE), "\n")
    cat("IQR:", quantile(x, probs = 0.25, na.rm = TRUE), "-", quantile(x, probs = 0.75, na.rm = TRUE), "\n")
    cat("Mean:", mean(x, na.rm = TRUE), "\n")
    cat("Median:", median(x, na.rm = TRUE), "\n")
    cat("SD:", sd(x, na.rm = TRUE), "\n")
    if (include_miss) {
        cat("Number missing:", sum(is.na(x)), "\t Fraction:", mean(is.na(x)))
    }
    if (plot) {
        plot(density(x, na.rm = TRUE), xlab = "", main = "")
    }
}

#' Display summary statistics for a categorical variable
#' @param x The categorical variable
summarize_categ_var <- function(x) {
    tab <- table(x, useNA = "always")
    names(dimnames(tab)) <- NULL
    prop_tab <- prop.table(tab)
    print(tab)
    print(prop_tab)
}

#' Display summary statistics for a dataset
#' @param dat The dataset (typically a form of colData)
#' @param plot If `TRUE`, show density plots.
display_summary_stats_baseline <- function(dat, plot=FALSE) {
    dat <- as.data.frame(dat)
    variables <- c("age", "female", "racecat", "diab", "renaldx", "cause_ckd_ckdbc", "smoke_status", "sys", "upro", "uun", "gfr", "bmi", "study", "bp", "diet", "cad", "pepulc", "cancer", "cerebrovasc", "pvd", "hyperlip", "hyperten", "alb", "phos", "hb")
    variable_types <- c("q", rep("c", 6), rep("q", 5), rep("c", 10), rep("q", 3))
    
    cat("Sample size:", nrow(dat), "\n")
    
    for (i in seq_along(variables)) {
        var_name <- variables[i]
        var_type <- variable_types[i]
        v <- dat[[var_name]]
        
        cat(var_name, paste(rep("-", 30), collapse = ""), "\n")
        
        if (var_name=="smoke_status") {
            v <- ifelse(dat$smoke_years > 0, "ever_smoker", "never_smoker")
        }
        if (var_name=="renaldx") {
            cat("renaldx among those with cause_ckd_bc = PKD\n")
            v <- dat %>% dplyr::filter(cause_ckd_ckdbc=="PKD") %>% pull(renaldx)
            summarize_categ_var(v)
            cat("\nrenaldx among those with cause_ckd_bc = GN\n")
            v <- dat %>% dplyr::filter(cause_ckd_ckdbc=="GN") %>% pull(renaldx)
            summarize_categ_var(v)
            cat("\nrenaldx among those with cause_ckd_bc = Others\n")
            v <- dat %>% dplyr::filter(cause_ckd_ckdbc=="Others") %>% pull(renaldx)
            summarize_categ_var(v)
            cat("\nOverall summary of renaldx variable\n")
        }
        
        if (var_type=="q") {
            summarize_quant_var(v, plot = plot)
        } else if (var_type=="c") {
            summarize_categ_var(v)
        }
        cat("\n\n")
    }
}

#####################################################################
## Functions for summarizing symptom information
#####################################################################

#' Get a `tibble`` of symptom information
#' @return A `tibble` containing symptom number, description
#'         grouping, and filename prefix
get_symptom_definitions <- function() {
    tibble(
        num = c(5:8, 14, 17, 18, 20, 23, 25, 26),
        descrip = c("Bad Taste", "Loss of Appetite", "Nausea", "Vomiting", "Itching", "Lack of Pep and Energy", "Tiring Easily/Weakness", "Numbness/Tingling of Hands/Feet", "Fall Asleep During Day", "Decreased Alertness", "Forgetfulness"),
        group = c(rep(1,4), 2, 4, 4, 3, 4, 4, 4),
        alt_descrip = c("05_bad_taste", "06_loss_appetite", "07_nausea", "08_vomiting", "14_itching", "17_lack_pep_energy", "18_tiring_easily_weakness", "20_numb_tingle_hand_feet", "23_fall_asleep_day", "25_decreased_alertness", "26_forgetfulness")
    )
}

#' Extract the days and severity column
#' @param symptom_num An integer giving the symptom number
#' @return A `tibble` with the extracted days and severity
get_symptom_days_severity <- function(symptom_num, col_data) {
    ## Extract days (col A) and severity (col B)
    num <- sprintf("%.2d", symptom_num)
    A_col <- paste0("F26Q", num, "A")
    B_col <- paste0("F26Q", num, "B")
    A <- col_data[[A_col]]
    B <- col_data[[B_col]] %>% as.character() %>% as.integer()
    tibble(days = A, severity = B)
}

#' Compute score for a given symptom
#' @param symptom_num An integer giving the symptom number
#' @param col_data The colData to be used
#' @param symptom_var_type Either `"quant"` for quantitative score
#'                         or `"cat"` for categorical "score"
#' @return A vector containing the Q or C symptom score
get_symptom_score <- function(symptom_num, col_data, symptom_var_type = c("quant", "cat")) {
    symptom_var_type <- match.arg(symptom_var_type)
    
    ## Obtain the day and severity columns for this symptom
    days_sev <- get_symptom_days_severity(symptom_num, col_data)
    
    if (symptom_var_type=="quant") {
        ## Create symptom scores as product: days*severity (A*B)
        symp <- days_sev$days * days_sev$severity
    } else if (symptom_var_type=="cat") {
        ## Define new categories to see sample sizes
        days_sev <- days_sev %>%
            mutate(symp_cat = case_when(
                severity==1 & days <= 15 ~ "cat1_less",
                severity==1 & days > 15 ~ "cat1_most",
                severity %in% c(2,3) & days <= 15 ~ "cat23_less",
                severity %in% c(2,3) & days > 15 ~ "cat23_most",
                days==0 ~ "cat0"
            ))
        symp <- days_sev$symp_cat
    }
    
    symp
}

#' Display summary statistics for symptoms for a dataset
#' @param dat The dataset (typically a form of colData)
display_summary_stats_symptoms <- function(dat) {
    symptom_df <- get_symptom_definitions()
    for (i in seq_len(nrow(symptom_df))) {
        cat(symptom_df$descrip[i], paste(rep("-", 20), collapse = ""), "\n")
        symp_score <- get_symptom_score(symptom_df$num[i], dat)
        
        ## Missing values
        cat("Number missing:", sum(is.na(symp_score)), "Fraction:", mean(is.na(symp_score)), "\n")
        
        ## Asymptomatic vs symptomatic
        cat("Number symptomatic (nonzero score):", sum(symp_score!=0, na.rm = TRUE), "\t Fraction:", mean(symp_score!=0, na.rm = TRUE), "\n")
        
        ## Among the symptomatic, get statistics
        cat("Among symptompatic:\n")
        symp_score_subs <- symp_score[!is.na(symp_score)]
        symp_score_subs <- symp_score_subs[symp_score_subs!=0]
        summarize_quant_var(symp_score_subs, include_miss = FALSE)
        cat("\n\n")
    }
}

#' Get summary statistics for symptoms for a dataset
#' @param dat The dataset (typically a form of colData)
get_summary_stats_symptoms <- function(dat) {
    dat <- as.data.frame(dat)
    symptom_df <- get_symptom_definitions()
    
    study_arms <- expand.grid(
        bp = c("Low", "Moderate"),
        diet = c("Very low", "Low", "Usual")
    )
    lapply(seq_len(nrow(symptom_df)), function(i) {
        lapply(seq_len(nrow(study_arms)), function(j) {
            study_arm_bp <- study_arms$bp[j]
            study_arm_diet <- study_arms$diet[j]
            dat_subs <- dat %>%
                dplyr::filter(bp==study_arm_bp, diet==study_arm_diet)
          
            ## Get symptom score vector
            symp_score <- get_symptom_score(symptom_df$num[i], dat_subs)
            
            ## Missing values
            num_missing <- sum(is.na(symp_score))
            frac_missing <- mean(is.na(symp_score))
            
            ## Asymptomatic vs symptomatic (nonzero score)
            num_symptomatic <- sum(symp_score!=0, na.rm = TRUE)
            frac_symptomatic <- mean(symp_score!=0, na.rm = TRUE)
            
            ## Among the symptomatic, get statistics
            symp_score_subs <- symp_score[!is.na(symp_score)]
            symp_score_subs <- symp_score_subs[symp_score_subs!=0]
            iqr_score <- stats::IQR(symp_score_subs, na.rm = TRUE)
            med_score <- median(symp_score_subs, na.rm = TRUE)
            
            tibble(
                symptom_descrip = symptom_df$descrip[i],
                study_bp = study_arm_bp,
                study_diet = study_arm_diet,
                num_symptomatic, frac_symptomatic,
                p25_score = iqr_score[1],
                med_score,
                p75_score = iqr_score[2]
            )
        }) %>% bind_rows()
    }) %>% bind_rows()
}

# Inference for symptom comparisons between groups
#' @param dat1 Dataset for first comparison group (colData)
#' @param dat2 Dataset for second comparison group (colData)
compare_symptomatic <- function(dat1, dat2) {
    symptom_df <- get_symptom_definitions()
    for (i in seq_len(nrow(symptom_df))) {
        cat(symptom_df$descrip[i], paste(rep("-", 20), collapse = ""), "\n")
        symp_score1 <- get_symptom_score(symptom_df$num[i], dat1)
        symp_score2 <- get_symptom_score(symptom_df$num[i], dat2)
        
        num_symptomatic1 <- sum(symp_score1!=0, na.rm = TRUE)
        num_symptomatic2 <- sum(symp_score2!=0, na.rm = TRUE)
        
        n1 <- sum(!is.na(symp_score1))
        n2 <- sum(!is.na(symp_score2))
        
        pt_result <- prop.test(
            x = c(num_symptomatic1, num_symptomatic2),
            n = c(n1, n2)
        )
        print(pt_result)
    }
}
