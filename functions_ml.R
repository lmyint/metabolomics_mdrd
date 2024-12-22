#####################################################################
## Functions for fitting individual models
#####################################################################

train_lasso <- function(data) {
    lambdas <- 10^seq(-3, 3, length.out = 100)
    lasso_mod <- train(
        symptom_score_res ~ .,
        data = data,
        method = "glmnet",
        trControl = trainControl(method = "cv", number = 10),
        tuneGrid = data.frame(alpha = 1, lambda = lambdas),
        metric = "RMSE"
    )
    lasso_mod
}

train_rand_forest <- function(data) {
    mtry_values <- seq(2, sqrt(nrow(data)), length.out = 5) %>% ceiling()
    rf_mod <- train(
        symptom_score_res ~ .,
        data = data,
        method = "rf",
        metric = "RMSE",
        trControl = trainControl(method = "oob"),
        tuneGrid = data.frame(mtry = mtry_values)
    )
    rf_mod
}

#####################################################################
## Functions for obtaining variable importance measures
#####################################################################

var_importance <- function(caret_mod, method) {
    mod <- caret_mod$finalModel
    if (method=="glmnet") {
        best_lambda <- mod$tuneValue$lambda
        coeffs_best <- tail(as.numeric(coef(mod, s = best_lambda)), -1)
        coef_mat <- mod$beta
        # first_zeros <- sapply(seq_len(nrow(coef_mat)), function(i) {
        #     is_zero <- coef_mat[i,]==0
        #     match(TRUE, is_zero)
        # })
        # first_zeros[is.na(first_zeros)] <- max(first_zeros, na.rm = TRUE)+1
        imp <- rowSums(coef_mat!=0)
        var_imp <- tibble(metab_id = rownames(coef_mat), importance_lasso = imp) %>%
            arrange(desc(importance_lasso))
    } else if (method=="rf") {
        var_imp <- randomForest::importance(mod, type = 2)
        var_imp <- tibble(
            metab_id = rownames(var_imp),
            importance_rf = var_imp[,"IncNodePurity"]
        ) %>%
            arrange(desc(importance_rf))
    }
    var_imp
}

#####################################################################
## Functions for fitting all models
#####################################################################

fit_ml_models <- function(col_data, se, log_abund = NULL, which_symptoms = NULL) {
    ## Store the appropriate log_abund
    if (is.null(log_abund)) {
        col_data <- prepare_col_data(col_data, primary_obj = se)
        
        log_abund <- assay(se, "log_abund")
        log_abund <- log_abund %>% t() %>% as_tibble()
    } else {
        col_data <- prepare_col_data(col_data, primary_obj = log_abund)
        log_abund <- log_abund %>% t() %>% as_tibble()
    }
    
    ## Keep only identified peaks
    ### First obtain them from se and rowData
    row_data <- rowData(se)
    bool_se_keep <- !is.na(row_data$SUPER_PATHWAY) & row_data$SUPER_PATHWAY != "Xenobiotics"
    met_ids_keep <- rownames(row_data)[bool_se_keep]
    ### Subset columns of log_abund
    bool_log_abund_keep <- colnames(log_abund) %in% met_ids_keep
    log_abund <- log_abund[,bool_log_abund_keep, drop=FALSE]
    
    ## Get symptom information
    symptom_df <- get_symptom_definitions()
    
    ## If no specific symptoms specified, choose all of them
    if (is.null(which_symptoms)) {
        which_symptoms <- symptom_df$num
    }
    
    ## Subset to chosen symptoms
    symptom_df <- symptom_df %>%
        filter(num %in% which_symptoms)
    
    ## Loop over symptoms
    results <- lapply(seq_len(nrow(symptom_df)), function(i) {
        ## Get symptom score
        col_data$symptom_score <- get_symptom_score(symptom_df$num[i], col_data)
        
        ## Fit model: symptom_score ~ lab+demographics
        lm_mod <- lm(symptom_score ~ age + female + racecat + diab + cause_ckd_ckdbc + smoke_status + sys + upro + uun + gfr + bmi + study + bp + diet, data = col_data)
        res_lm_mod <- residuals(lm_mod)
        stopifnot(length(res_lm_mod)==ncol(se))
        
        ## Use the residuals as the response variable in ML models
        ## Combine residuals and metabolite abundances in a dataset
        dat_ml <- bind_cols(tibble(symptom_score_res = res_lm_mod), log_abund)
        
        ## LASSO
        lasso_mod <- train_lasso(dat_ml)
        
        ## Random forest
        rf_mod <- train_rand_forest(dat_ml)
        
        list(lasso_mod, rf_mod)
    })
    
    results
}
