#' Compare differential analysis results for all symptom groups
#' @param results List of `topTable`'s from `run_analysis()`
#' @return A list of results from `compare_symptom_group()`
compare_all_symptom_groups <- function(results) {
    ## Get symptom information
    symptom_df <- get_symptom_definitions()
    
    ## Obtain the groups that have multiple symptoms
    groups_multiple_metabs <- symptom_df$group[duplicated(symptom_df$group)] %>% unique()
    
    ## Loop over these multi-symptom groups
    res <- lapply(groups_multiple_metabs, function(g) {
        bool <- symptom_df$group==g
        ## List of topTables for this group
        tt_list_group <- results[bool]
        ## Symptom names for this group
        symp_names <- symptom_df$descrip[bool]
        compare_symptom_group(tt_list_group, symp_names)
    })
    names(res) <- paste("Group", groups_multiple_metabs)
    res
}

#' Compare symptoms within a single group
#' @param tt_list A list of `topTable`'s from `run_analysis()`
#' @param symptom_names Names of the symptoms in the group
#' @return Summary dataset of all differential 
#'         metabolites in the group.
compare_symptom_group <- function(tt_list, symptom_names) {
    ## Filter down to differential metabolites
    ## and get metabolite name
    sub_tt_list <- lapply(tt_list, function(tt) {
        tt_subs <- tt %>%
            filter(adj_pval < 0.1) %>%
            select(BIOCHEMICAL, SUPER_PATHWAY)
    })
    combined_data <- bind_rows(sub_tt_list)
    combined_data <- combined_data %>%
        arrange(SUPER_PATHWAY, BIOCHEMICAL)
    
    ## Get union of differential metabolites
    all_diff_metabs <- unique(combined_data$BIOCHEMICAL)
    
    ## For each metab in the union, get logFC and adj_pval
    sub_tt_list <- lapply(tt_list, function(tt) {
        tt %>%
            filter(BIOCHEMICAL %in% all_diff_metabs) %>%
            mutate(
                signif = adj_pval < 0.1,
                FC = 2^(logFC*30)
            ) %>%
            select(BIOCHEMICAL, FC, signif)
    })
    names(sub_tt_list) <- symptom_names
    combined_data <- bind_rows(sub_tt_list, .id = "symptom")
    combined_data
}

#' Plot results of compare_symptom_group
#' @param data The result of `compare_symptom_group()`
plot_symptom_comparison <- function(data) {
    ggplot(data, aes(x = FC, y = BIOCHEMICAL, color = symptom, shape = signif)) +
        geom_point(size = 3) +
        geom_vline(xintercept = 1) +
        scale_shape_manual(values = c(4,19), name = "Associated at\nFDR = 0.1") +
        scale_color_manual(values = c("deepskyblue", "darkorange", "darkorchid1", "limegreen", "deeppink"), name = "Symptom") +
        labs(x = "Adjusted fold change in metabolite\n per 30 unit change in symptom score", y = "") +
        theme_classic()
}


plot_symptom_comparison_both_sets <- function(data_v1, data_v2) {
    ## Identify the metabolites in common
    metabs_v1 <- unique(data_v1$BIOCHEMICAL)
    metabs_v2 <- unique(data_v2$BIOCHEMICAL)
    common_metabs <- intersect(metabs_v1, metabs_v2)
    
    ## Indicate common metabolites in data
    data_v1 <- data_v1 %>%
        mutate(common = BIOCHEMICAL %in% common_metabs) %>%
        mutate(
            BIOCHEMICAL = ifelse(common,
                                 paste0("ZZZ", BIOCHEMICAL),
                                 BIOCHEMICAL
                                 )
        )
    data_v2 <- data_v2 %>%
        mutate(common = BIOCHEMICAL %in% common_metabs) %>%
        mutate(
            BIOCHEMICAL = ifelse(common,
                                 paste0("ZZZ", BIOCHEMICAL),
                                 BIOCHEMICAL
            )
        )
    
    xlim <- range(data_v1$FC, data_v2$FC)
    plot_v1 <- plot_symptom_comparison(data_v1) +
        theme(legend.position = "none") +
        labs(title = "Total Group") +
        coord_cartesian(xlim = xlim)
    plot_v2 <- plot_symptom_comparison(data_v2) +
        labs(title = "mGFR < 20") +
        coord_cartesian(xlim = xlim)
    
    grid.arrange(
        plot_v1, plot_v2,
        layout_matrix = matrix(1:2, nrow = 1),
        widths = c(1.7,2)
    )
}
