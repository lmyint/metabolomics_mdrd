contains_item <- function(x, item) {
    item %in% x
}

table_to_string <- function(x) {
    tab <- table(x)
    tab_str <- paste(names(tab), ":", tab)
    paste(tab_str, collapse = ".   ")
}