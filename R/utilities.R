

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#' @importFrom purrr as_mapper
#'
#' @param input.df A tibble
#' @param condition A boolean
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
  switch(.p %>% `!` %>% sum(1),
         as_mapper(.f1)(.x),
         if (.f2 %>% is.null %>% `!`)
           as_mapper(.f2)(.x)
         else
           .x)

}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @import tidyr
#'
#' @param input.df A tibble
#' @param condition A boolean
#' @return A tibble
ifelse2_pipe = function(.x, .p1, .p2, .f1, .f2, .f3 = NULL) {
  # Nested switch
  switch(# First condition
    .p1 %>% `!` %>% sum(1),

    # First outcome
    as_mapper(.f1)(.x),
    switch(
      # Second condition
      .p2 %>% `!` %>% sum(1),

      # Second outcome
      as_mapper(.f2)(.x),

      # Third outcome - if there is not .f3 just return the original data frame
      if (.f3 %>% is.null %>% `!`)
        as_mapper(.f3)(.x)
      else
        .x
    ))
}

#' Get matrix from tibble
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @return A matrix
as_matrix <- function(tbl,
                      rownames = NULL,
                      do_check = TRUE) {
  rownames = enquo(rownames)
  tbl %>%

    # Through warning if data frame is not numerical beside the rownames column (if present)
    ifelse_pipe(
      do_check &&
        tbl %>%
        # If rownames defined eliminate it from the data frame
        ifelse_pipe(!quo_is_null(rownames), ~ .x[, -1], ~ .x) %>%
        dplyr::summarise_all(class) %>%
        tidyr::gather(variable, class) %>%
        pull(class) %>%
        unique() %>%
        `%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
      ~ {
        warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
        .x
      }
    ) %>%
    as.data.frame() %>%

    # Deal with rownames column if present
    ifelse_pipe(
      !quo_is_null(rownames),
      ~ .x %>%
        magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
        select(-1)
    ) %>%

    # Convert to matrix
    as.matrix()
}

#' Check whether a numeric vector has been log transformed
#'
#' @param x A numeric vector
error_if_log_transformed <- function(x, .abundance) {
  .abundance = enquo(.abundance)

  if (x %>% nrow %>% `>` (0))
    if (x %>% summarise(m = !!.abundance %>% max) %>% pull(m) < 50)
      stop(
        "The input was log transformed, this algorithm requires raw (un-normalised) read counts"
      )
}

#' Check whether there are duplicated genes/transcripts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble of read counts
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the gene/transcript name column
#' @param .abundance A character name of the read count column
error_if_duplicated_genes <- function(input.df,
                                      .sample = `sample`,
                                      .transcript = `transcript`,
                                      .abundance = `read count`) {
  .sample = enquo(.sample)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)

  duplicates <-
    input.df %>%
    select(!!.sample,!!.transcript,!!.abundance) %>%
    distinct() %>%
    count(!!.sample,!!.transcript) %>%
    filter(n > 1) %>%
    arrange(n %>% desc())

  if (duplicates %>% nrow() > 0) {
    writeLines("Those are the duplicated genes")
    duplicates %>% print()
    stop(
      "Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding."
    )
  }

  input.df

}

#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble of read counts
#' @param .abundance A character name of the read count column
error_if_counts_is_na = function(input.df, .abundance) {
  .abundance = enquo(.abundance)

  # Do the check
  if (input.df %>% filter(!!.abundance %>% is.na) %>% nrow %>% `>` (0))
    stop("You have NA values in your counts")

  # If all good return original data frame
  input.df
}

#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param input.df A tibble of read counts
#' @param .abundance A character name of the read count column
error_if_wrong_input = function(input.df, list_input, expected_type) {




  # Do the check
  if (
    list_input %>%
    map(~ .x %>% class() %>% `[` (1)) %>%
    unlist %>%
    equals(expected_type) %>%
    `!`
  )
    stop("You have passed the wrong argument to the function. Please check again.")

  # If all good return original data frame
  input.df
}


#' Formula parser
#'
#' @param fm A formula
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
  if (attr(terms(fm), "response") == 1)
    stop("The formula must be of the kind \"~ covariates\" ")
  else
    as.character(attr(terms(fm), "variables"))[-1]
}

#' Scale design matrix
#'
#' @param df A tibble
#' @return A tibble
#'
#'
scale_design = function(df, formula) {
  df %>%
    setNames(c("sample_idx", "(Intercept)", parse_formula(formula))) %>%
    gather(cov, value,-sample_idx) %>%
    group_by(cov) %>%
    mutate(value = ifelse(
      !grepl("Intercept", cov) &
        length(union(c(0, 1), value)) != 2,
      scale(value),
      value
    )) %>%
    ungroup() %>%
    spread(cov, value) %>%
    arrange(as.integer(sample_idx)) %>%
    select(`(Intercept)`, one_of(parse_formula(formula)))
}

#' Add attribute to abject
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}

#' Remove class to abject
#'
#'
#' @param var A tibble
#' @param name A character name of the class
#'
#' @return A tibble with an additional attribute
drop_class = function(var, name) {
	class(var) <- class(var)[!class(var)%in%name]
	var
}

#' From rlang deprecated
prepend = function (x, values, before = 1)
{
  n <- length(x)
  stopifnot(before > 0 && before <= n)
  if (before == 1) {
    c(values, x)
  }
  else {
    c(x[1:(before - 1)], values, x[before:n])
  }
}

#' Add class to abject
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {

  class(var) <- prepend(class(var),name)

  var
}


#' Sub function of drop_redundant_elements_though_reduced_dimensions
#'
#' @param input.df A tibble
#' @return A tibble with pairs to drop
select_closest_pairs = function(df) {
  couples <- df %>% head(n = 0)

  while (df %>% nrow() > 0) {
    pair <- df %>%
      arrange(dist) %>%
      head(n = 1)
    couples <- couples %>% bind_rows(pair)
    df <- df %>%
      filter(
        !`sample 1` %in% (pair %>% select(1:2) %>% as.character()) &
          !`sample 2` %in% (pair %>% select(1:2) %>% as.character())
      )
  }

  couples

}

#' Set defauly assay of a seurat object
#'
#' @import Seurat
#'
#' @param var A seurat object
#' @param assay_name An character name of the assay
#'
#' @return A seurat object
set_default_assay_seurat = function(var, assay_name) {
  DefaultAssay(object = var) <- assay_name
  var
}

#' Add assay to a seurat object
#'
#' @import Seurat
#'
#' @param var A seurat object
#' @param assay A seurat assay
#' @param assay_name An character name of the assay
#'
#' @return A seurat object
add_assay = function(var, assay, assay_name) {
  var@assays[[assay_name]] = assay
  var
}

#' Drop assay of a seurat object
#'
#' @import Seurat
#'
#' @param var A seurat object
#' @param assay_name An character name of the assay
#'
#' @return A seurat object
drop_assay = function(var, assay_name) {
  var@assays[[assay_name]] = NULL
  var
}

#' Get column names either from user or from attributes
#'
#' @importFrom rlang quo_is_symbol
#'
#' @param input.df A tibble
#' @param elements_column A character name of the sample column
#'
#' @return A list of column enquo or error
get_cell = function(input.df, .cell) {
  # If setted by the user, enquo those
  if (.cell %>% quo_is_symbol())
    return(list(.cell = .cell))

  # Otherwise check if attribute exists
  else {
    # If so, take them from the attribute
    if (input.df %>% attr("parameters") %>% is.null %>% `!`)
      return(list(.cell =   unlist(attr(
        input.df, "parameters"
      ))$.cell))

    # Else through error
    else
      stop(
        "
        You might have altered the tt object and lost the attributes.
        The fucntion does not know what are your cell column (e.g., cell) are.
        You have to either enter those as symbols (e.g., `cell`),
        or use the funtion create_tt_from_tibble() to pass your column names that will be remembered.
        "
      )
  }
}

#' Update seurat attribute based on cell content of the data tt frame
#'
#' @import furrr
#' @importFrom purrr map
#' @importFrom dplyr pull
#' @import Seurat
#'
#' @param var A tt object
#' @param assay A seurat assay
#' @param assay_name An character name of the assay
#'
#' @return A tt object
update_object_sc = function(input.df, .cell = NULL) {
  # Get column names
  .cell = enquo(.cell)
  col_names = get_cell(input.df, .cell)
  .cell = col_names$.cell

  input.df %>%
    add_attr((.) %>%
               attr("seurat") %>%
               map(~ subset(
                 .x,
                 cells = input.df %>%
                   pull(!!.cell) %>%
                   as.character()
               )),
             "seurat")
}

update_metadata_sc = function(input.df, .cell = NULL) {
	# Get column names
	.cell = enquo(.cell)
	col_names = get_cell(input.df, .cell)
	.cell = col_names$.cell

	seurat_obj = input.df %>% attr("seurat")

	# update meta.data with tibble
	seurat_obj[[1]]@meta.data = input.df %>% as.data.frame(row.names = quo_name(.cell))

	input.df %>% add_attr(seurat_obj, "seurat")
}

#' Check if sample already set by the user, otherwise take sample information from Seurat object
#' I should check it with parameters attribute if .sample is present
#' For the moment is fine like this
#'
#' @import dplyr
#' @import Seurat
#'
#' @param input.df A seurat object
#'
#' @return A seurat object
rename_sample_if_samples_not_set_by_user = function(input.df) {
  input.df %>%
    ifelse_pipe(
      (.) %>% distinct(`orig.ident`) %>% as.character == 1,
      ~ .x %>% select(-`orig.ident`),
      ~ .x %>% rename(sample = `orig.ident`)
    )
}

error_if_parameters_not_match = function(par1, par2){

	# Covert enquo to strings
	par1$.sample = quo_name(par1$.sample)
	par1$.cell = quo_name(par1$.cell)
	par1$.transcript = quo_name(par1$.transcript)
	par1$.abundance = quo_name(par1$.abundance)

	par2$.sample = quo_name(par2$.sample)
	par2$.cell = quo_name(par2$.cell)
	par2$.transcript = quo_name(par2$.transcript)
	par2$.abundance = quo_name(par2$.abundance)

	par = foreach(n = par1 %>% names) %do% {
		c(par1[[n]], par2[[n]]) %>% grep("NULL", ., invert = T, value=T) %>% unique

	} %>% setNames(par1 %>% names)

	if(par %>% map_int(~ .x %>% length) %>% unique %>% equals(1) %>% `!`) {
		print(cbind(par1, par2))
		stop("the parameters of the two objects must match")
	}
}
