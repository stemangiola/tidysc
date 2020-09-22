#' Join datasets
#' 
#' @importFrom rlang dots_values
#' @importFrom rlang flatten_if
#' @importFrom rlang is_spliced
#' 
#' @export
bind_rows <- function(...) {
	UseMethod("bind_rows")
}

#' @export
bind_rows.default <-  function(...)
{
	dplyr::bind_rows(...)
}

#' @export
bind_rows.tidysc <- function(...)
{
  tts = flatten_if(dots_values(...), is_spliced)
	
	par1 = tts[[1]] %>% attr("parameters") %>% unlist
	par2 = tts[[2]] %>% attr("parameters") %>% unlist
	
	.sample_1 = tts[[1]] %>% attr("parameters") %$%
		# Parameters of the two objects must match
		error_if_parameters_not_match(par1, par2)
	
	par =
		unique(c(par1 %>% names, par2 %>% names)) %>%
		map(~ switch(par1[[.x]] %>% is.null %>% sum(1), par1[[.x]], par2[[.x]])) %>%
		setNames(par1 %>% names)
	
	# Check if cell with same name
	seurat_object = merge(
		tts[[1]] %>% attr("seurat") %>% `[[` (1),
		tts[[2]] %>% attr("seurat") %>% `[[` (1),
		add.cell.ids = 1:2
	)
	
	new.arguments = c(seurat_object = seurat_object %>% list %>% list, par)
	
	create_tt_from_seurat(
		seurat_object = new.arguments$seurat_object,
		min.transcripts = new.arguments$min.transcripts,
		min.cells = new.arguments$min.cells,
		high.mito.thresh = new.arguments$high.mito.thresh,
		high.umi.thresh = new.arguments$high.umi.thresh,
		genome = new.arguments$genome,
		species = new.arguments$species,
		.sample = !!new.arguments$.sample,
		.cell = !!new.arguments$.cell
	) %>%
		
		# Add parameters attribute
		add_attr((.) %>% attr("parameters") %>% c(
			list(
				.transcript = new.arguments$.transcript,
				.abundance = new.arguments$.abundance
			)
		),
		"parameters")
}


#' Mutate datasets
#' @export
mutate <- function(.data, ...) {
	UseMethod("mutate")
}

#' @export
mutate.default <-  function(.data, ...)
{
	dplyr::mutate(.data, ...)
}

#' @export
mutate.tidysc <- function(.data, ...)
{
	mutate_update_and_add_attr(.data, ...)
}

#' left_join datasets
#' @export
left_join <- function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),		 ...)  {
	UseMethod("left_join")
}

#' @export
left_join.default <-  function (x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),    ...)
{
	dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...)
}
#' @export
left_join.tidysc <- function(x, y, by = NULL, copy = FALSE, suffix = c(".x", ".y"),    ...)
{
	
	# Drop class to use dplyr
	class(x) = class(x)[-c(1:2)]
	
	dplyr::left_join(x, y, by = by, copy = copy, suffix = suffix, ...) %>%
		add_attr(x %>% attr("seurat"), "seurat") %>%
		add_attr(x %>% attr("parameters"), "parameters") %>%
		update_metadata_sc() %>%
		
		# Add tt class
		add_class("tt") %>%
		add_class("tidysc")
	
	
}

#' unite columns
#' @export
unite <- function(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE) {
	UseMethod("unite")
}

#' @export
unite.default <-  function(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE)
{
	tidyr::unite(data, col, ..., sep = sep, remove = remove, na.rm = na.rm)
}

#' @export
unite.tidysc <- function(data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE)
{
	unite_update_and_add_attr(data, col, ..., sep = sep, remove = remove, na.rm = na.rm)
}

#' Subset rows using column values
#'
#' `filter()` retains the rows where the conditions you provide a `TRUE`. Note
#' that, unlike base subsetting with `[`, rows where the condition evaluates
#' to `NA` are dropped.
#'
#' dplyr is not yet smart enough to optimise filtering optimisation
#' on grouped datasets that don't need grouped calculations. For this reason,
#' filtering is often considerably faster on [ungroup()]ed data.
#'
#' @section Useful filter functions:
#'
#' * [`==`], [`>`], [`>=`] etc
#' * [`&`], [`|`], [`!`], [xor()]
#' * [is.na()]
#' * [between()], [near()]
#'
#' @section Grouped tibbles:
#'
#' Because filtering expressions are computed within groups, they may
#' yield different results on grouped tibbles. This will be the case
#' as soon as an aggregating, lagging, or ranking function is
#' involved. Compare this ungrouped filtering:
#'
#'
#' The former keeps rows with `mass` greater than the global average
#' whereas the latter keeps rows with `mass` greater than the gender
#'
#' average.
#' @family single table verbs
#' @inheritParams arrange
#' @param ... <[`tidy-eval`][dplyr_tidy_eval]> Logical predicates defined in
#'   terms of the variables in `.data`.
#'   Multiple conditions are combined with `&`. Only rows where the
#'   condition evaluates to `TRUE` are kept.
#' @param .preserve when `FALSE` (the default), the grouping structure
#'   is recalculated based on the resulting data, otherwise it is kept as is.
#' @return
#' An object of the same type as `.data`.
#'
#' * Rows are a subset of the input, but appear in the same order.
#' * Columns are not modified.
#' * The number of groups may be reduced (if `.preserve` is not `TRUE`).
#' * Data frame attributes are preserved.
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' @seealso [filter_all()], [filter_if()] and [filter_at()].
#' @export
#' @examples
#'
#' # Learn more in ?dplyr_tidy_eval
############# START ADDED tidysc #####################################
#' @export
filter <- function (.data, ..., .preserve = FALSE)  {
	UseMethod("filter")
}

#' @export
filter.default <-  function (.data, ..., .preserve = FALSE)
{
	dplyr::filter(.data, ..., .preserve = .preserve)
}

#' @export
filter.tidysc <- function (.data, ..., .preserve = FALSE)
{
	
	.data %>%
		drop_class(c("tidysc", "tt")) %>%
		dplyr::filter( ..., .preserve = .preserve) %>%
		
		# Update seurat
		update_object_sc() %>%
		
		# Add class
		add_class("tt") %>%
		add_class("tidysc")
	
}



#' Subset columns using their names and types
#'
#' @description
#'
#' Select (and optionally rename) variables in a data frame, using a concise
#' mini-language that makes it easy to refer to variables based on their name
#' (e.g. `a:f` selects all columns from `a` on the left to `f` on the
#' right). You can also use predicate functions like [is.numeric] to select
#' variables based on their properties.
#'
#'
#' ## Overview of selection features
#'
#' ```{r, child = "man/rmd/overview.Rmd"}
#' ```
#'
#' @inheritParams arrange
#' @param ... <[`tidy-select`][dplyr_tidy_select]> One or more unquoted
#'   expressions separated by commas. Variable names can be used as if they
#'   were positions in the data frame, so expressions like `x:y` can
#'   be used to select a range of variables.
#' @return
#' An object of the same type as `.data`. The output has the following
#' properties:
#'
#' * Rows are not affected.
#' * Output columns are a subset of input columns, potentially with a different
#'   order. Columns will be renamed if `new_name = old_name` form is used.
#' * Data frame attributes are preserved.
#' * Groups are maintained; you can't select off grouping variables.
#'
#' @section Methods:
#' This function is a **generic**, which means that packages can provide
#' implementations (methods) for other classes. See the documentation of
#' individual methods for extra arguments and differences in behaviour.
#'
#' The following methods are currently available in loaded packages:
#' \Sexpr[stage=render,results=rd]{dplyr:::methods_rd("select")}.
#'
#' @section Examples:
#' ```{r, child = "man/rmd/setup.Rmd"}
#' ```
#'
#' Here we show the usage for the basic selection operators. See the
#' specific help pages to learn about helpers like [starts_with()].
#'
#' The selection language can be used in functions like
#' `dplyr::select()` or `tidyr::pivot_longer()`. Let's first attach
#' the tidyverse:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' library(tidyverse)
#'
#' # For better printing
#' iris <- as_tibble(iris)
#' ```
#'
#' Select variables by name:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(height)
#'
#' iris %>% pivot_longer(Sepal.Length)
#' ```
#'
#' Select multiple variables by separating them with commas. Note how
#' the order of columns is determined by the order of inputs:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(homeworld, height, mass)
#' ```
#'
#' Functions like `tidyr::pivot_longer()` don't take variables with
#' dots. In this case use `c()` to select multiple variables:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' iris %>% pivot_longer(c(Sepal.Length, Petal.Length))
#' ```
#'
#' ## Operators:
#'
#' The `:` operator selects a range of consecutive variables:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(name:mass)
#' ```
#'
#' The `!` operator negates a selection:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' starwars %>% select(!(name:mass))
#'
#' iris %>% select(!c(Sepal.Length, Petal.Length))
#'
#' iris %>% select(!ends_with("Width"))
#' ```
#'
#' `&` and `|` take the intersection or the union of two selections:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' iris %>% select(starts_with("Petal") & ends_with("Width"))
#'
#' iris %>% select(starts_with("Petal") | ends_with("Width"))
#' ```
#'
#' To take the difference between two selections, combine the `&` and
#' `!` operators:
#'
#' ```{r, comment = "#>", collapse = TRUE}
#' iris %>% select(starts_with("Petal") & !ends_with("Width"))
#' ```
#'
#' @family single table verbs
#' @export
select <- function(.data, ...) {
	UseMethod("select")
}

#' @export
select.default <-  function (.data, ...)
{
	dplyr::select(.data, ...)
}

#' @export
select.tidysc <- function (.data, ...)
{
	
	.data %>%
		drop_class(c("tidysc", "tt")) %>%
		dplyr::select( ...) %>%
		
		# Update seurat
		add_attr(.data %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters") %>%
		update_metadata_sc() %>%
		
		# Add class
		add_class("tt") %>%
		add_class("tidysc")
	
}
