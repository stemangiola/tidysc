#' Join datasets
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
	tts = dplyr:::flatten_bindable(rlang::dots_values(...))
	
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
