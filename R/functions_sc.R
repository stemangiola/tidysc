# library(furrr)
# plan(strategy = "multicore", workers = 20)
# options(future.globals.maxSize = 8000 * 1024 ^ 2)

#' Aggregates multiple counts from the same samples/cells (e.g., from isoforms)
#' This function aggregates counts over samples, concatenates other character columns, and averages other numeric columns
#'
#' @importFrom dplyr summarise_all
#' @importFrom dplyr bind_rows
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the gene/transcript name column
#' @param .abundance A character name of the count column
#' @param aggregation_function A function for counts aggregation (e.g., sum)
#' @return A tibble with aggregated genes and annotation
#'
#' @export
aggregate_duplicated_transcripts_sc =
	function(.data,
					 aggregation_function = sum,
					 .sample,
					 .cell ,
					 .transcript,
					 .abundance,
					 keep_integer = T) {

		# Get column names
		.sample = enquo(.sample)
		.cell = enquo(.cell)
		.transcript = enquo(.transcript)
		.abundance = enquo(.abundance)

		# Robust paste function that preserves NAs
		paste3 <- function(..., sep = ", ") {
			L <- list(...)
			L <- lapply(L, function(x) {
				x[is.na(x)] <- ""
				x
			})
			ret <- gsub(paste0("(^", sep, "|", sep, "$)"),
									"",
									gsub(paste0(sep, sep), sep,
											 do.call(paste, c(
											 	L, list(sep = sep)
											 ))))
			is.na(ret) <- ret == ""
			ret
		}

		# Through warning if there are logicals of factor in the data frame
		# because they cannot be merged if they are not unique
		if ((lapply(.data, class) %>% unlist %in% c("logical", "factor")) %>% any) {
			warning("for aggregation fctors and logical columns were converted to character")
			writeLines("Converted to characters")
			lapply(.data, class) %>% unlist %>% `[` (. %in% c("logical", "factor") %>% which) %>% print
		}

		# Select which are the numerical columns
		numerical_columns = .data %>% ungroup() %>% select_if(is.numeric) %>% select(-!!.abundance) %>% colnames() %>% c("n_aggr")

		# ggregates read .data over samples, concatenates other character columns, and averages other numeric columns
		.data %>%

			# Through error if some counts are NA
			error_if_counts_is_na(!!.abundance) %>%

			# transform logials and factors
			mutate_if(is.factor, as.character) %>%
			mutate_if(is.logical, as.character) %>%

			# Add the nuber of duplicates for each gene
			left_join(
				(.) %>% count(!!.sample, !!.cell, !!.transcript, name = "n_aggr"),
				by = c(quo_name(.sample), quo_name(.cell), quo_name(.transcript))
			) %>%

			# Anonymous function - binds the unique and the reduced genes,
			# in the way we have to reduce redundancy just for the duplicated genes
			# input: tibble
			# output tibble distinct
			{
				dplyr::bind_rows(
					# Unique symbols
					(.) %>%
						filter(n_aggr == 1),

					# Duplicated symbols
					(.) %>%
						filter(n_aggr > 1) %>%
						group_by(!!.sample, !!.cell, !!.transcript) %>%
						mutate(
							!!.abundance := !!.abundance %>% aggregation_function()
						) %>%
						mutate_at(vars(numerical_columns), aggregation_function) %>%
						mutate_at(
							vars(-group_cols(),-!!.abundance,-!!numerical_columns),
							list( ~ paste3(unique(.), collapse = ", "))
						) %>%
						distinct()
				)
			} %>%

			# Rename column of number of duplicates for each gene
			rename(`number of merged transcripts` = n_aggr)

	}

#' Aggregates multiple counts from the same samples/cells in a wide format (e.g., from isoforms)
#' This function aggregates counts over samples, concatenates other character columns, and averages other numeric columns
#'
#' @importFrom dplyr summarise_all
#' @importFrom dplyr bind_rows
#' @import Seurat
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the gene/transcript name column
#' @param .abundance A character name of the count column
#' @param aggregation_function A function for counts aggregation (e.g., sum)
#' @return A tibble with aggregated genes and annotation
#'
#' @export
aggregate_duplicated_transcripts_wide_sc =
	function(.data,
					 aggregation_function = sum,
					 .sample,
					 .transcript,
					 keep_integer = T) {

		# Get column names
		.sample = enquo(.sample)
		.transcript = enquo(.transcript)

		# Robust paste function that preserves NAs
		paste3 <- function(..., sep = ", ") {
			L <- list(...)
			L <- lapply(L, function(x) {
				x[is.na(x)] <- ""
				x
			})
			ret <- gsub(paste0("(^", sep, "|", sep, "$)"),
									"",
									gsub(paste0(sep, sep), sep,
											 do.call(paste, c(
											 	L, list(sep = sep)
											 ))))
			is.na(ret) <- ret == ""
			ret
		}

		# Through warning if there are logicals of factor in the data frame
		# because they cannot be merged if they are not unique
		if ((lapply(.data, class) %>% unlist %in% c("logical", "factor")) %>% any) {
			warning("for aggregation fctors and logical columns were converted to character")
			writeLines("Converted to characters")
			lapply(.data, class) %>% unlist %>% `[` (. %in% c("logical", "factor") %>% which) %>% print
		}

		# Select which are the numerical columns
		numerical_columns = .data %>% ungroup() %>% select_if(is.numeric) %>% colnames() %>% c("n_aggr")
		integer_columns = .data %>% ungroup() %>% select_if(is.integer) %>% colnames() %>% c("n_aggr")

		# ggregates read .data over samples, concatenates other character columns, and averages other numeric columns
		.data %>%

			# Through error if some counts are NA
			#error_if_counts_is_na(!!.abundance) %>%

			# transform logials and factors
			mutate_if(is.factor, as.character) %>%
			mutate_if(is.logical, as.character) %>%

			# Add the nuber of duplicates for each gene
			left_join(
				(.) %>% count(!!.sample, !!.transcript, name = "n_aggr"),
				by = c(quo_name(.sample), quo_name(.transcript))
			) %>%

			# Anonymous function - binds the unique and the reduced genes,
			# in the way we have to reduce redundancy just for the duplicated genes
			# input: tibble
			# output tibble distinct
			{
				dplyr::bind_rows(
					# Unique symbols
					(.) %>%
						filter(n_aggr == 1),

					# Duplicated symbols
					(.) %>%
						filter(n_aggr > 1) %>%
						group_by(!!.sample,  !!.transcript) %>%

						# Mutate integer
						mutate_at(vars(integer_columns), aggregation_function, na.rm = T) %>%

						# Mutate reals
						mutate_at(vars(numerical_columns), aggregation_function, na.rm = T) %>%

						# Mutate all the rest
						mutate_at(
							vars(-group_cols(),-!!numerical_columns, -!!integer_columns),
							list( ~ paste3(unique(.), collapse = ", "))
						) %>%

						# Take unique
						distinct()
				)
			} %>%

			# Rename column of number of duplicates for each gene
			rename(`number of merged transcripts` = n_aggr)

	}

# Add variable 500 genes to tibble
add_variable_genes_to_tt = function(.data, n = 500){

	seurat_object = .data %>% attr("seurat") %>% `[[` (1)

	seurat_object@assays$SCT@counts[xx@assays$SCT@var.features,] %>%
		as_tibble(rownames="transcript") %>%
		inner_join(

			# Select top features
			seurat_object@assays$SCT@meta.features %>%
				as_tibble(rownames = "transcript") %>%
				arrange(sct.residual_variance %>% desc) %>%
				head(n=n) %>%
				select(transcript) ,
			by = "transcript"
		)

}

#' Create tt object from seurat object
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom purrr map
#' @import Seurat
#'
#' @param seurat_object A seurat object
#' @param dir_names A character array with the output directories of cellRanger if run from internal function
#' @param min.transcripts An integer with the threshold of the minimum number of genes allowed for each cell
#' @param min.cells An integer with the threshold of the minimum anumber of cells allowed for each sample
#' @param high.mito.thresh An numeric with the threshold of the maximum fraction of reads from mitochondrion. Used for filtering
#' @param high.umi.thresh An numeric with the threshold of the maximum fraction of umis. Used for filtering
#' @param genome A character name of the mapping genome used
#' @param .sample A symbol for the sample column
#' @param .cell A symbol for the cell column
#' @param species A character name of the species
#'
#' @return A tt object
create_tt_from_seurat = function(seurat_object,
																 dir_names = list(""),
																 min.transcripts = 400,
																 min.cells = 5,
																 high.mito.thresh = 0.08,
																 high.umi.thresh = 10000,
																 .sample = `sample`,
																 .cell = `cell`,
																 species,
																 genome = ifelse(species == "Human", "hg38", "mm10")) {
	writeLines("Converting Seurat object back to tibble")

	# Parse column names
	.sample = enquo(.sample)
	.cell = enquo(.cell)

	# Create object
	seurat_object %>%
		map_dfr(
			~ .x@meta.data %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				mutate_if(is.factor, as.character)
		)  %>%
		mutate_if(is.character, as.factor) %>%
		select(-one_of(c("count total", "gene count"))) %>%
		rename(`count total` = nCount_RNA,
					 `gene count` = nFeature_RNA) %>%

		# Pass the sample column instead of origin.ident
		select(-one_of(quo_name(.sample))) %>%
		mutate(!!.sample := `orig.ident`) %>%
		select(-`orig.ident`) %>%

		# # Add transcript counts
		# get_transcript_nested_tibble = function(seurat_object){
		# 	map2(
		# 		seurat_object@assays, seurat_object@assays %>% names,
		# 		~ {
		# 			assay = .x
		# 			assay@counts %>%
		# 				colnames %>%
		# 				as_tibble() %>%
		# 				setNames(quo_name(.cell)) %>%
		# 				mutate(!!sprintf("counts %s", .y) := map(cell, ~ assay[,.x] %>% Matrix::Matrix(sparse = TRUE)  ))
		# 		}
		# 	)
		# }

		# Reorde columns logically
		select(!!.sample,!!.cell, everything()) %>%

		# Add Seurat object
		add_attr(seurat_object, "seurat") %>%

		# Add parameters object
		add_attr(list(
								 		min.transcripts = min.transcripts,
								 		min.cells = min.cells,
								 		high.mito.thresh = high.mito.thresh,
								 		high.umi.thresh = high.umi.thresh,
								 		genome = genome,
								 		species = species,
								 		.sample = enquo(.sample),
								 		.cell = enquo(.cell)
								 	),
						 "parameters") %>%

		# Add mitochondrion information if not present already
		ifelse_pipe(
			!("mito.fraction" %in% ((.) %>% attr("seurat") %>% `[[` (1) %>% `@` (meta.data) %>% colnames)),
			~ .x %>% add_mitochndrion_transcription_abundance_sc(.cell = !!.cell)
		) %>%

		# Add cell cycle information
		ifelse_pipe(
			!("Phase" %in% ((.) %>% attr("seurat") %>% `[[` (1) %>% `@` (meta.data) %>% colnames)),
			~ .x %>% add_cell_cycle_annotation_sc(.cell = !!.cell)
		) %>%

		# Add tt class
		add_class("tt") %>%
		add_class("tidysc") %>%

		# Filter dead cells
		mutate(low_quality =  !(`count total` > 200 & `mito.fraction` < 0.1))

}

#' Create seurat object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param .data A tibble
#' @param .sample A symbol for the sample column
#' @param .cell A symbol for the cell column
#' @param .transcript A symbol for the transcript name column
#' @param .abundance A symbol for counts
#'
#' @return A tt object
create_seurat_from_tibble = function(.data,
																		 .sample,
																		 .cell,
																		 .transcript,
																		 .abundance,
																		 min.transcripts = 400,
																		 min.cells = 5,
																		 ...) {
	writeLines("Creating seurat object")

	# Prepare column name enquo
	.sample = enquo(.sample)
	.cell = enquo(.cell)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	# Check that there are not multiple information for each cell type, that each cell type is not repeted row-wise
	if (.data %>%
			select(-!!.transcript,-!!.abundance) %>%
			distinct() %>%
			count(!!.cell) %>%
			pull(n) %>%
			max %>%
			`>` (1))
		stop(
			"The annotations in the input dataframe must be in one record for each cell.
      That is,
      YOUR_DATA %>%
      select(.transcript, .abundance) %>%
      distinct() %>%
      count(!!.cell) %>%
      pull(n) %>%
      max %>% `==` (1)"
		)

	# Check if data rectangular
	if(.data %>% count(!!.cell) %>% count(n) %>% nrow %>% `>` (1))
		stop("The input data frame does not represent a rectangular structure. Each transcript must be present in all cells")

	# Set sample names
	sample_names =
		.data %>%
		select(-!!.transcript,-!!.abundance) %>%
		distinct() %>%
		pull(!!.sample)

	# Create Seurat object
	seurat_object =
		CreateSeuratObject(
			# counts
			.data %>%
				distinct(!!.cell,!!.transcript,!!.abundance) %>%
				spread(!!.cell,!!.abundance) %>%
				data.frame(row.names = quo_name(.transcript)),

			# Sample/cell information
			meta.data =
				.data %>%
				mutate_if(is.factor, as.character) %>%
				select(-!!.transcript,-!!.abundance) %>%
				distinct() %>%
				data.frame(row.names = quo_name(.cell)),

			# Other parameters
			min.cells = min.cells,
			min.features = min.transcripts

		)

	# pass the sample column to the orig.ident column of meta data
	seurat_object@meta.data$orig.ident =
		colnames(seurat_object) %>%
		as_tibble() %>%
		rename(!!.cell := value) %>%
		left_join(
			.data %>%
				distinct(!!.sample,!!.cell),
			by = quo_name(.cell)
		)  %>%
		pull(!!.sample) # seurat_object@meta.data[[quo_name(.sample)]]

	# Return
	seurat_object

}

#' Create seurat object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param .data A tibble
#' @param .sample A symbol for the sample column
#' @param .cell A symbol for the cell column
#' @param .transcript A symbol for the transcript name column
#' @param .abundance A symbol for counts
#'
#' @return A tt object
create_seurat_from_tibble_wide = function(.data,
																					.sample,
																					transcript_col_names,
																					.cell,
																					min.transcripts = 400,
																					min.cells = 5,
																					...) {
	writeLines("Creating seurat object")

	# Prepare column name enquo
	.sample = enquo(.sample)
	.cell = enquo(.cell)

	# Set sample names
	sample_names =
		.data %>%
		pull(!!.sample)

	# Create Seurat object
	seurat_object =
		CreateSeuratObject(
			# counts
			.data %>%
				select(!!.cell, transcript_col_names) %>%
				data.frame(row.names = quo_name(.cell)) %>%
				t(),

			# Sample/cell information
			meta.data =
				.data %>%
				select(-transcript_col_names) %>%
				data.frame(row.names = quo_name(.cell)),

			# Other parameters
			min.cells = min.cells,
			min.features = min.transcripts

		)

	# pass the sample column to the orig.ident column of meta data
	seurat_object@meta.data$orig.ident =
		colnames(seurat_object) %>%
		enframe() %>%
		rename(!!.cell := value) %>%
		left_join(
			.data %>%
				distinct(!!.sample,!!.cell),
			by = quo_name(.cell)
		)  %>%
		pull(!!.sample) # seurat_object@meta.data[[quo_name(.sample)]]

	# Return
	seurat_object

}


#' Create tt object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map
#' @import Seurat
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .cell A character name of the cell name column
#' @param .transcript A character name of the transcript name column
#' @param .abundance A character name of the count column
#' @param species A character name of the species used (e.g., Mouse, Human)
#'
#' @return A tibble with an additional column
#'
#' @export
create_tt_from_tibble_sc = function(.data,
																		.sample,
																		.cell,
																		.transcript,
																		.abundance,
																		species,
																		min.transcripts = 400,
																		min.cells = 5,
																		high.mito.thresh = 0.08,
																		high.umi.thresh = 10000,
																		genome = ifelse(species == "Human", "hg38", "mm10"),
																		...) {
	writeLines("Start parsing the data frame")

	# Check is species is of the right type
	if (!(species %in% c("Human", "Mouse")))
		stop("Species must be \"Human\" or \"Mouse\"")

	# Prepare column name enquo
	.sample = enquo(.sample)
	.cell = enquo(.cell)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	convert_underscore_to_dash = function(.data, .transcript){

		.transcript = enquo(.transcript)

		wired_names =
			.data %>%
			distinct(!!.transcript) %>%
			filter(grepl("_", !!.transcript)) %>%
			pull(!!.transcript)

		.data %>%

			# Check if wired names are present
			ifelse_pipe(
				wired_names %>% length %>% `>` (0),
				~ {
					warning(sprintf("Genes with possibly a strange name are present %s. Converting \"_\" to \"-\" ", paste(wired_names, collapse=", ")))

					bind_rows(

						# Genes with normal names
						.x %>%
							filter(!(!!.transcript %in% wired_names)),

						# Genes with wired names converted
						.x %>%
							filter((!!.transcript %in% wired_names)) %>%
							mutate(!!.transcript := gsub("_", "-", !!.transcript))
					)
				}
			)
	}

	convert_dash_to_underscore = function(.data, .cell){

		.cell = enquo(.cell)

		wired_names =
			.data %>%
			distinct(!!.cell) %>%
			filter(grepl("-", !!.cell)) %>%
			pull(!!.cell)

		.data %>%

			# Check if wired names are present
			ifelse_pipe(
				wired_names %>% length %>% `>` (0),
				~ {
					warning(sprintf("Genes with possibly a strange name are present %s. Converting \"-\" to \"_\" ", paste(wired_names, collapse=", ")))

					bind_rows(

						# Genes with normal names
						.x %>%
							filter(!(!!.cell %in% wired_names)),

						# Genes with wired names converted
						.x %>%
							filter((!!.cell %in% wired_names)) %>%
							mutate(!!.cell := gsub("-", "_", !!.cell))
					)
				}
			)
	}


	# Create seurat object
	.data %>%

		# Eliminate gene statistics from genes and convert genes with wired name
		filter(!(!!.transcript %in% c("no_feature", "too_low_aQual", "not_aligned", "alignment_not_unique"))) %>%
		convert_underscore_to_dash(!!.transcript) %>%

		# Check if any cell name starts with a number
		# If so put a _ before, save original name
		# and through warning
		ifelse_pipe((.) %>% filter(grepl("^[0-9]",!!.cell)) %>% nrow %>% `>` (1),
								~ {
									warning(
										"
                    some cell names started with a number.
                    This is incompatible with Seurat object, the character _ had been added as prefix.
                    The orignal cell name as been kept as column in the data frame"
									)
									.x %>%
										mutate(!!as.symbol(sprintf("%s original", quo_name(.cell))) := !!.cell) %>%
										mutate(!!.cell := !!.cell %>% paste0("X_", .))
								}) %>%

		# Check if any cell name starts with a number
		# If so put a _ before, save original name
		# and through warning
		convert_dash_to_underscore(!!.cell) %>%

		# Convert characters to factors
		mutate_if(is.character, as.factor) %>%

		# If there is more than one sample add number to cell names
		ifelse_pipe(
			(.) %>% pull(!!.sample) %>% levels %>% length %>% `>` (1),
			~ .x %>%
				mutate(sample_idx = !!.sample %>% as.integer) %>%
				tidyr::unite(!!.cell, c(!!.cell, sample_idx), sep = "_")
		) %>%

		# Anonymous function - create Seurat object
		# input: tibble
		# output Seurat object
		create_seurat_from_tibble(!!.sample,
															!!.cell,
															!!.transcript,
															!!.abundance,
															min.transcripts = min.transcripts,
															min.cells = min.cells,
															...) %>%

		# Create tt object from seurat
		list %>%
		create_tt_from_seurat(
			min.transcripts = min.transcripts,
			min.cells = min.cells,
			high.mito.thresh = high.mito.thresh,
			high.umi.thresh = high.umi.thresh,
			genome = genome,
			.sample = !!.sample,
			.cell = !!.cell,
			species = species
		) %>%

		# Eliminate orig.ident column, because same as sample
		select(-contains("orig.ident")) %>%

		# Add parameters attribute
		add_attr((.) %>% attr("parameters") %>% c(
								 	list(
								 		.transcript = enquo(.transcript),
								 		.abundance = enquo(.abundance)
								 	)
								 ),
						 "parameters")

}

#' Create tt object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map
#' @import Seurat
#'
#' @param .data A tibble
#' @param .sample A character name of the sample column
#' @param .cell A character name of the cell name column
#' @param .transcript A character name of the transcript name column
#' @param .abundance A character name of the count column
#' @param species A character name of the species used (e.g., Mouse, Human)
#'
#' @return A tibble with an additional column
#'
#' @export
create_tt_from_tibble_wide_sc = function(.data,
																				 .sample,
																				 .transcript_position,
																				 .cell,
																				 species,
																				 min.transcripts = 400,
																				 min.cells = 5,
																				 high.mito.thresh = 0.08,
																				 high.umi.thresh = 10000,
																				 genome = ifelse(species == "Human", "hg38", "mm10"),
																				 ...) {
	writeLines("Start parsing the data frame")

	# Check is species is of the right type
	if (!(species %in% c("Human", "Mouse")))
		stop("Species must be \"Human\" or \"Mouse\"")

	# Prepare column name enquo
	.sample = enquo(.sample)
	.cell = enquo(.cell)

	convert_underscore_to_dash = function(.data, transcript_col_names){

		wired_names =
			grep("_", .data %>%
					 	select(transcript_col_names) %>%
					 	colnames, value = T)

		.data %>%

			# Check if wired names are present
			ifelse_pipe(
				wired_names %>% length %>% `>` (0),
				~ {
					warning(sprintf("Genes with possibly a strange name are present %s. Converting \"_\" to \"-\" ", paste(wired_names, collapse=", ")))

					new_col_names_idx = (.data %>% colnames %in% wired_names) %>% which
					new_col_names = .data %>% colnames
					new_col_names[new_col_names_idx] = new_col_names[new_col_names_idx] %>% gsub("_", "-", .)

					.data %>% setNames(new_col_names)

				}
			)
	}

	convert_dash_to_underscore = function(.data, .cell){

		.cell = enquo(.cell)

		wired_names =
			.data %>%
			distinct(!!.cell) %>%
			filter(grepl("-", !!.cell)) %>%
			pull(!!.cell)

		.data %>%

			# Check if wired names are present
			ifelse_pipe(
				wired_names %>% length %>% `>` (0),
				~ {
					warning(sprintf("Genes with possibly a strange name are present %s. Converting \"-\" to \"_\" ", paste(wired_names, collapse=", ")))

					bind_rows(

						# Genes with normal names
						.x %>%
							filter(!(!!.cell %in% wired_names)),

						# Genes with wired names converted
						.x %>%
							filter((!!.cell %in% wired_names)) %>%
							mutate(!!.cell := gsub("-", "_", !!.cell))
					)
				}
			)
	}

	# Update trancript column position
	wired_names =  c("no_feature", "too_low_aQual", "not_aligned", "alignment_not_unique")
	transcript_col_names = .data %>% colnames %>% `[` (.transcript_position)
	transcript_col_names = transcript_col_names[!transcript_col_names %in% wired_names]
	transcript_colnames_converted = transcript_col_names %>% gsub("_", "-", .)

	# Create seurat object
	.data %>%

		# Eliminate gene statistics from genes and convert genes with wired name
		select(-one_of(wired_names)) %>%

		convert_underscore_to_dash(transcript_col_names) %>%

		# Check if any cell name starts with a number
		# If so put a _ before, save original name
		# and through warning
		ifelse_pipe((.) %>% filter(grepl("^[0-9]",!!.cell)) %>% nrow %>% `>` (1),
								~ {
									warning(
										"
                    some cell names started with a number.
                    This is incompatible with Seurat object, the character _ had been added as prefix.
                    The orignal cell name as been kept as column in the data frame"
									)
									.x %>%
										mutate(!!as.symbol(sprintf("%s original", quo_name(.cell))) := !!.cell) %>%
										mutate(!!.cell := !!.cell %>% paste0("X_", .))
								}) %>%

		# Check if any cell name starts with a number
		# If so put a _ before, save original name
		# and through warning
		convert_dash_to_underscore(!!.cell) %>%

		# Convert characters to factors
		mutate_if(is.character, as.factor) %>%

		# If there is more than one sample add number to cell names
		ifelse_pipe(
			(.) %>% pull(!!.sample) %>% levels %>% length %>% `>` (1),
			~ .x %>%
				mutate(sample_idx = !!.sample %>% as.integer) %>%
				unite(!!.cell, c(!!.cell, sample_idx), sep = "_")
		) %>%

		# Anonymous function - create Seurat object
		# input: tibble
		# output Seurat object
		create_seurat_from_tibble_wide(!!.sample,
																	 transcript_colnames_converted,
																	 !!.cell,
																	 min.transcripts = min.transcripts,
																	 min.cells = min.cells,
																	 ...) %>%

		# Create tt object from seurat
		list %>%
		create_tt_from_seurat(
			min.transcripts = min.transcripts,
			min.cells = min.cells,
			high.mito.thresh = high.mito.thresh,
			high.umi.thresh = high.umi.thresh,
			genome = genome,
			.sample = !!.sample,
			.cell = !!.cell,
			species = species
		) %>%

		# Eliminate orig.ident column, because same as sample
		select(-contains("orig.ident")) %>%
	# %>%

		# Add parameters attribute

		# Attach attributes

		add_attr(
			(.) %>%
				attr("parameters") %>%
				c(	.transcript =	(function(x, v) enquo(v))(x, !!as.symbol("transcript"))	),
			"parameters"
		) %>%
		add_attr(
			(.) %>%
				attr("parameters") %>%
				c(	.abundance =	(function(x, v) enquo(v))(x, !!as.symbol("abundance"))	),
			"parameters"
		)



}

#' Create tt object from cellRanger results
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr pmap
#' @import Seurat
#'
#' @param dir_names A character array with the output directories of cellRanger if run from internal function
#' @param min.transcripts An integer with the threshold of the minimum number of genes allowed for each cell
#' @param min.cells An integer with the threshold of the minimum anumber of cells allowed for each sample
#' @param high.mito.thresh An numeric with the threshold of the maximum fraction of reads from mitochondrion. Used for filtering
#' @param high.umi.thresh An numeric with the threshold of the maximum fraction of umis. Used for filtering
#' @param genome A character name of the mapping genome used
#' @param species A character name of the species
#'
#' @return A tt object
#'
#' @export
create_tt_from_cellRanger_sc <- function(dir_names,
																				 min.transcripts = 400,
																				 min.cells = 5,
																				 high.mito.thresh = 0.08,
																				 high.umi.thresh = 10000,
																				 species,
																				 genome = ifelse(species == "Human", "hg38", "mm10")) {
	# n_cores <- system("nproc", intern = TRUE) %>%
	# 	as.integer() %>%
	# 	`-`(2)

	seurat_object =
		pmap(list(dir_names,
							1:length(dir_names), # cell suffix
							rep(length(dir_names), length(dir_names))), # if more than one sample),
				 ~ {
				 	# Create object
				 	my_obj =
				 		..1 %>%
				 		Read10X() %>%
				 		CreateSeuratObject(
				 			min.cells = min.cells,
				 			min.features = min.transcripts,
				 			project = ..1
				 		)

				 	# If there is more than one sample add number to cell names
				 	if (..3 %>% `>` (1))
				 		my_obj = my_obj %>% RenameCells((.) %>% colnames %>% paste(..2, sep = "_"))

				 	my_obj

				 })

		# Merge objects
	if(length(seurat_object) > 1)
	seurat_object =
			merge(x = seurat_object[[1]], y=seurat_object[2:length(seurat_object)]) %>%
			list()

	# If object empty through error
	if (seurat_object %>% length == 0)
		stop(
			"The directory provided seem to be empty. Are you sure you are in the right working directory?"
		)

	# create_tt_from_seurat
	seurat_object %>%
		create_tt_from_seurat(
			min.transcripts = min.transcripts,
			min.cells = min.cells,
			high.mito.thresh = high.mito.thresh,
			high.umi.thresh = high.umi.thresh,
			genome = genome,
			species = species) %>%

		# Eliminate orig.ident column, because same as sample
		select(-contains("orig.ident")) %>%

		# Attach attributes

		add_attr(
			(.) %>%
				attr("parameters") %>%
				c(	.transcript =	(function(x, v) enquo(v))(x, !!as.symbol("transcript"))	),
			"parameters"
		) %>%
		add_attr(
			(.) %>%
				attr("parameters") %>%
				c(	.abundance =	(function(x, v) enquo(v))(x, !!as.symbol("abundance"))	),
			"parameters"
		)

		# # Add parameters attribute
		# add_attr((.) %>% attr("parameters") %>% c(
		# 						 	list(
		# 						 		.transcript = NULL,
		# 						 		.abundance = NULL
		# 						 	)
		# 						 ),
		# 				 "parameters")

		# %>%
		#
		# # Add parameters attribute
		# add_attr(map((.) %>% attr("parameters"),
		# 						 ~ .x %>% c(
		# 						 	list(
		# 						 		.transcript = enquo(`.transcript`),
		# 						 		.abundance = enquo(`.abundance`)
		# 						 	)
		# 						 )),
		# 				 "parameters")

	# Rename orig.ident
	#rename(sample = `orig.ident`)

	# # Get counts
	# mutate(`count tibble` =
	# 			 	map2(`seurat object`, sample,
	# 			 			 ~ {
	# 			 			 	GetAssayData(object = .x, slot = "counts") %>%
	# 			 			 		as.data.frame %>%
	# 			 			 		as_tibble(rownames = "transcript") %>%
	# 			 			 		mutate(transcript = transcript %>% toupper()) %>%
	# 			 			 		gather(sample, `count`,-transcript) %>%
	# 			 			 		mutate(`count` = `count` %>% as.integer) %>%
	# 			 			 		mutate_if(is.character, as.factor)
	# 			 			 })) %>%

}

#' Calculate mitochondrion transcription abundance
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
get_mitochndrion_transcription_abundance_sc = function(.data, .cell) {
	writeLines("Calculating mitochondrion trancription")

	# Set column names
	.cell = enquo(.cell)

	# Cell column name
	# .cell_name = .data %>% attr("parameters")  %$% .cell %>% quo_name

	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	seurat_object =
		.data %>%
		attr("seurat") %>%
		map(~ {
			# calculate the mitichondrion sequences
			mito.genes = grep(pattern = "^MT-|^mt-",
												x = (rownames(.x)),
												value = TRUE)

			# Select mitochondrion transcription
			mito =
				.x %>%
				GetAssayData(slot = "counts") %>%
				`[` (mito.genes, ) %>%

				# Add percent mitochondrion
				Matrix::colSums() %>%
				enframe(name = quo_name(.cell), value = "count mitochondrion total") %>%
				mutate(`count mitochondrion total` = `count mitochondrion total` %>% as.integer) %>%
				#mutate(sample = .x@project.name) %>%
				left_join(
					.data %>%
						select(!!.cell, `count total`) %>%
						mutate_if(is.factor, as.character),
					by = quo_name(.cell)
				) %>%
				mutate(`fraction mitochondrion` = `count mitochondrion total` / `count total`) %>%
				select(!!.cell,
							 `count mitochondrion total`,
							 `fraction mitochondrion`)

			.x %>%
				AddMetaData(metadata = mito %>% pull(`fraction mitochondrion`),
										col.name = "mito.fraction") %>%
				AddMetaData(
					metadata = mito %>% pull(`count mitochondrion total`),
					col.name = "mito.tot"
				)
		})

	seurat_object %>%
		map_dfr(~ {
			# Select mitochondrion transcription
			.x %>%
				`@` (meta.data) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				#mutate(sample = .x@project.name) %>%
				# rename(
				#   `count mitochondrion total` = mito.tot,
				#   `fraction mitochondrion` = mito.fraction
				# ) %>%
				select(!!.cell, mito.fraction, mito.tot)
		}) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add mitochondrion transcription abundance
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
add_mitochndrion_transcription_abundance_sc = function(.data, .cell) {
	# Get column names
	.cell = enquo(.cell)

	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	# Cell column name
	# .cell_name = .data %>% attr("parameters")  %$% .cell %>% quo_name

	# Get now object
	.data.annotated =
		.data %>%
		get_mitochndrion_transcription_abundance_sc(.cell = !!.cell)

	# Merge
	.data %>%
		left_join(.data.annotated  ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Get normalised counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
get_normalised_counts_sc = function(.data, verbose = TRUE) {
	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell
	.sample = .data %>% attr("parameters")  %$% .sample
	# .transcript = .data %>% attr("parameters")  %$% .transcript
	# .abundance = .data %>% attr("parameters")  %$% .abundance

	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	new_assay_name = "normalised"

	seurat_object =
		.data %>%
		attr("seurat") %>%
		map(~ .x %>%	SCTransform(
			verbose = verbose,
			new.assay.name = new_assay_name
		))

	# Get normalised counts
	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (meta.data) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				select(
					!!.cell,
					sprintf("nCount_%s", new_assay_name),
					sprintf("nFeature_%s", new_assay_name)
				) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add normalised counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
add_normalised_counts_sc = function(.data, verbose = TRUE) {
	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	.data.normalised =
		.data %>%
		get_normalised_counts_sc(verbose = verbose)

	.data %>%
		dplyr::left_join(.data.normalised,  by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.normalised %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add doublet classification
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr pmap
#' @import Seurat
#'
#' @param .data A tt object
#' @param doublet.reps A numeric
#'
#' @return A tt object
#'
#' @export
add_doublet_classification_sc = function(.data, doublet.reps) {
	.data %>%
		mutate(`seurat object` =
					 	pmap(list(`seurat object`, dir_names, parameters), ~ {
					 		seurat_object = ..1
					 		file_doublets = ..3 %$% output.dir %>% file.path("/doublets.txt")
					 		file_filtered = ..2 %>% sprintf("%s/outs/filtered_feature_bc_matrix", .) %>%	file.path("matrix.mtx")
					 		barcodes =
					 			.y %>%
					 			file.path("outs/filtered_feature_bc_matrix", "barcodes.tsv.gz") %>%
					 			gzfile() %>%
					 			read.delim(header = F,
					 								 as.is = T,
					 								 sep = "\t") %>%
					 			pull(1) %>%
					 			gsub("-1", "", .)

					 		seurat_object$doublet =
					 			switch(
					 				file_doublets %>% file.exists %>% `!` %>% sum(1),

					 				# If file already exists
					 				file_doublets %>%
					 					read.delim(
					 						sep = "\t",
					 						header = F,
					 						as.is = T
					 					) %>%
					 					pull(1),

					 				# Otherwise calculate
					 				file_filtered %>%
					 					{
					 						if (!file.exists(file_filtered))
					 							system(
					 								command = paste(
					 									"gunzip -c",
					 									paste0(file_filtered, ".gz"),
					 									">",
					 									file_filtered
					 								),
					 								intern = FALSE
					 							)
					 						file_filtered %>% detect_doublets(file_doublets, doublet.reps) %>% setNames(barcodes) %>%
					 							{
					 								# Delete before returning
					 								if (file.exists(paste0(file_filtered, ".gz")))
					 									system(command = paste("rm", file_filtered),
					 												 intern = T)
					 								(.)
					 							}

					 					}
					 			) %>%

					 			# Match to Seurat object
					 			setNames(barcodes)

					 		seurat_object$doublet = seurat_object$doublet[match(colnames(seurat_object), barcodes)]

					 		# seurat_object@meta.data =
					 		# 	seurat_object@meta.data
					 		# %>%
					 		# 	mutate(barcodes = colnames(..1)) %>%
					 		# 	left_join(
					 		# 			seurat_object$doublet %>%
					 		# 			enframe %>%
					 		# 			rename(barcodes = name, doublet = value),
					 		# 		by = "barcodes"
					 		# 	) %>%
					 		#
					 		# 	select(-barcodes)
					 		# return

					 		seurat_object
					 	}))
}

#' Add variable gene annotation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
add_variable_genes_classification = function(.data) {
	# Update on tibble
	.data = .data %>% update_object_sc()

	.data %>%
		add_attr((.) %>%
						 	attr("seurat") %>%
						 	map(~
						 				.x %>%
						 				FindVariableFeatures(selection.method = "vst")),
						 "seurat")
}

#' Add variable gene annotation
#' This function is needed for a bug in Seurat about binning
#' https://github.com/satijalab/seurat/issues/1227
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A seurat object
#'
#' @return A seurat object
iterate_cell_cycle_scoreing = function(x, s.features, g2m.features, StartNBin = 24) {
	# This function is needed for a bug in Seurat about binning
	# https://github.com/satijalab/seurat/issues/1227

	object.cc = NULL

	while (class(object.cc) == "try-error" || is.null(object.cc[1])) {
		object.cc <- try(x %>%
										 	CellCycleScoring(
										 		s.features = s.features,
										 		g2m.features = g2m.features,
										 		set.ident = TRUE,
										 		nbin = StartNBin
										 	),
										 silent = T)

		StartNBin = round(StartNBin / 2)

	}


	# Return
	object.cc

}

#' Get cell cyle annotation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom magrittr equals
#' @importFrom magrittr %$%
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr map2_dfr
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
get_cell_cycle_annotation_sc = function(.data, .cell) {
	# Get column names
	.cell = enquo(.cell)

	# Progress
	writeLines("Classifying cells among cell cycle states")

	# Cell column name
	# .cell = .data %>% attr("parameters")  %$% .cell
	# .sample = .data %>% attr("parameters")  %$% .sample

	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	# Update Seurat object
	seurat_object =

		.data %>%
		attr("seurat") %>%
		map(
				 ~ {
				 	# Get cell cycle genes
				 	cc.genes =
				 		switch(
				 			grepl("mm", .data %>% attr("parameters") %$% genome) %>% `!` %>% sum(1),
				 			cc.genes_mouse,
				 			Seurat::cc.genes
				 		) %>%

				 		# Add both upper and lower case because Seurat likes lower cases for it's genes
				 		map( ~ .x %>% c(toupper(.x)))

				 	# Check that I have any matching genes (i.e., If I gave the wrong species)
				 	if (cc.genes %>% unlist %in% (.x %>% rownames) %>% which %>% length %>% equals(0))
				 		stop(
				 			sprintf(
				 				"
                 You don't have any transcript that is within the cell cycle gene signatures for %s.
                 Are you sure you selected the right species?
                 ",
				 				.data %>% attr("parameters") %$% species
				 			)
				 		)


				 	# Return
				 	iterate_cell_cycle_scoreing(.x,
				 															s.features = cc.genes$s.genes,
				 															g2m.features = cc.genes$g2m.genes)


				 })

	# Convert Seurat object to tibble
	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (meta.data) %>%
				as_tibble(rownames = quo_name(.cell)) %>%

				# Function that I have to sort out better
				select(-`orig.ident`) %>%

				# Prepare the data frame
				select(!!.cell, S.Score, G2M.Score, Phase) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add cell cyle annotation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt object
#'
#' @export
add_cell_cycle_annotation_sc = function(.data, .cell) {
	# Get column names
	.cell = enquo(.cell)


	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	# Cell column name
	# .cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_cell_cycle_annotation_sc(.cell = !!.cell)

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Get principal .dims
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @import Seurat
#'
#' @param .data A tt object
#' @param .dims An integer, the number of .dims to calculate
#'
#' @return A tt object
#'
#' @export
get_reduced_dimensions_PCA = function(.data,
																			.dims = 10) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	seurat_object =

		.data %>%
		add_variable_genes_classification() %>%
		attr("seurat") %>%

		# Scale data for PCA
		map(~
					.x %>%
					#ScaleData(display.progress = T, num.cores=4, do.par = TRUE) %>%
					RunPCA(npcs = .dims %>% max))


	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (reductions) %$% pca %>% `@` (cell.embeddings) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				setNames((.) %>% colnames %>% gsub("_", " ", .)) %>%
				#mutate(sample = .x@project.name) %>%
				select(!!.cell, everything()) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add principal .dims
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A tt object
#' @param .dims An integer, the number of .dims to calculate
#'
#' @return A tt object
#'
#' @export
add_reduced_dimensions_PCA = function(.data, .dims = 10) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_reduced_dimensions_PCA(.dims)

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Get UMAP dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @import Seurat
#'
#' @param .data A tt object
#' @param .dims An integer, the number of .dims to calculate
#'
#' @return A tt object
#'
#' @export
get_reduced_dimensions_UMAP = function(.data, .dims = 10) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	seurat_object =

		.data %>%
		add_variable_genes_classification() %>%
		add_reduced_dimensions_PCA(.dims) %>%

		attr("seurat") %>%

		# Scale data for UMAP
		map(~
					.x %>%
					#ScaleData(display.progress = T, num.cores=4, do.par = TRUE) %>%
					RunUMAP(reduction = "pca", dims = 1:.dims))

	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (reductions) %$% umap %>% `@` (cell.embeddings) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				setNames((.) %>% colnames %>% gsub("_", " ", .)) %>%
				#mutate(sample = .x@project.name) %>%
				select(!!.cell, everything()) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add UMAP dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A tt object
#' @param .dims An integer, the number of .dims to calculate
#'
#' @return A tt object
#'
#' @export
add_reduced_dimensions_UMAP = function(.data, .dims = 10) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_reduced_dimensions_UMAP(.dims = .dims)

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Get TSNE dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @import Seurat
#'
#' @param .data A tt object
#' @param .dims An integer, the number of .dims to calculate
#'
#' @return A tt object
#'
#' @export
get_reduced_dimensions_TSNE = function(.data, .dims = 10) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	seurat_object =

		.data %>%
		add_variable_genes_classification() %>%
		add_reduced_dimensions_PCA(.dims) %>%

		attr("seurat") %>%

		# Scale data for UMAP
		map(~
					.x %>%
					#ScaleData(display.progress = T, num.cores=4, do.par = TRUE) %>%
					RunTSNE(reduction = "pca", dims = 1:.dims))

	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (reductions) %$% tsne %>% `@` (cell.embeddings) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				setNames((.) %>% colnames %>% gsub("_", " ", .)) %>%
				#mutate(sample = .x@project.name) %>%
				select(!!.cell, everything()) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add TSNE dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import Seurat
#'
#' @param .data A tt object
#' @param .dims An integer, the number of .dims to calculate
#'
#' @return A tt object
#'
#' @export
add_reduced_dimensions_TSNE = function(.data, .dims = 10) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_reduced_dimensions_TSNE(.dims = .dims)

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Integrate many seurat objects
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map_int
#' @import Seurat
#'
#' @param .data A seurat list
#'
#' @return A tt seurat object
do_integration_seurat = function(seurat_list) {
	my_features = seurat_list %>% SelectIntegrationFeatures(nfeatures = 3000)

	#plan(strategy = "multicore", workers = 30)

	# Prepare for integration, whatever it means
	prep_integration =
		seurat_list %>%
		PrepSCTIntegration(anchor.features = my_features,
											 verbose = FALSE)

	# Max common cells, in case the common cells are < 30
	dims = map_int(prep_integration, ~ .x %>% ncol) %>% min %>% min(30) %>% sum(-1)

	# Another parameter that is impotant if I have small number of cells
	k.filter <- min(200, min(sapply(prep_integration, ncol)))

	# Integrate
	prep_integration %>%
		FindIntegrationAnchors(
			normalization.method = "SCT",
			anchor.features = my_features,
			verbose = TRUE,
			dims = 1:dims,
			k.filter = k.filter
		) %>%
		IntegrateData(normalization.method = "SCT", verbose = FALSE)
}

#' Get adjusted counts for unwanted variation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom furrr future_map
#' @import Seurat
#'
#' @param .data A seurat list
#' @param .formula A formula
#' @param do.scale A boolean
#' @param do.center A boolean
#'
#' @return A tt seurat object
#'
#' @export
get_adjusted_counts_for_unwanted_variation_sc = function(.data,
																												 .formula,
																												 do.scale = F,
																												 do.center = F,
																												 verbose = T,
																												 ...) {

	# Check if there are column names x or y and stop
	if(.data %>% colnames %in% c("x", "y") %>% any)
		stop("tidysc says: column named \"x\" or \"y\" are banned from seurat metadata because will crash the SCTransform function.")

	# Check if package is installed, otherwise install
	if ("benchmarkme" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing benchmarkme needed for benchmarkme")
		install.packages("benchmarkme", repos = "https://cloud.r-project.org")
	}

	if(Sys.info()[['sysname']] == "Windows")
		available_Gb_ram =
			as.numeric(gsub("\r","",gsub("FreePhysicalMemory=","",system('wmic OS get FreePhysicalMemory /Value',intern=TRUE)[3])))/1024/1024
	else 	available_Gb_ram =  benchmarkme::get_ram() / 10 * 9
	options(future.globals.maxSize = available_Gb_ram * 1000 * 1024 ^ 2)

	# Update on tibble
	.data = .data %>% update_object_sc

	# Sample column name
	#.sample_name = .data %>% attr("parameters")  %$% .sample %>% quo_name

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get integrate column
	.integrate_column = parse_formula(.formula) %>% grep("integrate(", ., fixed = T, value = T) %>% gsub("integrate\\(|\\)", "", .)

	# Get character array of variable to regress
	variables_to_regress = parse_formula(.formula) %>% gsub("integrate\\(|\\)", "", .)


	variables_to_regress_no_sample =
		variables_to_regress %>%
		ifelse_pipe(
			length(.integrate_column) > 0,
			~ .x %>% grep(.integrate_column, ., value = T, invert = T)
		) %>%
		ifelse_pipe( (.) %>% length %>% equals(0), ~ NULL)

	# Check if object normalised already
	# if(
	# 	variables_to_regress_no_sample %>% length %>% `>` (0) &
	# 	"SCT" %in% (.data %>% attr("seurat") %>% `[[` (1) %>% `@` (assays) %>% names)
	# ) stop("The object has been already normalised")

	# Evaluate ...
	arguments <- list(...)

	#plan(strategy = "multicore", workers = 30)

	# Get seurat object
	seurat_object =
		.data %>%
		attr("seurat") %>%

		# If Integration split the object before batch correct
		ifelse_pipe(
			.integrate_column %in% variables_to_regress,
			~ .x %>%
				map(~ .x %>% SplitObject(split.by = "orig.ident")) %>%
				unlist
		) %>%

		# Scale data for covariates other than sample
		map(~ {

			# library(future)
			# plan("multiprocess", workers = 25)

			SCTransform(.x,
									verbose = verbose,
									assay = "RNA",
									vars.to.regress = variables_to_regress_no_sample
									)
					# Needed for the use of ... for a function that has ... already for another thing
			# args = list(
			# 	object = .x,
			# 	verbose = verbose,
			# 	assay = "RNA",
			# 	vars.to.regress = variables_to_regress_no_sample
			# ) %>%
			# c(arguments)
			# exec(SCTransform, !!!args, .env = parent.frame())
			#
			#
			# 		do.call(SCTransform, c(
			# 			.x,
			# 			verbose = verbose,
			# 			assay = "RNA",
			# 			vars.to.regress = variables_to_regress_no_sample,
			# 			arguments
			# 		))
		}
			) %>%

		# INTEGRATION - If sample within covariates Eliminate sample variation with integration
		ifelse_pipe(
			.integrate_column %in% variables_to_regress ,
			~ .x %>%	do_integration_seurat %>%	list
		)

	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (meta.data) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				select(!!.cell, nCount_SCT, nFeature_SCT) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add addjust counts for unwanted variation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @import Seurat
#'
#' @param .data A seurat list
#' @param .formula A formula
#' @param do.scale A boolean
#' @param do.center A boolean
#'
#' @return A tt seurat object
#'
#' @export
add_adjusted_counts_for_unwanted_variation_sc = function(.data,
																												 .formula,
																												 do.scale = F,
																												 do.center = F,
																												 verbose = T,
																												 ...) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_adjusted_counts_for_unwanted_variation_sc(.formula = .formula,
																									do.scale = do.scale,
																									do.center = do.center,
																									verbose = verbose,
																									...)

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Get cluster information of single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#'
#' @param .data A tt object
#'
#' @return A tt tt object
#'
#' @export
get_cluster_annotation_SNN_sc = function(.data, ...) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Add PCA if not present
	.data  = .data %>% add_reduced_dimensions_PCA(.dims = 50)

	# Evaluate ...
	arguments <- list(...)

	# Calculate the new Seurat object
	seurat_object =
		.data %>%
		attr("seurat") %>%
		map(~ {

			# Find neighbours
			my_obj = (.) %>% FindNeighbors()

			# Needed for the use of ... for a function that has ... already for another thing
			my_obj = do.call(FindClusters, c(my_obj, method = "igraph", arguments))

			# Add renamed annotation
			my_obj@meta.data$cluster = my_obj@meta.data$seurat_clusters

			my_obj
		})

	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (meta.data) %>%
				as_tibble(rownames = quo_name(.cell)) %>%
				select(!!.cell, cluster) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add cluster information of single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt tt object
#'
#' @export
add_cluster_annotation_SNN_sc = function(.data, ...) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_cluster_annotation_SNN_sc(...)

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Get cell type information on single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @import SingleR
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom magrittr equals
#'
#' @param .data A tt object
#'
#' @return A tt tt object
#'
#' @export
get_cell_type_annotation_sc = function(.data) {
	# Check if package is installed, otherwise install
	if ("SingleR" %in% rownames(installed.packages()) == FALSE) {
		writeLines("Installing SingleR")
		devtools::install_github('dviraran/SingleR')

	}

	library(SingleR)

	run_singleR = function(seurat,
												 clusters,
												 name = NULL,
												 species,
												 min.transcripts = 200) {
	f = !is.na(clusters)

	# Get references
	hpca.se <- SingleR::HumanPrimaryCellAtlasData()
	blueprint <- SingleR::BlueprintEncodeData()
	MouseRNAseq = MouseRNAseqData()

	# Create ref list
	ref =
		species %>%
			ifelse2_pipe(
				.p1 = (.) %>% tolower() %>% equals("human"),
				.p2 = (.) %>% tolower() %>% equals("mouse"),
				~ list(hpca = hpca.se, blueprint = blueprint),
				~ list(MouseRNAseq = MouseRNAseq),
				stop("Species should be either human or mouse, case insensitive")
			)

	# my assay
	my_assay =
		seurat@assays %>%
		names %>%
		base::grep("integrated", ., invert = TRUE, value=TRUE) %>%
		`[` (length(.))


	# Test data frame
	test = seurat@assays[[my_assay]]@counts %>% `+` (1) %>% log %>% Matrix::Matrix(sparse = TRUE)

	# Calculate and return
	ref %>%
		map2_dfr(
			names(ref),
			~ 	{
				singler = SingleR::SingleR(
					test =  test,
					ref = .x,
					labels = .x$label.main,
					clusters = clusters,
					method = ifelse(clusters %>% is.null, "single", "cluster")
				)

			singler$scores %>% as_tibble %>%
				dplyr::mutate(ID = singler %>% rownames) %>%
				mutate(label = singler$first.labels) %>%
				#gather(ct, score, -ID, -label) %>%
				nest(`Cell type scores` =  -c(ID, label)) %>%
				mutate(SingleR_DB = .y)

			}
		) %>%

	# Integrate inference
	pivot_wider(values_from = c(label, `Cell type scores`), names_from = SingleR_DB)

	}

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Update on tibble
	.data = .data %>% update_object_sc(!!.cell)

	.x = .data %>%
		attr("seurat") %>% `[[` (1)
	.y = .data %>% attr("parameters")

	 	my_clusters = switch(
	 		"cluster" %in% (.x@meta.data %>% colnames) %>% `!` %>% sum(1),
	 		.x@meta.data$cluster %>% as.character,
	 		NULL
	 	)

	 	# Save cell type annotation
	 	ct_class =
	 		.x %>%
	 		run_singleR(my_clusters,
	 								species =  .y$species) %>%

	 		# Change ID name if I have clusters of not
	 		ifelse_pipe(
	 			my_clusters %>% is.null,
	 			~ .x %>% rename(!!.cell := ID),
	 			~ .x %>% rename(cluster = ID)
	 		)

	 	# make dataframe with no nest
	 	ct_class_simple = ct_class %>% select(-contains("Cell type scores"))

	 	# Make writable variable
	 	my_obj = .x

	 	# Add renamed annotation
	 	my_obj@meta.data =
	 		my_obj@meta.data %>%
	 		ifelse_pipe(
	 			my_clusters %>% is.null,

	 			# Bind based on cell name
	 			~ .x %>%
	 				cbind(
	 					ct_class_simple[
	 						match(
	 							.x %>% rownames ,
	 							ct_class_simple %>% pull(!!.cell)),
	 						] %>%
	 						mutate_if(is.character, as.factor) %>%
	 						select(-!!.cell)
	 				),

	 			# Else, bind based on cluster
	 			~ .x %>%
	 				cbind(
	 					ct_class_simple[
	 						match(
	 							.x %>% pull(cluster) ,
	 							ct_class_simple %>% pull(cluster)),
	 						] %>%
	 						mutate_if(is.character, as.factor) %>%
	 						select(-cluster)
	 				)
	 		)

	 	# Return
	 	seurat_object = my_obj %>% list()

	# Make tibble
	seurat_object %>%
		map_dfr(
			~ .x %>%
				`@` (meta.data) %>%
				as_tibble(rownames = quo_name(.cell)) %>%

				# Grub last two columns
				select(!!.cell, one_of("cluster"), (ncol(.) - 1):ncol(.)) %>%
				mutate_if(is.factor, as.character)
		) %>%
		mutate_if(is.character, as.factor) %>%

		# Join back the nested data
		left_join(ct_class ) %>%
		select(-one_of("cluster")) %>%

		# Add back the attributes objects
		add_attr(seurat_object, "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' Add cell type information on single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt tt object
#'
#' @export
add_cell_type_annotation_sc = function(.data) {
	# Update on tibble
	.data = .data %>% update_object_sc

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell

	# Get now object
	.data.annotated =
		.data %>%
		get_cell_type_annotation_sc()

	# Merge
	.data %>%
		left_join(.data.annotated ,
							by = quo_name(.cell)) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}


#' Get rotated dimensions of two principal .dims or MDS dimension of choice, of an angle
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param .data A tibble
#' @param dimension_1_column A column symbol. The column of the dimension 1
#' @param dimension_2_column   A column symbol. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param of_samples A boolean
#' @param dimension_1_column_rotated A column symbol. The column of the dimension 1 rotated
#' @param dimension_2_column_rotated   A column symbol. The column of the dimension 2 rotated
#'
#' @return A tibble with additional rotated columns
#'
#'
get_rotated_dimensions_sc =
	function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 of_samples = T,
					 dimension_1_column_rotated = NULL,
					 dimension_2_column_rotated = NULL) {
		# Get column names
		.element = ifelse(of_samples, .data %>% attr("parameters")  %$% sample)

		# .element = enquo(.element)
		# col_names = get_elements(.data, .element)
		# .element = col_names$.element

		# Parse other colnames
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)

		if (.data %>%
				distinct(!!.element,!!dimension_1_column,!!dimension_2_column) %>%
				count(!!.element,!!dimension_1_column,!!dimension_2_column) %>%
				pull(n) %>%
				max %>%
				`>` (1))
			stop(sprintf(
				"%s must be unique for each row for the calculation of rotation",
				quo_name(.element)
			))

		# Set default col names for rotated dimensions if not set
		if (quo_is_null(dimension_1_column_rotated))
			dimension_1_column_rotated = as.symbol(sprintf(
				"%s rotated %s",
				quo_name(dimension_1_column),
				rotation_degrees
			))
		if (quo_is_null(dimension_2_column_rotated))
			dimension_2_column_rotated = as.symbol(sprintf(
				"%s rotated %s",
				quo_name(dimension_2_column),
				rotation_degrees
			))

		# Function that rotates a 2D space of a arbitrary angle
		rotation = function(m, d) {
			r = d * pi / 180
			((dplyr::bind_rows(
				c(`1` = cos(r), `2` = -sin(r)),
				c(`1` = sin(r), `2` = cos(r))
			) %>% as_matrix) %*% m)
		}

		# Sanity check of the angle selected
		if (rotation_degrees %>% between(-360, 360) %>% `!`)
			stop("rotation_degrees must be between -360 and 360")

		# Return
		.data %>%
			distinct(!!.element,!!dimension_1_column,!!dimension_2_column) %>%
			as_matrix(rownames = !!.element) %>% t %>%
			rotation(rotation_degrees) %>%
			as_tibble() %>%
			mutate(`rotated dimensions` =
						 	c(
						 		quo_name(dimension_1_column_rotated),
						 		quo_name(dimension_2_column_rotated)
						 	)) %>%
			gather(!!.element, value, -`rotated dimensions`) %>%
			spread(`rotated dimensions`, value)

	}

#' Add Rotated dimensions of two principal .dims or MDS dimensions, of an angle
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom rlang quo_is_null
#'
#'
#' @param .data A tibble
#' @param dimension_1_column A column symbol. The column of the dimension 1
#' @param dimension_2_column   A column symbol. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param .element A column symbol. The column that is used to calculate distance (i.e., normally samples)
#' @param of_samples A boolean
#' @param dimension_1_column_rotated A column symbol. The column of the dimension 1 rotated
#' @param dimension_2_column_rotated   A column symbol. The column of the dimension 2 rotated
#'
#' @return A tibble with additional rotated columns
#'
#'
add_rotated_dimensions_sc =
	function(.data,
					 dimension_1_column,
					 dimension_2_column,
					 rotation_degrees,
					 .element = NULL,
					 of_samples = T,
					 dimension_1_column_rotated = NULL,
					 dimension_2_column_rotated = NULL) {
		# Get column names
		.element = enquo(.element)
		col_names = get_elements(.data, .element)
		.element = col_names$.element

		# Parse other colnames
		dimension_1_column = enquo(dimension_1_column)
		dimension_2_column = enquo(dimension_2_column)
		dimension_1_column_rotated = enquo(dimension_1_column_rotated)
		dimension_2_column_rotated = enquo(dimension_2_column_rotated)

		# Set default col names for rotated dimensions if not set
		if (quo_is_null(dimension_1_column_rotated))
			dimension_1_column_rotated = as.symbol(sprintf(
				"%s rotated %s",
				quo_name(dimension_1_column),
				rotation_degrees
			))
		if (quo_is_null(dimension_2_column_rotated))
			dimension_2_column_rotated = as.symbol(sprintf(
				"%s rotated %s",
				quo_name(dimension_2_column),
				rotation_degrees
			))

		.data %>%
			left_join(
				(.) %>%
					get_rotated_dimensions(
						dimension_1_column = !!dimension_1_column,
						dimension_2_column = !!dimension_2_column,
						rotation_degrees = rotation_degrees,
						.element = !!.element,
						of_samples = of_samples,
						dimension_1_column_rotated = !!dimension_1_column_rotated,
						dimension_2_column_rotated = !!dimension_2_column_rotated
					),
				by = quo_name(.element)
			)
	}

#' Add cell type information on single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param .data A tt object
#'
#' @return A tt tt object
#'
#' @export
merged_tt_object = function(...){

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
	)

}

#' @export
mutate_update_and_add_attr = function(.data, ...){

	# Drop class to use dplyr
	class(.data) = class(.data)[-c(1:2)]

	dplyr::mutate(.data, ...) %>%
		add_attr(.data %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters") %>%
		update_metadata_sc() %>%

		# Add tt class
		add_class("tt") %>%
		add_class("tidysc")

}


#' @export
unite_update_and_add_attr = function(.data, col, ..., sep = "_", remove = TRUE, na.rm = FALSE){

	# Drop class to use dplyr
	class(.data) = class(.data)[-c(1:2)]

	tidyr::unite(.data, col, ..., sep = sep, remove = remove, na.rm = na.rm) %>%
		add_attr(.data %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters") %>%
		update_metadata_sc() %>%

		# Add tt class
		add_class("tt") %>%
		add_class("tidysc")

}

filter_update_and_add_attr = function(.data, ...){

	# Drop class to use dplyr
	class(.data) = class(.data)[-c(1:2)]

	dplyr::filter(.data, ...) %>%
		update_object_sc() %>%

		# Add tt class
		add_class("tt") %>%
		add_class("tidysc")

	# %>%
	# 	add_attr(.data %>% attr("seurat"), "seurat") %>%
	# 	add_attr(.data %>% attr("parameters"), "parameters")

}

#' get abundance long
#' @importFrom magrittr "%$%"
#'
#' @export
get_abundance_sc_long = function(.data, transcripts = NULL, all = F){

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell
	.transcript = .data %>% attr("parameters")  %$% .transcript
	.abundance = .data %>% attr("parameters")  %$% .abundance

	# Check if output would be too big without forcing
	if(
		"var.features" %in% (
			.data %>%
					 attr("seurat") %>%
					 `[[` (1) %>%
					 attributes() %>%
					 names
				)  &
		transcripts %>% is.null &
		all == F
	) stop("
				 Your object do not contain variable trancript labels,
				 transcript argument is empty and all argument is set to FALSE.
				 Either:
				 1. use detect_variable_features() to select variable feature
				 2. pass an array of transcripts names
				 3. set all = TRUE (this will output a very large object, do you computer have enough RAM?)
				 ")

	# Get variable features if existing
	if(
		transcripts %>% is.null &
		all == F &
		any(c("normalised", "SCT") %in% (.data %>% attr("seurat") %>% `[[` (1) %>% `@` (assays) %>% names))
	) variable_genes = .data %>% attr("seurat") %>% `[[` (1) %>% `@` (assays) %>% `[[` (length(.)) %>% `@` (var.features)
	else variable_genes = NULL

	assay_names =
		.data %>%
		attr("seurat") %>%
		`[[` (1) %>%
		`@` (assays) %>%
		names


	.data %>%
		attr("seurat") %>%
		`[[` (1) %>%
		`@` (assays) %>%

		# Take active assay
		map2(assay_names,

			~ .x %>%
				ifelse2_pipe(
					variable_genes %>% is.null %>% `!`,
				transcripts %>% is.null %>% `!`,
				~ .x@data[variable_genes,],
				~ .x@data[transcripts[transcripts %in% rownames(.x@data)],],
				~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
			) %>%
				as_tibble(rownames = quo_name(.transcript)) %>%
				gather(!!.cell, !!(quo_name(.abundance) %>% paste(.y, sep="_")), -!!.transcript) %>%
				mutate_if(is.character, as.factor) %>%

				# Add back class
				# Add tt class
				add_class("tt") %>%
				add_class("tidysc") %>%

				# Add back the attributes objects
				add_attr(.data %>% attr("seurat"), "seurat") %>%
				add_attr(.data %>% attr("parameters"), "parameters")
		) %>%
		Reduce(function(...) left_join(..., by=c(quo_name(.transcript), quo_name(.cell))), .)

}

#' @export
add_abundance_sc_long = function(.data, transcripts = NULL, all = F){

	# Update on tibble
	.data = .data %>% update_object_sc()

	# Cell column name
	.cell_name = .data %>% attr("parameters")  %$% .cell %>% quo_name

	# Get now object
	.data.annotated =
		.data %>%
		get_abundance_sc_long(transcripts = transcripts, all = all)

	# Merge
	.data %>%
		left_join(.data.annotated,	by = .cell_name) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

#' get abundance long
#' @importFrom magrittr "%$%"
#'
#' @export
get_abundance_sc_wide = function(.data, transcripts = NULL, all = F){

	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell
	.transcript = .data %>% attr("parameters")  %$% .transcript
	.abundance = .data %>% attr("parameters")  %$% .abundance

	# Check if output would be too big without forcing
	if(
		"var.features" %in% (
			.data %>%
			attr("seurat") %>%
			`[[` (1) %>%
			attributes() %>%
			names
		)  &
		transcripts %>% is.null &
		all == F
	) stop("
				 Your object do not contain variable trancript labels,
				 transcript argument is empty and all argument is set to FALSE.
				 Either:
				 1. use detect_variable_features() to select variable feature
				 2. pass an array of transcripts names
				 3. set all = TRUE (this will output a very large object, do you computer have enough RAM?)
				 ")

	# Get variable features if existing
	if(
		transcripts %>% is.null &
		all == F &
		any(c("normalised", "SCT") %in% (.data %>% attr("seurat") %>% `[[` (1) %>% `@` (assays) %>% names))
	) variable_genes = .data %>% attr("seurat") %>% `[[` (1) %>% `@` (assays) %>% `[[` (length(.)) %>% `@` (var.features)
	else variable_genes = NULL

	# Just grub last assay
	my_assay =
		.data %>%
		attr("seurat") %>%
		`[[` (1) %>%
		`@` (assays) %>%
		`[[` (length(.))

	assay_names =	my_assay %>% names

 my_assay %>%
 	ifelse2_pipe(
 		variable_genes %>% is.null %>% `!`,
 		transcripts %>% is.null %>% `!`,
 		~ .x@counts[variable_genes,],
 		~ .x@counts[transcripts,],
 		~ stop("It is not convenient to extract all genes, you should have either variable features or transcript list to extract")
 	) %>%
 	as.matrix() %>%
 	t %>%
 	as_tibble(rownames = quo_name(.transcript)) %>%
 	mutate_if(is.character, as.factor) %>%

 	# Add back class
 	# Add tt class
 	add_class("tt") %>%
 	add_class("tidysc") %>%

 	# Add back the attributes objects
 	add_attr(.data %>% attr("seurat"), "seurat") %>%
 	add_attr(.data %>% attr("parameters"), "parameters")

}

#' @export
add_abundance_sc_wide= function(.data, transcripts = NULL, all = F){

	# Update on tibble
	.data = .data %>% update_object_sc()

	# Cell column name
	.cell_name = .data %>% attr("parameters")  %$% .cell %>% quo_name

	# Get now object
	.data.annotated =
		.data %>%
		get_abundance_sc_wide(transcripts = transcripts, all = all)

	# Merge
	.data %>%
		left_join(.data.annotated,	by = .cell_name) %>%

		# Add back the attributes objects
		add_attr(.data.annotated %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")

}

