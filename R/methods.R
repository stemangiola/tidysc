setOldClass("tidysc")

#' Creates a `tt` object from a `tbl``
#'
#' \lifecycle{experimental}
#' 
#' @description tidysc_long() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name tidysc
#' @rdname tidysc
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .cell The name of the cell column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @details This function created a tidysc object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The tidysc object have an attribute called parameters where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "parameters").
#'
#' @return A `tidysc` object
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' my_tt =
#'     tidysc::counts %>%
#'     tidysc(sample, transcript, counts)
#'
#' class(my_tt)
#'
#' }
#'
#' @export
tidysc_long <- function(.data,
                      .sample,
                      .cell,
                      .transcript,
                      .abundance,
                      species,
                      min.transcripts = 400,
                      min.cells = 5,
                      ...) {
  UseMethod("tidysc_long", .data)
}
#' @export
tidysc_long.default <- function(.data,
                              .sample,
                              .cell,

                              .transcript,
                              .abundance,
                              species,
                              min.transcripts = 400,
                              min.cells = 5,
                              ...)
{
  print("This function cannot be applied to this object")
}
#' @export
tidysc_long.tbl_df <- function(.data,
                             .sample,
                             .cell,

                             .transcript,
                             .abundance,
                             species,
                             min.transcripts = 400,
                             min.cells = 5,
                             ...)
{
  # Make col names
  .sample = enquo(.sample)
  .cell = enquo(.cell)
  .transcript = enquo(.transcript)
  .abundance = enquo(.abundance)

  create_tt_from_tibble_sc(
    .data,
    .sample = !!.sample,
    .cell = !!.cell,
    .transcript = !!.transcript,
    .abundance = !!.abundance,
    species = species,
    min.transcripts = min.transcripts,
    min.cells = min.cells,
    ...
  )

}

#' Creates a `tt` object from a `tbl``
#'
#' \lifecycle{experimental}
#' 
#' @description tidysc_wide() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name tidysc
#' @rdname tidysc
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .cell The name of the cell column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @details This function created a tidysc object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The tidysc object have an attribute called parameters where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "parameters").
#'
#' @return A `tidysc` object
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' my_tt =
#'     tidysc::counts %>%
#'     tidysc(sample, transcript, counts)
#'
#' class(my_tt)
#'
#' }
#'
#' @export
tidysc_wide <- function(.data,
                      .sample,
                      .cell,
                      .transcript_position,
                      species,
                      min.transcripts = 400,
                      min.cells = 5,
                      ...) {
  UseMethod("tidysc_wide", .data)
}
#' @export
tidysc_wide.default <- function(.data,
                              .sample,
                              .cell,
                              .transcript_position,
                              species,
                              min.transcripts = 400,
                              min.cells = 5,
                              ...)
{
  print("This function cannot be applied to this object")
}
#' @export
tidysc_wide.tbl_df <- function(.data,
                             .sample,
                             .cell,
                             .transcript_position,
                             species,
                             min.transcripts = 400,
                             min.cells = 5,
                             ...)
{
  # Make col names
  .sample = enquo(.sample)
  .cell = enquo(.cell)

  create_tt_from_tibble_wide_sc(
    .data =.data,
    .sample = !!.sample,
    .cell = !!.cell,
    .transcript_position = .transcript_position,
    species = species,
    min.transcripts = min.transcripts,
    min.cells = min.cells,
    genome = ifelse(species == "Human", "hg38", "mm10"),
    ...)

}


#' Creates a `tt` object from a `tbl``
#'
#' \lifecycle{experimental}
#'
#' @description tidysc_cell_ranger() creates a `tt` object from a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name tidysc
#' @rdname tidysc
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @details This function created a tidysc object and is useful if you want
#' to avoid to specify .sample, .transcript and .abundance arguments all the times.
#' The tidysc object have an attribute called parameters where these three
#' arguments are stored as metadata. They can be extracted as attr(<object>, "parameters").
#'
#' @return A `tidysc` object
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' my_tt =
#'     tidysc::counts %>%
#'     tidysc(sample, transcript, counts)
#'
#' class(my_tt)
#'
#' }
#'
#' @export
tidysc_cell_ranger <- function(dir_names,
                             min.transcripts = 400,
                             min.cells = 5,
                             high.mito.thresh = 0.08,
                             high.umi.thresh = 10000,
                             species,
                             genome = ifelse(species == "Human", "hg38", "mm10")) {
  UseMethod("tidysc_cell_ranger", .data)
}
#' #' @export
#' tidysc_cell_ranger.default <- function(dir_names,
#'                                      min.transcripts = 400,
#'                                      min.cells = 5,
#'                                      high.mito.thresh = 0.08,
#'                                      high.umi.thresh = 10000,
#'                                      species,
#'                                      genome = ifelse(species == "Human", "hg38", "mm10"))
#' {
#'   print("This function cannot be applied to this object")
#' }
#' @export
tidysc_cell_ranger.default <- function(dir_names,
                                    min.transcripts = 400,
                                    min.cells = 5,
                                    high.mito.thresh = 0.08,
                                    high.umi.thresh = 10000,
                                    species,
                                    genome = ifelse(species == "Human", "hg38", "mm10"))
{

  create_tt_from_cellRanger_sc(
    dir_names = dir_names,
   min.transcripts = min.transcripts,
   min.cells = min.cells,
   high.mito.thresh = high.mito.thresh,
   high.umi.thresh = high.umi.thresh,
   species = species,
   genome = genome
  )
}

#' Creates a `tt` object 
#'
#' \lifecycle{maturing}
#'
#' @description aggregate_cells() creates a speudo-bulk `tbl` object formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | 
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name aggregate_cells
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#'
#' @details ...
#'
#' @return A `tbl_df` object
#'
#'
#' @examples
#'
#'
#'
#'
#' my_tt =  aggregate_cells(aggregate_cells::counts_mini, sample, transcript, count)
#'
#'
#' @docType methods
#' @rdname aggregate_cells-methods
#' @export
#'
setGeneric("aggregate_cells", function(.data,
                                .sample = NULL, slot = "data", assays = NULL, aggregation_function = Matrix::rowSums
                             )
  standardGeneric("aggregate_cells"))


#' aggregate_cells
#'
#'
#' @inheritParams aggregate_cells
#' @return A `aggregate_cells` object
#'
setMethod("aggregate_cells", "Seurat",  function(.data, .sample = NULL, slot = "data", assays = NULL, aggregation_function = Matrix::rowSums) {

	.sample = enquo(.sample)
	
	# Subset only wanted assays
	if(!is.null(assays)){
		DefaultAssay(.data) = assays[1]
		.data@assays = .data@assays[assays]
	}
	
	.data %>%
	
		tidyseurat::nest(data = -!!.sample) %>%
		mutate(.aggregated_cells = map_int(data, ~ ncol(.x))) %>% 
		mutate(data = map(data, ~ 
				
			# loop over assays
			map2(
				.x@assays, names(.x@assays),
					 
					 # Get counts
					 ~ GetAssayData(.x, slot = slot) %>%
					 	aggregation_function(na.rm = T) %>%
					 	tibble::enframe(
					 		name  = "feature",
					 		value = sprintf("%s", .y)
					 	) %>%
					 	mutate(feature = as.character(feature)) 
			) %>%
			Reduce(function(...) full_join(..., by=c("feature")), .)
			
		)) %>%
		left_join(.data %>% tidyseurat::as_tibble() %>% subset(!!.sample)) %>%
		tidyseurat::unnest(data) %>%
		
		drop_class("tidyseurat_nested")
	
})

#' aggregate_cells
#' 
#' @importFrom tidybulk as_SummarizedExperiment
#' @importFrom tidySingleCellExperiment nest
#' @importFrom tidySingleCellExperiment mutate
#' @importFrom tidySingleCellExperiment as_tibble
#'
#'
#' @inheritParams aggregate_cells
#' @return A `aggregate_cells` object
#'
setMethod("aggregate_cells", "SingleCellExperiment",  function(.data, .sample = NULL, slot = "data", assays = NULL, aggregation_function = Matrix::rowSums) {

	.sample = enquo(.sample)
	
	# Subset only wanted assays
	if(!is.null(assays)){
		.data@assays@data = .data@assays@data[assays]
	}
	
	.data %>%
		
		tidySingleCellExperiment::nest(data = -!!.sample) %>%
		mutate(.aggregated_cells = map_int(data, ~ ncol(.x))) %>% 
		mutate(data = map(data, ~ 
												# loop over assays
												map2(
													as.list(assays(.x)), names(.x@assays),
													
													# Get counts
													~  .x %>%
														aggregation_function(na.rm = T) %>%
														tibble::enframe(
															name  = "feature",
															value = sprintf("%s", .y)
														) %>%
														mutate(feature = as.character(feature)) 
												) %>%
												Reduce(function(...) full_join(..., by=c("feature")), .)
											
		)) %>%
		left_join(.data %>% tidySingleCellExperiment::as_tibble() %>% subset(!!.sample), by = quo_names(.sample)) %>%
		tidySingleCellExperiment::unnest(data) %>%
		
		drop_class("tidySingleCellExperiment_nested") |> 

		as_SummarizedExperiment(.sample = !!.sample, .transcript = feature, .abundance = !!as.symbol(names(.data@assays)))
	
})


#' Normalise the counts of transcripts/genes
#'
#' \lifecycle{experimental}
#'
#' @description scale_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and normalises the data for the library size (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name scale_abundance
#' @rdname scale_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param cpm_threshold A real positive number. It is the threshold of count per million that is used to filter transcripts/genes out from the normalisation procedure. The normalisation inference is then applied back to all unfiltered data.
#' @param prop A real positive number between 0 and 1. It is the threshold of proportion of samples for each transcripts/genes that have to be characterised by a cmp bigger than the threshold to be included for normalisation procedure.
#' @param method A character string. The methods used by the function. The method is passed to the function `calcNormFactors` from edgeR package.
#' @param reference_selection_function A fucntion that is used to selecting the reference sample for normalisation. It could be max (default), which choose the sample with maximum library size; or median, which chooses the sample with median library size.
#' @param action A character string between "add" (default) and "get". "add" joins the new information to the input tbl (default), "get" return a non-redundant tbl with the just new information.
#'
#' @details normalises the data for the library size
#' (e.g., with TMM algorithm, Robinson and Oshlack doi.org/10.1186/gb-2010-11-3-r25).
#' Lowly transcribed transcripts/genes (defined with cpm_threshold and prop parameters)
#' are filtered out from the normalisation procedure.
#' The normalisation inference is then applied back to all unfiltered data.
#'
#' @return A tbl object with additional columns with normalised data as `<NAME OF COUNT COLUMN> normalised`
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' counts %>%
#'     scale_abundance(sample, transcript, `count`)
#'
#'
#'}
#'
#' @export
scale_abundance <- function(.data,
                                verbose = TRUE,
                                action = "add") {
  UseMethod("scale_abundance", .data)
}
#' @export
scale_abundance.default <-  function(.data,
                                         verbose = TRUE,
                                         action = "add")
{
  print("This function cannot be applied to this object")
}
#' @export
scale_abundance.tbl_df = scale_abundance.tidysc <-
  function(.data,
           verbose = TRUE,
           action = "add")
  {
    if (action == "add")
      add_normalised_counts_sc(.data, verbose = verbose)
    else if (action == "get")
      get_normalised_counts_sc(.data, verbose = verbose)
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )

  }

#' @export
scale_abundance.tidyseurat <-
	function(.data,
					 verbose = TRUE,
					 action = "add")
	{

		.data %>%	SCTransform(verbose = verbose,	new.assay.name = "normalised"	)
		
	}
#' Get clusters of elements (e.g., samples or transcripts)
#'
#' \lifecycle{experimental}
#'
#' @description cluster_elements() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and identify clusters in the data. It uses the function Seurat::FindClusters to identify clusters.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name cluster_elements
#' @rdname cluster_elements
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .value The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a tidysc object, it indicates Whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function kmeans

#' @details identifies clusters in the data, normally of samples.
#' This function returns a tibble with additional columns for the cluster annotation.
#' At the moment only k-means clustering is supported, the plan is to introduce more clustering methods.
#'
#' @return A tbl object with additional columns with cluster labels
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' counts %>%
#'     cluster_elements(sample, transcript, count,	centers = 2)
#'
#'     }
#'
#' @export
#'
cluster_elements <- function(.data,
                              action = "add", ...) {
  UseMethod("cluster_elements", .data)
}

#' @export
cluster_elements.default <-  function(.data,
                                       action = "add", ...)
{
  print("This function cannot be applied to this object")
}

#' @export
cluster_elements.tidysc <-  function(.data,
																		 resolution = 0.8,
																		 automatic_resolution = FALSE,
                                    action = "add", ...)
{
	
	# Cell column name
	.cell = .data %>% attr("parameters")  %$% .cell
	
	.data_processed = 
		.data %>%
		update_object_sc %>%
		when(
			automatic_resolution ~ get_cluster_annotation_SNN_automatic(.data, ...),
			~	 get_cluster_annotation_SNN_sc(.data, resolution = resolution, ...)
		) 
	
	if(action=="get"){
		.data_processed
	}
	else if(action=="add"){
	
	# Merge
	.data %>%
		left_join(.data_processed ,	by = quo_name(.cell)) %>%
		
		# Add back the attributes objects
		add_attr(.data_processed %>% attr("seurat"), "seurat") %>%
		add_attr(.data %>% attr("parameters"), "parameters")
	}
  else
    stop(
      "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
    )
}

#' @export
cluster_elements.Seurat <-  function(.data,
																		 resolution = 0.8,
																		 automatic_resolution = FALSE,
																		 action = "add", ...)
{
	
	arguments <- list(...)
	
	.data %>%
		when(
			automatic_resolution ~ IKAP(
						.,
						out.dir = tempfile(pattern = "IKAP_", tmpdir = ".", fileext = ""),
						k.max = 15
					),
			~	RunPCA(., npcs = 15) %>% FindNeighbors() %>% FindClusters(method = "igraph", resolution = resolution )
		) 
	
}


#' Dimension reduction of the transcript abundance data
#'
#' \lifecycle{experimental}
#'
#' @description reduce_dimensions() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and calculates the reduced dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name reduce_dimensions
#' @rdname reduce_dimensions
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .value The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The dimension reduction algorithm to use (PCA, MDS, tSNE).
#' @param top An integer. How many top genes to select for dimensionality reduction
#' @param of_samples A boolean. In case the input is a tidysc object, it indicates Whether the element column will be sample or transcript column
#' @param .dims A list of integer vectors corresponding to principal .dims of interest (e.g., list(1:2, 3:4, 5:6))

#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param scale A boolean for method="PCA", this will be passed to the `prcomp` function. It is not included in the ... argument because although the default for `prcomp` if FALSE, it is advisable to set it as TRUE.
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function prcomp if you choose method="PCA" or Rtsne if you choose method="tSNE"
#'
#' @details This function reduces the dimensions of the transcript abundances.
#' It can use multi-dimensional scaling (MDS) of principal component analysis (PCA).
#'
#' @return A tbl object with additional columns for the reduced dimensions
#'
#'
#' @examples
#' \donttest{
#'
#'
#' library(GGally)
#'
#' counts.MDS =
#'     counts %>%
#'     reduce_dimensions(sample, transcript, count, method="MDS", .dims = 3)
#'
#' counts.MDS %>%
#'     select(contains("Dim"), sample, `Cell type`) %>%
#'     distinct() %>%
#'     GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))
#'
#' counts.PCA =
#'     counts.norm %>%
#'     reduce_dimensions(sample, transcript, count, method="PCA", .dims = 3)
#'
#'counts.PCA %>%
#'    select(contains("PC"), sample, `Cell type`) %>%
#'    distinct() %>%
#'    GGally::ggpairs(columns = 1:3, ggplot2::aes(colour=`Cell type`))
#'
#'}
#'
#' @export
#'
#'
reduce_dimensions <- function(.data,
                              method,
                              .dims = 10,
                              action = "add") {
  UseMethod("reduce_dimensions", .data)
}

#' @export
reduce_dimensions.default <-
  function(.data,
           method,
           .dims = 10,
           action = "add")
  {
    print("This function cannot be applied to this object")
  }

#' @export
reduce_dimensions.tidysc <-
  function(.data,
           method,
           .dims = 10,
           action = "add")
  {
    if (method == "PCA") {
      if (action == "add")
        add_reduced_dimensions_PCA(.data,
                                   .dims = .dims)
      else if (action == "get")
        get_reduced_dimensions_PCA(.data,
                                   .dims = .dims)
      else
        stop(
          "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
        )
    }
    else if (method == "tSNE") {
      if (action == "add")
        add_reduced_dimensions_TSNE(.data,
                                    .dims = .dims)
      else if (action == "get")
        get_reduced_dimensions_TSNE(.data,
                                    .dims = .dims)
      else
        stop(
          "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
        )

    }
    else if (method == "UMAP") {
      if (action == "add")
        add_reduced_dimensions_UMAP(.data,
                                    .dims = .dims)
      else if (action == "get")
        get_reduced_dimensions_UMAP(.data,
                                    .dims = .dims)
      else
        stop(
          "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
        )

    }
    else
      stop("method must be either \"PCA\", \"tSNE\" or \"UMAP\"")

  }
#' @importFrom purrr when
#' @export
reduce_dimensions.tidyseurat <-
	function(.data,
					 method,
					 .dims = 10,
					 action = "add")
	{

		.data  %>%
			FindVariableFeatures(selection.method = "vst") %>%
			RunPCA(npcs = .dims %>% max) %>%
			
			when(
				method == "PCA" ~ (.),
				method == "tSNE" ~ RunTSNE(., reduction = "pca", dims = 1:.dims),
				method == "UMAP" ~ RunUMAP(., reduction = "pca", dims = 1:.dims, n.components = 3L),
				~ stop("method must be either \"PCA\", \"tSNE\" or \"UMAP\"")
			)
		
	}


#' Rotate two dimensions (e.g., principal .dims) of an arbitrary angle
#'
#' \lifecycle{experimental}
#'
#' @description rotate_dimensions() takes as imput a `tbl` formatted as | <DIMENSION 1> | <DIMENSION 2> | <...> | and calculates the rotated dimensional space of the transcript abundance.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name rotate_dimensions
#' @rdname rotate_dimensions
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#'
#' @param dimension_1_column A character string. The column of the dimension 1
#' @param dimension_2_column  A character string. The column of the dimension 2
#' @param rotation_degrees A real number between 0 and 360
#' @param of_samples A boolean. In case the input is a tidysc object, it indicates Whether the element column will be sample or transcript column
#' @param dimension_1_column_rotated A character string. The column of the rotated dimension 1 (optional)
#' @param dimension_2_column_rotated A character string. The column of the rotated dimension 2 (optional)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This function to rotate two dimensions such as the reduced dimensions.
#'
#' @return A tbl object with additional columns for the reduced dimensions. additional columns for the rotated dimensions. The rotated dimensions will be added to the original data set as `<NAME OF DIMENSION> rotated <ANGLE>` by default, or as specified in the input arguments.
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' counts.MDS.rotated =
#'     counts.MDS %>%
#'     rotate_dimensions(`tSNE 1`, `tSNE 2`, rotation_degrees = 45, .element = sample)
#'
#' counts.MDS.rotated %>%
#'     distinct(sample, `tSNE 1`,`tSNE 2`, `Cell type`) %>%
#'     ggplot(aes(x=`tSNE 1`, y=`tSNE 2`, color=`Cell type` )) +
#'     geom_point()
#'}
#'
#' @export
#'
rotate_dimensions <- function(.data,
                              dimension_1_column,
                              dimension_2_column,
                              rotation_degrees,
                              of_samples = T,
                              dimension_1_column_rotated = NULL,
                              dimension_2_column_rotated = NULL,
                              action = "add") {
  UseMethod("rotate_dimensions", .data)
}

#' @export
rotate_dimensions.default <-  function(.data,
                                       dimension_1_column,
                                       dimension_2_column,
                                       rotation_degrees,
                                       of_samples = T,
                                       dimension_1_column_rotated = NULL,
                                       dimension_2_column_rotated = NULL,
                                       action = "add")
{
  print("This function cannot be applied to this object")
}

#' @export
rotate_dimensions.tbl_df = rotate_dimensions.tidysc <-
  function(.data,
           dimension_1_column,
           dimension_2_column,
           rotation_degrees,
           of_samples = T,
           dimension_1_column_rotated = NULL,
           dimension_2_column_rotated = NULL,
           action =
             "add")
  {
    # Make col names
    dimension_1_column = enquo(dimension_1_column)
    dimension_2_column = enquo(dimension_2_column)
    dimension_1_column_rotated = enquo(dimension_1_column_rotated)
    dimension_2_column_rotated = enquo(dimension_2_column_rotated)


    if (action == "add")
      add_rotated_dimensions_sc(
        .data,
        dimension_1_column = !!dimension_1_column,
        dimension_2_column = !!dimension_2_column,
        rotation_degrees = rotation_degrees,
        of_samples = of_samples,
        dimension_1_column_rotated = !!dimension_1_column_rotated,
        dimension_2_column_rotated = !!dimension_2_column_rotated
      )
    else if (action == "get")
      get_rotated_dimensions_sc(
        .data,
        dimension_1_column = !!dimension_1_column,
        dimension_2_column = !!dimension_2_column,
        rotation_degrees = rotation_degrees,
        of_samples = of_samples,
        dimension_1_column_rotated = !!dimension_1_column_rotated,
        dimension_2_column_rotated = !!dimension_2_column_rotated
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }


#' Drop redundant elements (e.g., samples) for which feature (e.g., transcript/gene) aboundances are correlated
#'
#' \lifecycle{experimental}
#'
#' @description remove_redundancy() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | for correlation method or | <DIMENSION 1> | <DIMENSION 2> | <...> | for reduced_dimensions method, and returns a `tbl` with dropped elements (e.g., samples).
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name remove_redundancy
#' @rdname remove_redundancy
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .element The name of the element column (normally samples).
#' @param .feature The name of the feature column (normally transcripts/genes)
#' @param .value The name of the column including the numerical value the clustering is based on (normally transcript abundance)
#'
#' @param method A character string. The cluster algorithm to use, ay the moment k-means is the only algorithm included.
#' @param of_samples A boolean. In case the input is a tidysc object, it indicates Whether the element column will be sample or transcript column
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param correlation_threshold A real number between 0 and 1. For correlation based calculation.
#' @param Dim_a_column A character string. For reduced_dimension based calculation. The column of one principal component
#' @param Dim_b_column A character string. For reduced_dimension based calculation. The column of another principal component
#'
#' @details This function removes redundant elements from the original data set (e.g., samples or transcripts). For example, if we want to define cell-type specific signatures with low sample redundancy. This function returns a tibble with dropped recundant elements (e.g., samples). Two redundancy estimation approaches are supported: (i) removal of highly correlated clusters of elements (keeping a representative) with method="correlation"; (ii) removal of most proximal element pairs in a reduced dimensional space.
#'
#' @return A tbl object with with dropped recundant elements (e.g., samples).
#'
#' @examples
#' \donttest{
#'
#'
#'
#' counts %>%
#'     remove_redundancy(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.value =  count,
#' 	   	method = "correlation"
#' 	   	)
#'
#' counts %>%
#'     remove_redundancy(
#' 	   .element = sample,
#' 	   .feature = transcript,
#' 	   	.value = count,
#' 	   	method = "reduced_dimensions"
#' 	   	)
#'}
#'
#' @export
#'
#'
remove_redundancy <- function(.data,
                           .element = NULL,
                           .feature = NULL,
                           .value,
                           method,

                           of_samples = T,



                           correlation_threshold = 0.9,
                           log_transform = F,

                           Dim_a_column,
                           Dim_b_column) {
  UseMethod("remove_redundancy", .data)
}
#' @export
remove_redundancy.default <-  function(.data,
                                    .element = NULL,
                                    .feature = NULL,
                                    .value,
                                    method,

                                    of_samples = T,



                                    correlation_threshold = 0.9,
                                    log_transform = F,

                                    Dim_a_column,
                                    Dim_b_column)
{
  print("This function cannot be applied to this object")
}
#' @export
remove_redundancy.tbl_df = remove_redundancy.tidysc <-  function(.data,
                                                         .element = NULL,
                                                         .feature = NULL,
                                                         .value,
                                                         method,

                                                         of_samples = T,



                                                         correlation_threshold = 0.9,
                                                         log_transform = F,

                                                         Dim_a_column,
                                                         Dim_b_column)
{
  # Make col names
  .value = enquo(.value)
  .element = enquo(.element)
  .feature = enquo(.feature)

  Dim_a_column = enquo(Dim_a_column)
  Dim_b_column = enquo(Dim_b_column)

  if (method == "correlation")
    remove_redundancy_elements_through_correlation(
      .data,
      .value = !!.value,
      .element = !!.element,
      .feature = !!.feature,
      correlation_threshold = correlation_threshold,
      of_samples = of_samples,
      log_transform = log_transform
    )
  else if (method == "reduced_dimensions")
    remove_redundancy_elements_though_reduced_dimensions(
      .data,
      Dim_a_column = !!Dim_a_column,
      Dim_b_column = !!Dim_b_column,
      .element = !!.element,
      of_samples = of_samples
    )
  else
    stop(
      "method must be either \"correlation\" for dropping correlated elements or \"reduced_dimension\" to drop the closest pair according to two dimensions (e.g., PCA)"
    )

}



#' Adjust transcript abundance for unwanted variation
#'
#' \lifecycle{experimental}
#'
#' @description adjust_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with an edditional adjusted abundance column..
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name adjust_abundance
#' @rdname adjust_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_intrest + batch)
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param log_transform A boolean, whether the value should be log-transformed (e.g., TRUE for RNA sequencing data)
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function sva::ComBat
#'
#' @details This function adjusts the abundance for (known) unwanted variation. At the moment just an unwanted covariated is allowed at a time.
#'
#' @return A `tbl` with additional columns for the adjusted counts as `<COUNT COLUMN> adjusted`
#'
#'
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' adjust_abundance( ~ factor_of_interest + batch )
#'
#'}
#'
#' @export
#'
#'
adjust_abundance <- function(.data,
                             .formula,
                             do.scale = F,
                             do.center = F,
                             verbose = T,
                             action = "add",
														 reference_samples = NULL,
														 assay = counts@active.assay,
														 
                             ...) {
  UseMethod("adjust_abundance", .data)
}
#' @export
adjust_abundance.default <-  function(.data,
                                      .formula,
                                      do.scale = F,
                                      do.center = F,
                                      verbose = T,
                                      action = "add",
																			reference_samples = NULL,
																			assay = counts@active.assay,
																			
                                      ...)
{
  print("This function cannot be applied to this object")
}
#' @export
adjust_abundance.tbl_df = adjust_abundance.tidysc <-
  function(.data,
           .formula,
           do.scale = F,
           do.center = F,
           verbose = T,
           action = "add",
  				 reference_samples = NULL,
  				 assay = counts@active.assay,
  				 
           ...)
  {


    if (action == "add")
      add_adjusted_counts_for_unwanted_variation_sc(
        .data = .data,
        .formula = .formula,
        do.scale = do.scale,
        do.center = do.center,
        verbose = verbose,
        ...
      )
    else if (action == "get")
      get_adjusted_counts_for_unwanted_variation_sc(
        .data = .data,
        .formula = .formula,
        do.scale = do.scale,
        do.center = do.center,
        verbose = verbose,
        ...
      )
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }
#' @export
adjust_abundance.Seurat <-
	function(.data,
					 .formula = NULL,
					 do.scale = F,
					 do.center = F,
					 verbose = T,
					 action = "add",
					 reference_samples = NULL,
					 assay = .data@active.assay,
					 ...)
	{
		
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
		
		
		# # Evaluate ...
		# arguments <- list(...)
		
		.data %>%
		
			# If Integration split the object before batch correct
			when(
				length(.integrate_column) > 0 && .integrate_column %in% variables_to_regress ~ 
					SplitObject(., split.by = .integrate_column),
				~ (.) %>% list
			) %>%
			
			# Scale data for covariates other than sample
			map(~ SCTransform(
				.x,
				verbose = verbose,
				assay = assay,
				vars.to.regress = variables_to_regress_no_sample,
				return.only.var.genes = FALSE
			)) %>%
			
			# INTEGRATION - If sample within covariates Eliminate sample variation with integration
			when(
				length(.integrate_column) > 0 && .integrate_column %in% variables_to_regress  ~ 
					do_integration_seurat(., reference = reference_samples),
				~ (.)[[1]]
			) 
}

#' Aggregates multiple counts from the same samples (e.g., from isoforms), concatenates other character columns, and averages other numeric columns
#'
#' \lifecycle{experimental}
#'
#' @description aggregate_duplicates() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with aggregated transcripts that were duplicated.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name aggregate_duplicates
#' @rdname aggregate_duplicates
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .cell The name of the cell column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param aggregation_function A function for counts aggregation (e.g., sum,  median, or mean)
#' @param keep_integer A boolean. Whether to force the aggregated counts to integer
#'
#' @details This function aggregates duplicated transcripts (e.g., isoforms, ensembl).
#' For example, we often have to convert ensembl symbols to gene/transcript symbol,
#'  but in doing so we have to deal with duplicates. `aggregate_duplicates` takes a tibble
#'  and column names (as symbols; for `sample`, `transcript` and `count`) as arguments and
#'  returns a tibble with aggregate transcript with the same name. All the rest of the column
#'  are appended, and factors and boolean are appended as characters.
#'
#' @return A `tbl` object with aggregated transcript abundance and annotation
#'
#'
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' counts %>%
#'     aggregate_duplicates(
#'     sample,
#'     cell,
#'     transcript,
#'     `count`,
#'     aggregation_function = sum
#'     )
#'
#'}
#'
#' @export
#'
#'
aggregate_duplicates <- function(.data,
                                 aggregation_function = sum,
                                 .sample = NULL,
                                 .cell = NULL,
                                 .transcript = NULL,
                                 .abundance = NULL,
                                 shape = "long",
                                 keep_integer = T) {
  UseMethod("aggregate_duplicates", .data)
}

#' @export
aggregate_duplicates.default <-  function(.data,
                                          aggregation_function = sum,
                                          .sample = NULL,
                                          .cell = NULL,
                                          .transcript = NULL,
                                          .abundance = NULL,
                                          shape = "long",

                                          keep_integer = T)
{
  print("This function cannot be applied to this object")
}

#' @export
aggregate_duplicates.tbl_df = aggregate_duplicates.tidysc <-
  function(.data,
           aggregation_function = sum,
           .sample = NULL,
           .cell = NULL,
           .transcript = NULL,
           .abundance = NULL,
           shape = "long",
           keep_integer = T)  {

    # Make col names
    .sample = enquo(.sample)
    .cell = enquo(.cell)
    .transcript = enquo(.transcript)
    .abundance = enquo(.abundance)

    if (shape == "long")
      aggregate_duplicated_transcripts_sc(
        .data = .data,
        aggregation_function = aggregation_function,
        .sample = !!.sample,
        .cell  = !!.cell,
        .transcript = !!.transcript,
        .abundance = !!.abundance,
        keep_integer = keep_integer
      )
    else if (shape == "wide")
      aggregate_duplicated_transcripts_wide_sc(
        .data,
        aggregation_function = aggregation_function,
        .sample = !!.sample,
        .transcript = !!.transcript,
        keep_integer = keep_integer
      )
    else
      stop(
        "action must be either \"long\" for sample | cell | transcript | abundance data frame, or \"wide\" for sample | cell | gene1 | gene2 | .. data frame "
      )


  }


#' Get cell type proportions from samples
#'
#' \lifecycle{experimental}
#'
#' @description deconvolve_cellularity() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with the estimated cell type abundance for each sample
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name deconvolve_cellularity
#' @rdname deconvolve_cellularity
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#' @param ... Further parameters passed to the function Cibersort
#'
#' @details This function infers the cell type composition of our samples (with the algorithm Cibersort; Newman et al., 10.1038/nmeth.3337).
#'
#' @return A `tbl` object including additional columns for each cell type estimated
#'
#'
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' 	counts %>%
#' 	    deconvolve_cellularity(sample, transcript, `count`)
#'
#'}
#'
#' @export
#'
deconvolve_cellularity <- function(.data,
																	 species,
																	 clusters = NULL,
                               action = "add") {
  UseMethod("deconvolve_cellularity", .data)
}
#' @export
deconvolve_cellularity.default <-  function(.data,
																						species,
																						clusters = NULL,
                                        action = "add")
{
  print("This function cannot be applied to this object")
}
#' @export
deconvolve_cellularity.tbl_df = deconvolve_cellularity.tidysc <-
  function(.data,
  				 species,
  				 clusters = NULL,
           action = "add")  {

    if (action == "add")
      add_cell_type_annotation_sc(.data )
    else if (action == "get")
      get_cell_type_annotation_sc(.data)
    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )
  }

#' @export
deconvolve_cellularity.tidyseurat <-
	function(.data,
					 species,
					 clusters = NULL,
					 action = "add")  {
		
		clusters = enquo(clusters)
		
		# Check if package is installed, otherwise install
		if ("SingleR" %in% rownames(installed.packages()) == FALSE) {
			writeLines("Installing SingleR")
			devtools::install_github('dviraran/SingleR')
		}
		
		library(SingleR)
		
		my_clusters = switch(
			quo_name(clusters) %in% (.data@meta.data %>% colnames) %>% `!` %>% sum(1),
			.data %>% tidyseurat::pull(!!clusters) %>% as.character,
			NULL
		)
		
		# Save cell type annotation
		ct_class =
			.data %>%
			run_singleR(my_clusters, species =  species) %>%
			
			# Change ID name if I have clusters of not
			when(
				my_clusters %>% is.null	~ rename(., cell = ID),
				~ rename(., cluster = ID)
			)
		
		# make dataframe with no nest
		ct_class_simple = ct_class %>% select(-contains("Cell type scores"))
		
		# Add renamed annotation
		.data@meta.data =
			.data@meta.data %>%
			ifelse_pipe(
				my_clusters %>% is.null,
				
				# Bind based on cell name
				~ .x %>%
					cbind(
						ct_class_simple[
							match(
								.x %>% rownames ,
								ct_class_simple %>% pull(cell)),
						] %>%
							mutate_if(is.character, as.factor) %>%
							select(-cell)
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
		
		.data
		

	}

#' Add transcript symbol column from ensembl id
#'
#' \lifecycle{experimental}
#'
#' @description annotate_symbol() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the the additional transcript symbol column
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name annotate_symbol
#' @rdname annotate_symbol
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> |
#' @param .ensembl A character string. The column that is represents ensembl gene id
#'
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details This is useful since different resources use ensembl IDs while others use gene symbol IDs.
#'
#' @return A `tbl` object including additional columns for transcript symbol
#'
#'
#'
#'
#' @examples
#' \donttest{
#'
#'
#'
#' 	counts_ensembl %>% annotate_symbol(ens)
#'
#'}
#'
#' @export
#'
#'
annotate_symbol <- function(.data,
                            .ensembl,
                            action = "add") {
  UseMethod("annotate_symbol", .data)
}
#' @export
annotate_symbol.default <-  function(.data,
                                     .ensembl,
                                     action = "add")
{
  print("This function cannot be applied to this object")
}
#' @export
annotate_symbol.tbl_df = annotate_symbol.tidysc <-
  function(.data,
           .ensembl,
           action = "add")
  {
    # Make col names
    .ensembl = enquo(.ensembl)


    if (action == "add")
      add_symbol_from_ensembl(.data,!!.ensembl)

    else if (action == "get")
      get_symbol_from_ensembl(.data,!!.ensembl)

    else
      stop(
        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
      )

  }



#' Add differential transcription information to a tbl using Seurat
#'
#' \lifecycle{experimental}
#'
#' @description test_differential_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name test_differential_abundance
#' @rdname test_differential_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details At the moment this function uses edgeR only, but other inference algorithms will be added in the near future.
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'\donttest{
#'
#'
#'
#' 	test_differential_abundance(
#' 	    ~ condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'}
#'
#' @export
#'
test_differential_abundance <- function(.data,
                                            .formula,

                                            significance_threshold = 0.05) {
  UseMethod("test_differential_abundance", .data)
}
#' @export
test_differential_abundance.default <-  function(.data,
                                                     .formula,

                                                     significance_threshold = 0.05)
{
  print("This function cannot be applied to this object")
}
#' @export
test_differential_abundance.tbl_df = test_differential_abundance.tidysc <-
  function(.data,
           .formula,

           significance_threshold = 0.05)
  {


      get_differential_transcript_abundance_sc(
        .data,
        .formula,

        significance_threshold = significance_threshold
      )



  }

#' Add differential transcription information to a tbl using edgeR.
#'
#' \lifecycle{experimental}
#'
#' @description extract_abundance() takes as imput a `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> | and returns a `tbl` with additional columns for the statistics from the hypothesis test.
#'
#' @importFrom rlang enquo
#' @importFrom magrittr "%>%"
#'
#' @name extract_abundance
#' @rdname extract_abundance
#'
#' @param .data A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |
#' @param .formula A formula with no response variable, representing the desired linear model
#' @param .sample The name of the sample column
#' @param .transcript The name of the transcript/gene column
#' @param .abundance The name of the transcript/gene abundance column
#'
#' @param significance_threshold A real between 0 and 1 (usually 0.05).
#' @param action A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).
#'
#' @details At the moment this function uses edgeR only, but other inference algorithms will be added in the near future.
#'
#' @return A `tbl` with additional columns for the statistics from the hypothesis test (e.g.,  log fold change, p-value and false discovery rate).
#'
#'
#'
#'
#' @examples
#'\donttest{
#'
#'
#'
#' 	extract_abundance(
#' 	    ~ condition,
#' 	    sample,
#' 	    transcript,
#' 	    `count`
#' 	)
#'
#'}
#'
#' @export
#'
extract_abundance <- function(.data,
                              transcripts = NULL,
                              all = F,
															exclude_zeros = F,
															shape = "long",
                              action = "add") {
  UseMethod("extract_abundance", .data)
}
#' @export
extract_abundance.default <-
  function(.data,
           transcripts = NULL,
           all = F,
  				 exclude_zeros = F,
  				 shape = "long",
           action = "add")
  {
    print("This function cannot be applied to this object")
  }
#' @export
extract_abundance.tidysc <-
  function(.data,
           transcripts = NULL,
           all = F,
  				 exclude_zeros = F,
  				 shape = "long",
  				 action = "add")
  {

  	if(shape == "long"){
	    if (action == "add")
	      add_abundance_sc_long(
	        .data = .data,
	        transcripts = transcripts,
	        all = all,
	        exclude_zeros = exclude_zeros
	      )
	    else if (action == "get")
	    	get_abundance_sc_long(
	        .data = .data,
	        transcripts = transcripts,
	        all = all,
	        exclude_zeros = exclude_zeros
	      )
	    else
	      stop(
	        "action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
	      )
  	}
  	else if(shape == "wide"){
  		if (action == "add")
  			add_abundance_sc_wide(
  				.data = .data,
  				transcripts = transcripts,
  				all = all
  			)
  		else if (action == "get")
  			get_abundance_sc_wide(
  				.data = .data,
  				transcripts = transcripts,
  				all = all
  			)
  		else
  			stop(
  				"action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information"
  			)
  	}
  	else
  		stop(
  			"shape must be either \"long\" or \"wide\" "
  		)
  }

#' @export
#'
score_gene_set <- function(.data, signature = "conservative") {
	UseMethod("score_tcell_exhaustion", .data)
}

#' @importFrom Seurat DefaultAssay
#' @export
score_gene_set.tidyseurat = function(.data, signature = "exhaustion_robust"){
	
	# https://github.com/asmagen/robustSingleCell
	markers_small <- c('Pdcd1', 'Cd244', 'Havcr2', 'Ctla4', 'Cd160', 'Lag3', 'Tigit', 'Cd96')
	
	# https://github.com/sunnyzwu/stromal_subclasses/tree/master/05_gene_signature_analysis
	markers_big = c("PDCD1", "CCL3", "PTPN13", "CASP3", "CD244", "GP49A", "NR4A2", "TRG", "EEA1", "GPD2", "Sep-04", "FASLG", "CD160", "IFIH1", "TANK", "MDFIC", "WBP5", "PTGER4", "SH2D2A", "4631408O11RIK", "NRP1", "ISG20", "GPR56", "CD7", "1110067D22RIK", "CCL4", "GAS2", "ENTPD1", "GPR65", "EOMES", "LITAF", "SERPINA3G", "CCRL2", "PTGER2", "FYN", "NR1I4", "NFATC1", "1810035L17RIK", "RSAD2", "FGL2", "1300007C21RIK", "COCH", "D17H6S56E-5", "PBX3", "LO", "SFRS2IP", "PERP", "LAG3", "SPP1", "TNFRSF1B", "LOC381765", "BC039093", "CD200", "CXCL10", "H2-T23", "BCL2A1B", "IIGP1", "AHR", "2010100O12RIK", "C330007P06RIK", "RNF11", "LYCAT", "ADAM19", "CD9", "ART3", "LILRB4", "TRIM17", "5830471E12RIK", "HIST3H2A", "IFI204", "PGLYRP1", "GM1066", "IFIT1", "TCRB-J", "TNFRSF7", "GZMK", "MLLT3", "RGS16", "IFIT3", "CBX6", "ITM2C", "ALCAM", "BCL2A1A", "IRF8", "RBL2", "ITGAV", "CASP1", "CTLA4", "EIF3S1", "SPG21", "CAPZB", "ACADL", "MX1", "NDUFA5", "RASA1", "SERPINB6A", "RCN1", "TNFSF11", "HCCS", "CCR5", "2310015N07RIK", "CPEB2", "NUCB1")
	
	# Choose my signature
	my_markers = 
		signature %>%
		when(
			. == "exhaustion_blackburn" ~ markers_big,
			. == "exhaustion_robust" ~ markers_small,
			~ signature
		)
	
	# Get counts
	my_counts = .data@assays[[DefaultAssay(.data)]]@data
	
	# Calculate mean score
	score_df = 
		my_counts[tolower(rownames(my_counts)) %in% tolower(my_markers),] %>%
		as.matrix() %>%
		
		# If is exp log it
		when(max(.) > 10 ~ log1p(.), ~ (.)) %>%
		colMeans() %>%
		enframe(name = "cell", value = signature)
	
	# Integrate
	.data %>%
		tidyseurat::left_join(score_df, by = "cell")
	
}

#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import purrr
#' @import DropletUtils
#' @import EnsDb.Hsapiens.v86
#' @import AnnotationDbi
#' @import tidySingleCellExperiment
#' @import scuttle
#' @importFrom SingleCellExperiment counts
#' 
#' 
#' @export
drop_empty_and_dead = function(.data){
	
	.data = .data %>% as.SingleCellExperiment()
	
	barcode_ranks = 
		.data %>%
		nest(sce = -sample) %>% 
		mutate(barcodeRanks = map(
			sce, 
			~ {
				
				.x %>% counts() %>% DropletUtils::barcodeRanks()
			}
			
		)) %>% 
		mutate(
			barcodeRanks_df = map(
				barcodeRanks,
				~ .x %>% as_tibble(rownames = "cell") %>% mutate(
					knee =  metadata(.x)$knee,
					inflection =  metadata(.x)$inflection
				)
			)
		) %>%
		
		dplyr::select(-sce, -barcodeRanks) %>% 
		unnest(barcodeRanks_df) %>%
		arrange(rank) 
	
	# # Plot bar-codes ranks
	# plot_barcode_ranks =
	# 	barcode_ranks%>%
	# 	sample_frac(0.005) %>%
	# 	ggplot2::ggplot(aes(rank, total)) +
	# 	geom_point() +
	# 	facet_wrap(~sample, nrow=2) +
	# 	geom_line(aes(rank, fitted), color="red") +
	# 	geom_hline(aes(yintercept = knee), color="dodgerblue") +
	# 	geom_hline(aes(yintercept = inflection), color="forestgreen") +
	# 	scale_x_log10() +
	# 	scale_y_log10()
	
	counts_real_cells = 
		.data %>%
		inner_join(
			barcode_ranks %>%
				tidySingleCellExperiment::filter(total>inflection),
			by = c("cell", "sample")
		)
	
	drop_dead(counts_real_cells)


}

#' @export
drop_dead = function(counts_real_cells){
	
	is_seurat = is(counts_real_cells, "Seurat")
	
	if(is_seurat) counts_real_cells = as.SingleCellExperiment(counts_real_cells)
	
	which_genes_are_mitochondrion = 
		rownames(counts_real_cells) %>% 
		when(
			length(grep("^MT-", .)) > 1 ~ grep("^MT-", .),
			~ {
				# Annotate with mito - Gene-product location
				location <- mapIds(
					EnsDb.Hsapiens.v86, 
					keys=rownames(counts_real_cells), 
					column="SEQNAME",
					keytype="GENEID"
				)
				
				which(location=="MT")
			}
		)
	
	
	counts_mito = 
		counts_real_cells %>%
		
		nest(data = -sample) %>%
		mutate(data = map(
			data, 
			~ .x %>%
				
				# Join mitochondrion statistics
				left_join(
					scuttle::perCellQCMetrics(., subsets=list(Mito=which_genes_are_mitochondrion)) %>%
						as_tibble(rownames="cell"),
					by="cell"
				) %>%
				
				# Label cells
				tidySingleCellExperiment::mutate(high_mitochondrion = isOutlier(subsets_Mito_percent, type="higher")) %>%
				
				# Join mitochondrion statistics
				tidySingleCellExperiment::mutate(mito_RPS = PercentageFeatureSet(
					as.Seurat(., data = NULL) %>% 
						tidyseurat::tidy() %>% 
						tidyseurat::mutate(nCount_RNA = sum), 
					pattern = "^RPS|^RPL")[,1]
				) %>%
				
				# Label cells
				tidySingleCellExperiment::mutate(high_RPS = isOutlier(mito_RPS, type="higher"))
			
		)) %>%
		unnest(data) 
	
	return = 
		counts_mito %>% 
		tidySingleCellExperiment::filter(!high_mitochondrion) 
	
	
	if(is_seurat) return = as.Seurat(return)
	
	return
	
}

