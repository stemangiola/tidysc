# library(furrr)
# plan(strategy = "multicore", workers = 20)
# options(future.globals.maxSize = 8000 * 1024 ^ 2)



#' Create tt object from seurat object
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map_dfr
#' @importFrom purrr map
#'
#' @param seurat_object A seurat object
#' @param dir_names A character array with the output directories of cellRanger if run from internal function
#' @param min.transcripts An integer with the threshold of the minimum number of genes allowed for each cell
#' @param min.cells An integer with the threshold of the minimum anumber of cells allowed for each sample
#' @param high.mito.thresh An numeric with the threshold of the maximum fraction of reads from mitochondrion. Used for filtering
#' @param high.umi.thresh An numeric with the threshold of the maximum fraction of umis. Used for filtering
#' @param genome A character name of the mapping genome used
#' @param sample_column A symbol for the sample column
#' @param cell_column A symbol for the cell column
#' @param species A character name of the species
#'
#' @return A tt object
create_tt_from_seurat = function(seurat_object,
                                 dir_names = list(""),
                                 min.transcripts = 400,
                                 min.cells = 5,
                                 high.mito.thresh = 0.08,
                                 high.umi.thresh = 10000,
                                 sample_column = `sample`,
                                 cell_column = `cell`,
                                 species,
                                 genome = ifelse(species == "Human", "hg38", "mm10")) {
  writeLines("Converting Seurat object back to tibble")

  # Parse column names
  sample_column = enquo(sample_column)
  cell_column = enquo(cell_column)

  # Create object
  seurat_object %>%
    map_dfr(
      ~ .x@meta.data %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        mutate_if(is.factor, as.character) %>%
        mutate_if(is.numeric, as.integer)
    )  %>%
    mutate_if(is.character, as.factor) %>%
    rename(`read count total` = nCount_RNA,
           `gene count` = nFeature_RNA) %>%

    # Pass the sample column instead of origin.ident
    select(-one_of(quo_name(sample_column))) %>%
    mutate(!!sample_column := `orig.ident`) %>%
    select(-`orig.ident`) %>%

    # Reorde columns logically
    select(!!sample_column,!!cell_column, everything()) %>%

    # Add Seurat object
    add_attr(seurat_object, "seurat") %>%

    # Add parameters object
    add_attr(map(dir_names,
                 ~ {
                   list(
                     min.transcripts = min.transcripts,
                     min.cells = min.cells,
                     high.mito.thresh = high.mito.thresh,
                     high.umi.thresh = high.umi.thresh,
                     output.dir = file.path("Results", .x %>% basename()),
                     input.dir = .x,
                     genome = genome,
                     species = species
                   )
                 }),

             "parameters") %>%

    # Add mitochondrion information
    add_mitochndrion_transcription_abundance_sc(cell_column = !!cell_column) %>%

    # Filter dead cells
    filter(`read count total` > 200 & `mito.fraction` < 0.1) %>%
  	update_object_sc %>%

    # Add cell cycle information
    add_cell_cycle_annotation_sc(cell_column = !!cell_column) %>%

    # Normalise for observation
    add_attr((.) %>% attr("seurat") %>% map(~ .x %>% SCTransform(verbose = TRUE)),
             "seurat") %>%

    # Add tt class
    add_class("tt")

}

#' Create tt object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param input.df A tibble
#' @param sample_column A symbol for the sample column
#' @param cell_column A symbol for the cell column
#' @param transcript_column A symbol for the transcript name column
#' @param counts_column A symbol for read counts
#'
#' @return A tt object
create_seurat_from_tibble = function(input.df,
                                     sample_column,
                                     cell_column,
                                     transcript_column,
                                     counts_column,
                                     min.transcripts = 400,
                                     min.cells = 5,
                                     ...) {
  writeLines("Creating seurat object")

  # Prepare column name enquo
  sample_column = enquo(sample_column)
  cell_column = enquo(cell_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  # Check that there are not multiple information for each cell type, that each cell type is not repeted row-wise
  if (input.df %>%
      select(-!!transcript_column,-!!counts_column) %>%
      distinct() %>%
      count(!!cell_column) %>%
      pull(n) %>%
      max %>%
      `>` (1))
    stop(
      "The annotations in the input dataframe must be in one record for each cell.
      That is,
      YOUR_DATA %>%
      select(transcript_column, counts_column) %>%
      distinct() %>%
      count(!!cell_column) %>%
      pull(n) %>%
      max %>% `==` (1)"
    )

  # Set sample names
  sample_names =
    input.df %>%
    select(-!!transcript_column,-!!counts_column) %>%
    distinct() %>%
    pull(!!sample_column)

  # Create Seurat object
  seurat_obj =
    CreateSeuratObject(
      # Read counts
      input.df %>%
        distinct(!!cell_column,!!transcript_column,!!counts_column) %>%
        spread(!!cell_column,!!counts_column) %>%
        data.frame(row.names = quo_name(transcript_column)),

      # Sample/cell information
      meta.data =
        input.df %>%
        select(-!!transcript_column,-!!counts_column) %>%
        distinct() %>%
        data.frame(row.names = quo_name(cell_column)),

      # Other parameters
      min.cells = min.cells,
      min.features = min.transcripts

    )

  # pass the sample column to the orig.ident column of meta data
  seurat_obj@meta.data$orig.ident =
    colnames(seurat_obj) %>%
    as_tibble() %>%
    rename(cell = value) %>%
    left_join(
      input.df %>%
        distinct(!!sample_column,!!cell_column),
      by = "cell"
    )  %>%
    pull(!!sample_column) # seurat_obj@meta.data[[quo_name(sample_column)]]

  # Return
  seurat_obj

}

#' Create tt object from tibble
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom purrr map
#'
#' @param input.df A tibble
#' @param sample_column A character name of the sample column
#' @param cell_column A character name of the cell name column
#' @param transcript_column A character name of the transcript name column
#' @param counts_column A character name of the read count column
#' @param species A character name of the species used (e.g., Mouse, Human)
#'
#' @return A tibble with an additional column
#'
#' @export
create_tt_from_tibble_sc = function(input.df,
                                    sample_column,
                                    cell_column,
                                    transcript_column,
                                    counts_column,
                                    species,
                                    min.transcripts = 400,
                                    min.cells = 5,
                                    ...) {
  writeLines("Start parsing the data frame")

  # Check is species is of the right type
  if (!(species %in% c("Human", "Mouse")))
    stop("Species must be \"Human\" or \"Mouse\"")

  # Prepare column name enquo
  sample_column = enquo(sample_column)
  cell_column = enquo(cell_column)
  transcript_column = enquo(transcript_column)
  counts_column = enquo(counts_column)

  convert_underscore_to_dash = function(input.df, transcript_column){

    transcript_column = enquo(transcript_column)

    wired_names =
      input.df %>%
      distinct(!!transcript_column) %>%
      filter(grepl("_", !!transcript_column)) %>%
      pull(!!transcript_column)

    input.df %>%

      # Check if wired names are present
      ifelse_pipe(
        wired_names %>% length %>% `>` (0),
        ~ {
          warning("Genes with possibly a strange name are present %s. Converting \"_\" to \"-\" ", paste(wired_names, collapse=", "))

          bind_rows(

            # Genes with normal names
            .x %>%
              filter(!(!!transcript_column %in% wired_names)),

            # Genes with wired names converted
            .x %>%
              filter(!(!!transcript_column %in% wired_names)) %>%
              mutate(!!transcript_column := gsub("_", "-", !!transcript_column))
          )
        }
      )
  }

  # Create seurat object
  input.df %>%

    # Eliminate gene statistics from genes and convert genes with wired name
    filter(!(!!transcript_column %in% c("no_feature", "too_low_aQual", "not_aligned", "alignment_not_unique"))) %>%
    convert_underscore_to_dash(!!transcript_column) %>%

    # Check if any cell name starts with a number
    # If so put a _ before, save original name
    # and through warning
    ifelse_pipe((.) %>% filter(grepl("^[0-9]",!!cell_column)) %>% nrow %>% `>` (1),
                ~ {
                  warning(
                    "
                    some cell names started with a number.
                    This is incompatible with Seurat object, the character _ had been added as prefix.
                    The orignal cell name as been kept as column in the data frame"
                  )
                  .x %>%
                    mutate(!!as.symbol(sprintf("%s original", quo_name(cell_column))) := !!cell_column) %>%
                    mutate(!!cell_column := !!cell_column %>% paste0("X_", .))
                }) %>%

    # Convert characters to factors
    mutate_if(is.character, as.factor) %>%

    # If there is more than one sample add number to cell names
    ifelse_pipe(
      (.) %>% pull(!!sample_column) %>% levels %>% length %>% `>` (1),
      ~ .x %>%
        mutate(sample_idx = !!sample_column %>% as.integer) %>%
        unite(!!cell_column, c(!!cell_column, sample_idx), sep = "_")
    ) %>%

    # Anonymous function - create Seurat object
    # input: tibble
    # output Seurat object
    create_seurat_from_tibble(!!sample_column,
                              !!cell_column,
                              !!transcript_column,
                              !!counts_column,
                              min.transcripts = min.transcripts,
                              min.cells = min.cells,
                              ...) %>%

    # Create tt object from seurat
    list %>%
    create_tt_from_seurat(
      sample_column = !!sample_column,
      cell_column = !!cell_column,
      species = species
    ) %>%

    # Eliminate orig.ident column, because same as sample
    select(-contains("orig.ident")) %>%

    # Add parameters attribute
    add_attr(map((.) %>% attr("parameters"),
                 ~ .x %>% c(
                   list(
                     sample_column = enquo(sample_column),
                     cell_column = enquo(cell_column),
                     transcript_column = enquo(transcript_column),
                     counts_column = enquo(counts_column),
                     species = species
                   )
                 )),
             "parameters")

}

#' Create tt object from cellRanger results
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr pmap
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
             sprintf("%s/outs/filtered_feature_bc_matrix", ..1) %>%
             Read10X() %>%
             CreateSeuratObject(
               min.cells = min.cells,
               min.features = min.transcripts,
               project = ..1 %>% basename()
             )

           # If there is more than one sample add number to cell names
           if (..3 %>% `>` (1))
             my_obj = my_obj %>% RenameCells((.) %>% colnames %>% paste(..2, sep = "_"))

           my_obj

         })

  # If object empty through error
  if (seurat_object %>% length == 0)
    stop(
      "The directory provided seem to be empty. Are you sure you are in the right working directory?"
    )

  # create_tt_from_seurat
  seurat_object %>%
    create_tt_from_seurat(species = species)

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
  # 			 			 		gather(sample, `read count`,-transcript) %>%
  # 			 			 		mutate(`read count` = `read count` %>% as.integer) %>%
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
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
get_mitochndrion_transcription_abundance_sc = function(input.df, cell_column) {
  writeLines("Calculating mitochondrion trancription")

  # Set column names
  cell_column = enquo(cell_column)

  # Cell column name
  # cell_column_name = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column %>% quo_name

  # Update on tibble
  input.df = input.df %>% update_object_sc(!!cell_column)

  seurat_obj =
    input.df %>%
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
        enframe(name = quo_name(cell_column), value = "read count mitochondrion total") %>%
        mutate(`read count mitochondrion total` = `read count mitochondrion total` %>% as.integer) %>%
        #mutate(sample = .x@project.name) %>%
        left_join(
          input.df %>%
            select(cell, `read count total`) %>%
            mutate_if(is.factor, as.character),
          by = quo_name(cell_column)
        ) %>%
        mutate(`fraction mitochondrion` = `read count mitochondrion total` / `read count total`) %>%
        select(cell,
               `read count mitochondrion total`,
               `fraction mitochondrion`)

      .x %>%
        AddMetaData(metadata = mito %>% pull(`fraction mitochondrion`),
                    col.name = "mito.fraction") %>%
        AddMetaData(
          metadata = mito %>% pull(`read count mitochondrion total`),
          col.name = "mito.tot"
        )
    })

  seurat_obj %>%
    map_dfr(~ {
      # Select mitochondrion transcription
      .x %>%
        `@` (meta.data) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        #mutate(sample = .x@project.name) %>%
        # rename(
        #   `read count mitochondrion total` = mito.tot,
        #   `fraction mitochondrion` = mito.fraction
        # ) %>%
        select(cell, mito.fraction, mito.tot)
    }) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add mitochondrion transcription abundance
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
add_mitochndrion_transcription_abundance_sc = function(input.df, cell_column) {
  # Get column names
  cell_column = enquo(cell_column)

  # Update on tibble
  input.df = input.df %>% update_object_sc(!!cell_column)

  # Cell column name
  # cell_column_name = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column %>% quo_name

  # Get now object
  input.df.annotated =
    input.df %>%
    get_mitochndrion_transcription_abundance_sc(cell_column = !!cell_column)

  # Merge
  input.df %>%
    left_join(input.df.annotated  ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Get normalised counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
get_normalised_counts_sc = function(input.df) {
  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column
  sample_column = input.df %>% attr("parameters") %>% `[[` (1) %$% sample_column
  transcript_column = input.df %>% attr("parameters") %>% `[[` (1) %$% transcript_column
  counts_column = input.df %>% attr("parameters") %>% `[[` (1) %$% counts_column

  # Update on tibble
  input.df = input.df %>% update_object_sc(!!cell_column)

  seurat_object =
    input.df %>%
    attr("seurat") %>%
    map(~ {
      .x %>%

        # Add raw counts assay
        # add_assay((.)@assays$raw_counts, "normalised_counts") %>%
        # set_default_assay_seurat("normalised_counts") %>%

        NormalizeData(normalization.method = "LogNormalize",
                      scale.factor = 10000)

    })

  # Get normalised counts
  seurat_object %>%
    map_dfr(
      ~ .x@meta.data %>%
        as_tibble(rownames = quo_name(cell_column)) %>%

        mutate_if(is.factor, as.character) %>%
        mutate_if(is.numeric, as.integer)
    )  %>%
    select(-nFeature_RNA) %>%
    mutate_if(is.character, as.factor) %>%
    rename(!!sample_column := orig.ident,
           `read count normalised total` = nCount_RNA) %>%
    select(!!sample_column, everything()) %>%

    # Add back the attributes objects
    add_attr(seurat_object, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")
}

#' Add normalised counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
add_normalised_counts_sc = function(input.df) {
  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Update on tibble
  input.df = input.df %>% update_object_sc(!!cell_column)

  input.df.normalised =
    input.df %>%
    get_normalised_counts_sc()

  input.df %>%
    left_join(input.df.normalised,  by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.normalised %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add doublet classification
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr pmap
#'
#' @param input.df A tt object
#' @param doublet.reps A numeric
#'
#' @return A tt object
#'
#' @export
add_doublet_classification_sc = function(input.df, doublet.reps) {
  input.df %>%
    mutate(`seurat object` =
             pmap(list(`seurat object`, dir_names, parameters), ~ {
               seurat_obj = ..1
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

               seurat_obj$doublet =
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

               seurat_obj$doublet = seurat_obj$doublet[match(colnames(seurat_obj), barcodes)]

               # seurat_obj@meta.data =
               # 	seurat_obj@meta.data
               # %>%
               # 	mutate(barcodes = colnames(..1)) %>%
               # 	left_join(
               # 			seurat_obj$doublet %>%
               # 			enframe %>%
               # 			rename(barcodes = name, doublet = value),
               # 		by = "barcodes"
               # 	) %>%
               #
               # 	select(-barcodes)
               # return

               seurat_obj
             }))
}

#' Add variable gene annotation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
add_variable_genes_classification = function(input.df) {
  # Update on tibble
  input.df = input.df %>% update_object_sc()

  input.df %>%
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
#'
#' @param input.df A seurat object
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
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
get_cell_cycle_annotation_sc = function(input.df, cell_column) {
  # Get column names
  cell_column = enquo(cell_column)

  # Progress
  writeLines("Classifying cells among cell cycle states")

  # Cell column name
  # cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column
  # sample_column = input.df %>% attr("parameters") %>% `[[` (1) %$% sample_column

  # Update on tibble
  input.df = input.df %>% update_object_sc(!!cell_column)

  # Update Seurat object
  seurat_obj =

    input.df %>%
    attr("seurat") %>%
    map2(input.df %>% attr("parameters"),
         ~ {
           # Get cell cycle genes
           cc.genes =
             switch(
               grepl("mm", .y$genome) %>% `!` %>% sum(1),
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
                 .y %$% species
               )
             )


           # Return
           iterate_cell_cycle_scoreing(.x,
                                       s.features = cc.genes$s.genes,
                                       g2m.features = cc.genes$g2m.genes)


         })

  # Convert Seurat object to tibble
  seurat_obj %>%
    map2_dfr(
      input.df %>% attr("parameters"),
      ~ .x %>%
        `@` (meta.data) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%

        # Function that I have to sort out better
        select(-`orig.ident`) %>%

        # Prepare the data frame
        select(!!cell_column, S.Score, G2M.Score, Phase) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add cell cyle annotation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#'
#' @return A tt object
#'
#' @export
add_cell_cycle_annotation_sc = function(input.df, cell_column) {
  # Get column names
  cell_column = enquo(cell_column)


  # Update on tibble
  input.df = input.df %>% update_object_sc(!!cell_column)

  # Cell column name
  # cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_cell_cycle_annotation_sc(cell_column = !!cell_column)

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Get principal components
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#'
#' @param input.df A tt object
#' @param components An integer, the number of components to calculate
#'
#' @return A tt object
#'
#' @export
get_reduced_dimensions_PCA = function(input.df,
                                      components = 10) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  seurat_obj =

    input.df %>%
    add_variable_genes_classification() %>%
    attr("seurat") %>%

    # Scale data for PCA
    map(~
          .x %>%
          #ScaleData(display.progress = T, num.cores=4, do.par = TRUE) %>%
          RunPCA(npcs = components %>% max))


  seurat_obj %>%
    map_dfr(
      ~ .x %>%
        `@` (reductions) %$% pca %>% `@` (cell.embeddings) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        setNames((.) %>% colnames %>% gsub("_", " ", .)) %>%
        #mutate(sample = .x@project.name) %>%
        select(!!cell_column, everything()) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add principal components
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#' @param components An integer, the number of components to calculate
#'
#' @return A tt object
#'
#' @export
add_reduced_dimensions_PCA = function(input.df, components = 10) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_reduced_dimensions_PCA(components)

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Get UMAP dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#'
#' @param input.df A tt object
#' @param components An integer, the number of components to calculate
#'
#' @return A tt object
#'
#' @export
get_reduced_dimensions_UMAP = function(input.df, components = 10) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  seurat_obj =

    input.df %>%
    add_variable_genes_classification() %>%
    add_reduced_dimensions_PCA(components) %>%

    attr("seurat") %>%

    # Scale data for UMAP
    map(~
          .x %>%
          #ScaleData(display.progress = T, num.cores=4, do.par = TRUE) %>%
          RunUMAP(reduction = "pca", dims = 1:components))

  seurat_obj %>%
    map_dfr(
      ~ .x %>%
        `@` (reductions) %$% umap %>% `@` (cell.embeddings) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        setNames((.) %>% colnames %>% gsub("_", " ", .)) %>%
        #mutate(sample = .x@project.name) %>%
        select(!!cell_column, everything()) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add UMAP dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#' @param components An integer, the number of components to calculate
#'
#' @return A tt object
#'
#' @export
add_reduced_dimensions_UMAP = function(input.df, components = 10) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_reduced_dimensions_UMAP(components = components)

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Get TSNE dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#'
#' @param input.df A tt object
#' @param components An integer, the number of components to calculate
#'
#' @return A tt object
#'
#' @export
get_reduced_dimensions_TSNE = function(input.df, components = 10) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  seurat_obj =

    input.df %>%
    add_variable_genes_classification() %>%
    add_reduced_dimensions_PCA(components) %>%

    attr("seurat") %>%

    # Scale data for UMAP
    map(~
          .x %>%
          #ScaleData(display.progress = T, num.cores=4, do.par = TRUE) %>%
          RunTSNE(reduction = "pca", dims = 1:components))

  seurat_obj %>%
    map_dfr(
      ~ .x %>%
        `@` (reductions) %$% tsne %>% `@` (cell.embeddings) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        setNames((.) %>% colnames %>% gsub("_", " ", .)) %>%
        #mutate(sample = .x@project.name) %>%
        select(!!cell_column, everything()) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add TSNE dimensions
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#' @param components An integer, the number of components to calculate
#'
#' @return A tt object
#'
#' @export
add_reduced_dimensions_TSNE = function(input.df, components = 10) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_reduced_dimensions_TSNE(components = components)

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Integrate many seurat objects
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map_int
#'
#' @param input.df A seurat list
#'
#' @return A tt seurat object
do_integration_seurat = function(seurat_list) {
  my_features = seurat_list %>% SelectIntegrationFeatures(nfeatures = 3000)

  # Prepare for integration, whatever it means
  prep_integration =
    seurat_list %>%
    PrepSCTIntegration(anchor.features = my_features,
                       verbose = FALSE,
                       assay =)

  # Max common cells, in case the common cells are < 30
  dims = map_int(prep_integration, ~ .x %>% ncol) %>% min %>% min(30) %>% sum(-1)

  # Another parameter that is impotant if I have small number of cells
  k.filter <- min(200, min(sapply(prep_integration, ncol)))

  # Integrate
  prep_integration %>%
    FindIntegrationAnchors(
      normalization.method = "SCT",
      anchor.features = my_features,
      verbose = FALSE,
      dims = 1:dims,
      k.filter = k.filter
    ) %>%
    IntegrateData(normalization.method = "SCT", verbose = FALSE)
}

#' Get adjusted read counts for unwanted variation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#'
#' @param input.df A seurat list
#' @param formula A formula
#' @param do.scale A boolean
#' @param do.center A boolean
#'
#' @return A tt seurat object
#'
#' @export
get_adjusted_counts_for_unwanted_variation_sc = function(input.df,
                                                         formula,
                                                         do.scale = F,
                                                         do.center = F) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Sample column name
  sample_column_name = input.df %>% attr("parameters") %>% `[[` (1) %$% sample_column %>% quo_name

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get character array of variable to regress
  variables_to_regress = parse_formula(formula)

  # Get seurat object
  seurat_obj =
    input.df %>%
    attr("seurat") %>%

    # If Integration split the object before batch correct
    ifelse_pipe(
      sample_column_name %in% variables_to_regress,
      ~ .x %>%
        map(~ .x %>% SplitObject(split.by = "orig.ident")) %>%
        unlist
    ) %>%

    # Scale data for covariates other than sample
    map(~
          .x %>%
          SCTransform(
            verbose = TRUE,
            vars.to.regress = variables_to_regress %>% grep(sample_column_name, ., value = T, invert = T)
          )) %>%

    # INTEGRATION - If sample within covariates Eliminate sample variation with integration
    ifelse_pipe(
      sample_column_name %in% variables_to_regress ,
      ~ .x %>%
        do_integration_seurat %>%
        list
    )

  seurat_obj %>%
    map_dfr(
      ~ .x %>%
        `@` (meta.data) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        select(!!cell_column, nCount_SCT, nFeature_SCT) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add addjust read counts for unwanted variation
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map
#'
#' @param input.df A seurat list
#' @param formula A formula
#' @param do.scale A boolean
#' @param do.center A boolean
#'
#' @return A tt seurat object
#'
#' @export
add_adjusted_counts_for_unwanted_variation_sc = function(input.df,
                                                         formula,
                                                         do.scale = F,
                                                         do.center = F) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_adjusted_counts_for_unwanted_variation_sc(formula = formula,
                                                  do.scale = do.scale,
                                                  do.center = do.center)

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

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
#' @param input.df A tt object
#'
#' @return A tt tt object
#'
#' @export
get_cluster_annotation_SNN_sc = function(input.df) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Add PCA if not present
  input.df  = input.df %>% add_reduced_dimensions_PCA()

  # Calculate the new Seurat object
  seurat_obj =
    input.df %>%
    attr("seurat") %>%
    map(~ {
      my_obj =
        (.) %>%
        FindNeighbors() %>%
        FindClusters(method = "igraph")

      # Add renamed annotation
      my_obj@meta.data$cluster = my_obj@meta.data$seurat_clusters

      my_obj
    })

  seurat_obj %>%
    map_dfr(
      ~ .x %>%
        `@` (meta.data) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%
        select(!!cell_column, cluster) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add cluster information of single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#'
#' @return A tt tt object
#'
#' @export
add_cluster_annotation_SNN_sc = function(input.df) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_cluster_annotation_SNN_sc()

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Get cell type information on single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#' @importFrom purrr map2
#' @importFrom purrr map
#' @importFrom purrr map_dfr
#' @importFrom magrittr equals
#'
#' @param input.df A tt object
#'
#' @return A tt tt object
#'
#' @export
get_cell_type_annotation_sc = function(input.df) {
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
    singler = SingleR::CreateSinglerObject(
      as.matrix(
        GetAssayData(
         seurat, slot = "counts", assay = "RNA"
        )
      ),
      annot = NULL,
      seurat@project.name,
      min.genes = min.transcripts,
      technology = "10X",
      species = species,
      citation = "",
      ref.list = list(),
      normalize.gene.length = F,
      variable.genes = "de",
      fine.tune = F,
      do.signatures = F,
      clusters = clusters[f]
    )
    return(singler)
  }

  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  seurat_obj =

    input.df %>%
    attr("seurat") %>%
    map2(input.df %>% attr("parameters"),
         ~ {
           my_clusters = switch(
             "cluster" %in% (.x@meta.data %>% colnames) %>% `!` %>% sum(1),
             .x@meta.data$cluster %>% as.character,
             NULL
           )

           ct_class =
             .x %>%
             run_singleR(my_clusters,
                         species =  .y$species) %$%
             singler %>%
             map(# Get cell-clusters
               ~ {
                 clusters_tibble =
                   .x$SingleR.single$clusters$cl %>%
                   as_tibble(rownames = quo_name(cell_column)) %>%
                   rename(cluster = value)

                 cell_types_tibble =
                   .x$SingleR.clusters$labels %>%

                   # Get rownames only if there are rownames
                   ifelse_pipe(
                     (.) %>% nrow %>% `>` (1),
                     ~ .x %>% as_tibble(rownames = "cluster") ,
                     ~ .x %>% as_tibble()
                   ) %>%
                   mutate(reference = sprintf("Cell type %s", .x$about$RefData))

                 left_join(clusters_tibble, cell_types_tibble)

               }) %>%
             do.call("bind_rows", .) %>%

             # If I don't provide the cluster
             # eliminate it from the table
             # as it will be different for
             # the different data sets of SingleR
             ifelse_pipe(my_clusters %>% is.null,
                         ~ .x %>% select(-cluster)) %>%
             spread(reference, V1)

           #cn = ct_class %>% select(-one_of("cluster")) %>% colnames
           my_obj = .x

           # Add renamed annotation
           my_obj@meta.data =
             my_obj@meta.data %>%
             cbind(
               ct_class[
                 match(
                   my_obj@meta.data %>% rownames ,
                   ct_class %>% pull(!!cell_column)),
                 ] %>%
                 mutate_if(is.character, as.factor) %>%
                 select(-!!cell_column)
              )

           # mutate(cluster = cluster %>% as.character) %>%
           # left_join(ct_class,  by = "cluster") %>%
           # mutate(cluster = seurat_clusters)

           # Return
           my_obj
         })

  seurat_obj %>%
    map_dfr(
      ~ .x %>%
        `@` (meta.data) %>%
        as_tibble(rownames = quo_name(cell_column)) %>%

        # Grub last two columns
        select(!!cell_column, (ncol(.) - 1):ncol(.)) %>%
        mutate_if(is.factor, as.character)
    ) %>%
    mutate_if(is.character, as.factor) %>%

    # Add back the attributes objects
    add_attr(seurat_obj, "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}

#' Add cell type information on single cells
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import Seurat
#'
#' @param input.df A tt object
#'
#' @return A tt tt object
#'
#' @export
add_cell_type_annotation_sc = function(input.df) {
  # Update on tibble
  input.df = input.df %>% update_object_sc

  # Cell column name
  cell_column = input.df %>% attr("parameters") %>% `[[` (1) %$% cell_column

  # Get now object
  input.df.annotated =
    input.df %>%
    get_cell_type_annotation_sc()

  # Merge
  input.df %>%
    left_join(input.df.annotated ,
              by = quo_name(cell_column)) %>%

    # Add back the attributes objects
    add_attr(input.df.annotated %>% attr("seurat"), "seurat") %>%
    add_attr(input.df %>% attr("parameters"), "parameters")

}
