context('Single cell functions')

# Build tt object
tt =
  create_tt_from_tibble_sc(
    tidysc:: counts,
    .sample = sample,
    .cell = cell,
    .transcript = transcript,
    .abundance = `count`,
    species = "Human"
  )

test_that("Test data frame",{ expect_equal( ncol(tidysc::counts), 6 ) })

test_that("Create tt object from tibble",{

  expect_equal( ncol(tt), 11 )

  expect_equal( typeof(attr(tt, "parameters")), "list")

})

# test_that("Create tt object from Cell ranger",{
#
#   my_tt =
#     create_tt_from_cellRanger_sc(
#       grep(
#         "-tc",
#         list.dirs("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/triple-therapy-new-run/cell_ranger_input", recursive=F, full.names = T),
#         invert = T,
#         value = T
#       )[1:2],
#       species = "Mouse"
#     )
#
#   expect_equal( ncol(my_tt), 11 )
#
#   expect_equal( typeof(attr(my_tt, "parameters")), "list")
#
# })

test_that("Get mitochondrial abundance",{

  my_tt = get_mitochndrion_transcription_abundance_sc(tt, cell)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$cell[1:4], c("D101_43_1", "D101_5_1" , "D101_50_1", "D101_51_1"))

})

test_that("Add mitochondrial abundance",{

  my_tt =  add_mitochndrion_transcription_abundance_sc(tt, cell)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$cell[1:4], c("D101_43_1", "D101_5_1" , "D101_50_1", "D101_51_1"))

})

test_that("Get normalised counts",{

  my_tt = get_normalised_counts_sc(tt)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$cell[1:4], c("D101_43_1", "D101_5_1" , "D101_50_1", "D101_51_1"))

})

test_that("Add normalised counts",{

  my_tt = add_normalised_counts_sc(tt)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$cell[1:4], c("D101_43_1", "D101_5_1" , "D101_50_1", "D101_51_1"))

})

test_that("Add variable genes classification",{

  my_tt =  add_variable_genes_classification(tt)

  expect_equal( ncol(my_tt), 11 )
  expect_equal( my_tt$`count total`[1:4], c(8839,  3613, 11841, 14665))

})

test_that("Get cell cycle annotation",{

  my_tt =  get_cell_cycle_annotation_sc(tt, cell)

  expect_equal( ncol(my_tt), 4 )
  expect_equal( my_tt$S.Score[1:4], c(-0.2367882, -0.0902214, -0.2489621, -0.5640438), tolerance=1e-7)

})

test_that("Add cell cycle annotation",{

  my_tt =  add_cell_cycle_annotation_sc(tt, cell)

  expect_equal( ncol(my_tt), 14 )
  expect_equal( my_tt$S.Score.x[1:4], c(-0.4892333, -0.1579930, -0.3676411, -0.8946705), tolerance=1e-7)

})

test_that("Get reduced dimensions PCA",{

  .dims = 10

  my_tt =  get_reduced_dimensions_PCA(tt, .dims = 10)

  expect_equal( ncol(my_tt), .dims + 1 )
  expect_equal( my_tt$`PC 1`[1:4], c(6.279411 ,  4.115689, -27.747529, -32.714527), tolerance=1e-7)

})

test_that("Add reduced dimensions PCA",{

  .dims = 10

  my_tt =  add_reduced_dimensions_PCA(tt, .dims = 10)

  expect_equal( ncol(my_tt), .dims + 11 )
  expect_equal( my_tt$`PC 1`[1:4], c(21.75515, 14.47006, 17.83543, 16.21051), tolerance=1e-7)

})

test_that("Get reduced dimensions UMAP",{


  my_tt =  get_reduced_dimensions_UMAP(tt)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$`UMAP 1`[1:4], c(1.306986, 1.476799, 1.620057, 1.563925), tolerance=1e-7)

})

test_that("Add reduced dimensions UMAP",{


  my_tt =  add_reduced_dimensions_UMAP(tt)

  expect_equal( ncol(my_tt), .dims + 3 )
  expect_equal( my_tt$`UMAP 1`[1:4], c(1.306986, 1.476799, 1.620057, 1.563925), tolerance=1e-7)

})

test_that("Get reduced dimensions TSNE",{


  my_tt =  get_reduced_dimensions_TSNE(tt)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$`tSNE 1`[1:4], c(7.509182, 8.228692, 7.877330, 7.921603), tolerance=1e-7)

})

test_that("Add reduced dimensions TSNE",{


  my_tt =  add_reduced_dimensions_TSNE(tt)

  expect_equal( ncol(my_tt), 13 )
  expect_equal( my_tt$`tSNE 1`[1:4], c(7.509182, 8.228692, 7.877330, 7.921603), tolerance=1e-7)

})

test_that("Get adjusted counts for unwanted variation",{


  my_tt =   get_adjusted_counts_for_unwanted_variation_sc(tt, ~ integrate(sample) + S.Score + G2M.Score + mito.fraction)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$nCount_SCT[1:4], c(11936, 11471, 12267, 12848))

})

test_that("Add adjusted counts for unwanted variation",{


  my_tt =   add_adjusted_counts_for_unwanted_variation_sc(tt, ~ integrate(sample)  + S.Score + G2M.Score + mito.fraction)

  expect_equal( ncol(my_tt), 13 )
  expect_equal( my_tt$nCount_SCT[1:4], c(11936, 11471, 12267, 12848))

})

test_that("Get cluster annotation SNN",{


  my_tt =   get_cluster_annotation_SNN_sc(tt)

  expect_equal( ncol(my_tt), 2 )
  expect_equal( my_tt$cluster[1:4], c(2, 2, 1, 1))

})

test_that("Add cluster annotation SNN",{


  my_tt =   add_cluster_annotation_SNN_sc(tt)

  expect_equal( ncol(my_tt), 12 )
  expect_equal( my_tt$cluster[1:4], c(2, 2, 1, 1))

})

test_that("Get cell type annotation",{


  my_tt =   get_cell_type_annotation_sc(tt)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$`Cell type Blueprint_Encode`[1:4], c("Neurons", "Neurons", "Epithelial cells", "Epithelial cells"))

})

test_that("Add cell type annotation",{


  my_tt =   add_cell_type_annotation_sc(tt)

  expect_equal( ncol(my_tt), 3 )
  expect_equal( my_tt$`Cell type Blueprint_Encode`[1:4], c("Neurons", "Neurons", "Epithelial cells", "Epithelial cells"))

})
