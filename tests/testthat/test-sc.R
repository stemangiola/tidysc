context('Single cell functions')

`%>%` = magrittr::`%>%`
  
# Build tt object
tt =
  tidysc_long(
    tidysc::counts,
    .sample = sample,
    .cell = cell,
    .transcript = transcript,
    .abundance = count,
    species = "Human"
  ) %>%
  filter(!low_quality)

test_that("Test data frame",{ expect_equal( ncol(tidysc::counts), 6 ) })

test_that("Create tt object from tibble",{

  expect_equal( ncol(tt), 12 )

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

tt_scaled = scale_abundance(tt)

test_that("Get scaled counts",{

  expect_equal( ncol(tt_scaled), 14 )
  expect_equal( as.character(tt_scaled$cell[1:4]), c("D101_43_1", "D101_5_1" , "D101_50_1", "D101_51_1"))

})

test_that("Get reduced dimensions PCA",{

  my_tt =  tt_scaled %>% reduce_dimensions(.dims = 10, method = "PCA")

  expect_equal( ncol(my_tt), 24 )
  expect_equal( my_tt$`PC 1`[1:4], c(-6.301245, -4.106271, 27.878208, 32.859274), tolerance=1e-7)

})

test_that("Get reduced dimensions UMAP",{


  my_tt =  tt_scaled %>% reduce_dimensions( method = "UMAP")
  
  expect_equal( ncol(my_tt), 16 )
  expect_equal( my_tt$`UMAP 1`[1:4], c( 11.811723, 11.702131, -7.572128, -6.526407), tolerance=1e-7)

})


test_that("Get reduced dimensions TSNE",{


  my_tt =  tt_scaled %>% reduce_dimensions( method = "tSNE")
  
  expect_equal( ncol(my_tt), 16 )
  expect_equal( my_tt$`tSNE 1`[1:4], c( 9.304803 ,  9.459389, -22.754654, -16.733519), tolerance=1e-7)

})

test_that("Get adjusted counts for unwanted variation",{


  my_tt =   tt %>% adjust_abundance(~ integrate(sample) + S.Score + G2M.Score + mito.fraction)

  expect_equal( ncol(my_tt), 14 )
  expect_equal( my_tt$nCount_SCT[1:4], c(11630 ,11112,12046 ,12611))

})


test_that("Get cluster annotation SNN",{


  my_tt =   cluster_elements(tt_scaled)

  expect_equal( ncol(my_tt), 15 )
  expect_equal( as.integer(as.character(my_tt$cluster[1:4])), c(2, 2, 1, 1))

})



test_that("Get cell type annotation",{


  my_tt =   deconvolve_cellularity(tt_scaled)

  expect_equal( ncol(my_tt), 18 )
  expect_equal( my_tt$label_blueprint[1:4], c("B-cells"  ,       "Neurons"   ,      "Keratinocytes" ,  "Mesangial cells"))

})

