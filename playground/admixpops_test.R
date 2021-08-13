# setwd("/home/schmid/agora/paagen/playground") 
source("technical_helpers.R")

# prepare data
nd("admixpops_test_data")
s('trident fetch -d admixpops_test_data -f "*2012_PattersonGenetics*"')

#### one large population test ####

# run admixpops
s(paste0('paagen admixpops -d admixpops_test_data/2012_PattersonGenetics -a "[1:HanDom](Han=100);[2:HanDom](Han=100);[3:HanDom](Han=100);[4:HanDom](Han=100);[5:HanDom](Han=100)" -o admixpops_test_data/han'))

# create data subset
dd("admixpops_test_data/han_merged")
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/han -f "HanDom,Han" -n han_merged -o admixpops_test_data/han_merged')

# mds
nd("admixpops_test_data/han_mds")
s('plink1.9 --bfile admixpops_test_data/han_merged/han_merged --genome --out admixpops_test_data/han_mds/pairwise_stats')
s('plink1.9 --bfile admixpops_test_data/han_merged/han_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/han_mds/pairwise_stats.genome --out admixpops_test_data/han_mds/mds')

# plot
mds_raw <- readr::read_delim(
  "admixpops_test_data/han_mds/mds.mds", " ", trim_ws = T,
  col_types = "ccddd_"
)

library(ggplot2)
mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

#### one small population test ####

s(paste0('paagen admixpops -d admixpops_test_data/2012_PattersonGenetics -a "[1:BantuSADom](BantuSA=100);[2:BantuSADom](BantuSA=100);[3:BantuSADom](BantuSA=100);[4:BantuSADom](BantuSA=100);[5:BantuSADom](BantuSA=100)" -o admixpops_test_data/BantuSA'))

# create data subset
dd("admixpops_test_data/BantuSA_merged")
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/BantuSA -f "BantuSADom,BantuSA" -n BantuSA_merged -o admixpops_test_data/BantuSA_merged')

# mds
nd("admixpops_test_data/BantuSA_mds")
s('plink1.9 --bfile admixpops_test_data/BantuSA_merged/BantuSA_merged --genome --out admixpops_test_data/BantuSA_mds/pairwise_stats')
s('plink1.9 --bfile admixpops_test_data/BantuSA_merged/BantuSA_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/BantuSA_mds/pairwise_stats.genome --out admixpops_test_data/BantuSA_mds/mds')

# plot
mds_raw <- readr::read_delim(
  "admixpops_test_data/BantuSA_mds/mds.mds", " ", trim_ws = T,
  col_types = "ccddd_"
)

library(ggplot2)
mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

#### two populations test ####

ind_admixpops <- partitions::compositions(n = 10, m = 2, include.zero = T) |>
  {\(x) x*10}() |>
  as.matrix() |>
  t() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    id = 1:dplyr::n(),
    unit = dplyr::case_when(
      V1 > V2 ~ "HanDom",
      V2 > V1 ~ "FrenchDom",
      TRUE ~ "Center"
    )
  ) |> {\(x) { 
    purrr::pmap_chr(
      list(x$id, x$unit, x$V1, x$V2),
      \(a,b,c,d) {
        paste0("[",a,":",b,"]","(Han=",c,"+French=",d,")")
      }
    )
  }}() |>
  {\(x) paste(x, collapse = ";")}()

# run admixpops
s(paste0('paagen admixpops -d admixpops_test_data/2012_PattersonGenetics -a \"', ind_admixpops, '\" -o admixpops_test_data/hanfrench'))

# create data subset
dd("admixpops_test_data/hanfrench_merged")
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/hanfrench -f "HanDom,FrenchDom,Center,Han,French" -n hanfrench_merged -o admixpops_test_data/hanfrench_merged')

# mds
nd("admixpops_test_data/mds")
s('plink1.9 --bfile admixpops_test_data/hanfrench_merged/hanfrench_merged --genome --out admixpops_test_data/mds/pairwise_stats')
s('plink1.9 --bfile admixpops_test_data/hanfrench_merged/hanfrench_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/mds/pairwise_stats.genome --out admixpops_test_data/mds/mds')

# plot
mds_raw <- readr::read_delim(
  "admixpops_test_data/mds/mds.mds", " ", trim_ws = T,
  col_types = "ccddd_"
)

library(ggplot2)
mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

#### three populations test ####

ind_admixpops2_table <- partitions::compositions(n = 10, m = 3, include.zero = T) |>
  {\(x) x*10}() |>
  as.matrix() |>
  t() |>
  tibble::as_tibble() |>
  stats::setNames(c("Mbuti", "Han", "French")) |>
  dplyr::mutate(
    id = 1:dplyr::n(),
    unit = dplyr::case_when(
      Mbuti > Han & Mbuti > French ~ "MbutiDom",
      Han > Mbuti & Han > French ~ "HanDom",
      French > Mbuti & French > Han ~ "FrenchDom",
      TRUE ~ "Center"
    )
  )

ind_admixpops2 <- ind_admixpops2_table |> {\(x) { 
  purrr::pmap_chr(
      list(x$id, x$unit, x$Mbuti, x$Han, x$French),
      \(a,b,c,d,e) {
        paste0("[",a,":",b,"]","(Mbuti=",c,"+Han=",d,"+French=",e,")")
      }
    )
  }}() |>
  {\(x) paste(x, collapse = ";")}()

# run admixpops
s(paste0('paagen admixpops -d admixpops_test_data/2012_PattersonGenetics -a \"', ind_admixpops2, '\" -o admixpops_test_data/mbutihanfrench'))

# create data subset
dd("admixpops_test_data/mbutihanfrench_merged")
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/mbutihanfrench -f "MbutiDom,HanDom,FrenchDom,Center,Mbuti,Han,French" -n mbutihanfrench_merged -o admixpops_test_data/mbutihanfrench_merged')

# mds
nd("admixpops_test_data/mbutihanfrench_mds")
s('plink1.9 --bfile admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged --genome --out admixpops_test_data/mbuti_mds/pairwise_stats')
s('plink1.9 --bfile admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/mbuti_mds/pairwise_stats.genome --out admixpops_test_data/mbuti_mds/mds')

# plot
mds_raw <- readr::read_delim(
  "admixpops_test_data/mbuti_mds/mds.mds", " ", trim_ws = T,
  col_types = "ccddd_"
)

library(ggplot2)
mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

#### admixture analysis ####

eva.cluster::cluster_down(pw,
  "/mnt/archgen/users/schmid/paagen/playground/admixpops_test_data/admixture_test" ~
    "~/agora/paagen/playground/admixpops_test_data/admixture_test")

merged_admixture_results_wide <- list.files(
  "~/agora/paagen/playground/admixpops_test_data/admixture_test/3",
  recursive = T,
  pattern = ".Q",
  full.names = T
) |> 
  (\(x) Map(\(y) {
    num_chimeras <- length(ind_admixpops2_table)
    raw_out <- readr::read_delim(y, col_names = F, delim = " ")
    sorted_dims <- raw_out[raw_out |> colSums() |> sort() |> names()] |> stats::setNames(c("OMbuti", "OFrench", "OHan"))
    sorted_dims |>
      (\(x) x[(nrow(x) - nrow(ind_admixpops2_table) + 1):nrow(x),])() |>
      dplyr::mutate(run = y) |>
      dplyr::bind_cols(ind_admixpops2_table)
  }, x))() |>
  dplyr::bind_rows() |>
  dplyr::group_by(id) |>
  dplyr::summarise(
    mean_OMbuti = mean(OMbuti),
    mean_OHan = mean(OHan),
    mean_OFrench = mean(OFrench),
    # sd is trivially small!
    IMbuti = dplyr::first(Mbuti),
    IHan = dplyr::first(Han),
    IFrench = dplyr::first(French),
    unit = dplyr::first(unit)
  )

merged_admixture_results_long <- merged_admixture_results_wide |>
  tidyr::pivot_longer(
    tidyselect::starts_with(c("mean_", "I"), ignore.case = F)
  ) |>
  dplyr::mutate(
    type = dplyr::case_when(
      grepl("mean", name) ~ "out_paagen+admixture",
      TRUE ~ "in_theoretical"
    ),
    name = dplyr::case_when(
      grepl("Mbuti", name) ~ "Mbuti",
      grepl("Han", name) ~ "Han",
      grepl("French", name) ~ "French"
    )
  ) |>
  tidyr::pivot_wider(
    id_cols = c(id, unit, type),
    names_from = name,
    values_from = value
  )

library(ggtern)
ggtern() +
  geom_segment(
    data = merged_admixture_results_wide,
    aes(mean_OMbuti, mean_OHan, mean_OFrench, xend = IMbuti, yend = IHan, zend = IFrench), alpha = 0.5
  ) +
  geom_point(
    data = merged_admixture_results_long,
    aes(Mbuti, Han, French, color = unit, shape = type),
    size = 2
  ) +
  theme_nomask() +
  xlab("Mbuti") +
  ylab("Han") +
  zlab("French")
