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
      TRUE ~ "FrenchDom"
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
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/hanfrench -f "HanDom,FrenchDom,Han,French" -n hanfrench_merged -o admixpops_test_data/hanfrench_merged')

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

ind_admixpops2_raw <- partitions::compositions(n = 10, m = 3, include.zero = T) |>
  {\(x) x*10}() |>
  as.matrix() |>
  t() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    id = 1:dplyr::n(),
    unit = dplyr::case_when(
      V1 > V2 & V1 > V3 ~ "MbutiDom",
      V2 > V1 & V2 > V3 ~ "HanDom",
      V3 > V1 & V3 > V2 ~ "FrenchDom",
      TRUE ~ "Center"
    )
  ) |> {\(x) { 
    purrr::pmap_chr(
      list(x$id, x$unit, x$V1, x$V2, x$V3),
      \(a,b,c,d,e) {
        paste0("[",a,":",b,"]","(Mbuti=",c,"+Han=",d,"+French=",e,")")
      }
    )
  }}()

ind_admixpops2 <- ind_admixpops2_raw |>
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
  geom_point(aes(x = C1, y = C2, colour = FID)) +
  ggthemes::scale_colour_colorblind()

#### admixture analysis ####

eva.cluster::cluster_down(pw,
  "/mnt/archgen/users/schmid/paagen/playground/admixpops_test_data/admixture_test" ~
    "~/agora/paagen/playground/admixpops_test_data/admixture_test")

hu <- list.files(
  "~/agora/paagen/playground/admixpops_test_data/admixture_test/3",
  recursive = T,
  pattern = ".Q",
  full.names = T
) |> 
  (\(x) Map(\(y) {
    num_chimeras <- length(ind_admixpops2_raw)
    readr::read_delim(y, col_names = F, delim = " ") |> 
      dplyr::mutate(
        run = y,
        ind = c(rep(NA, (dplyr::n() - num_chimeras)), ind_admixpops2_raw)
      )
  }, x))() |>
  dplyr::bind_rows() |>
  dplyr::filter(
    !is.na(ind),
    run == dplyr::first(run)
  ) |>
  tidyr::pivot_longer(
    tidyselect::starts_with("X"),
    names_to = "K",
    values_to = "proportion"
  )

library(ggplot2)
hu |>
  ggplot() +
  geom_bar(
    aes(x = ind, y = proportion, fill = K),
    stat = "identity"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

