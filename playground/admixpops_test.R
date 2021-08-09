source("technical_helpers.R")

# prepare data
nd("admixpops_test_data")
s('trident fetch -d admixpops_test_data -f "*2012_PattersonGenetics*"')

ind_admixpops <- partitions::compositions(n = 5, m = 2, include.zero = T) |>
  {\(x) x*20}() |>
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

####

ind_admixpops2 <- partitions::compositions(n = 5, m = 3, include.zero = T) |>
  {\(x) x*20}() |>
  as.matrix() |>
  t() |>
  tibble::as_tibble() |>
  dplyr::mutate(
    id = 1:dplyr::n(),
    unit = dplyr::case_when(
      V1 > V2 & V1 > V3 ~ "MbutiDom",
      V2 > V1 & V2 > V3 ~ "HanDom",
      V3 > V1 & V3 > V2 ~ "FrenchDom",
      TRUE ~ "D"
    )
  ) |> {\(x) { 
    purrr::pmap_chr(
      list(x$id, x$unit, x$V1, x$V2, x$V3),
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
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/mbutihanfrench -f "MbutiDom,HanDom,FrenchDom,Mbuti,Han,French" -n mbutihanfrench_merged -o admixpops_test_data/mbutihanfrench_merged')

# mds
nd("admixpops_test_data/mbuti_mds")
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
