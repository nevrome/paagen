source("technical_helpers.R")
library(ggplot2)

# prepare data
#nd("admixpops_test_data")
#nd("admixpops_test_data/plots")

# download
s('trident fetch -d admixpops_test_data -f "*2012_PattersonGenetics-2.1.3*"')

# pruning
nd("admixpops_test_data/plink_patterson_pruned")
s("~/software/plink --bfile admixpops_test_data/2012_PattersonGenetics-2.1.3/2012_PattersonGenetics --exclude pruning_ranges.txt --range --maf --make-bed --out admixpops_test_data/plink_patterson_pruned/pruned")

s("trident init --snpSet Other -p admixpops_test_data/plink_patterson_pruned/pruned.bed -o admixpops_test_data/2012_PattersonGenetics_pruned -n 2012_PattersonGenetics_pruned")

#### one large population test ####

# run admixpops
s('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a "[1:HanDom](Han=100);[2:HanDom](Han=100);[3:HanDom](Han=100);[4:HanDom](Han=100);[5:HanDom](Han=100)" -o admixpops_test_data/han')
# here: tangent to test vcf and zipped output
s('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a "[1c:HanDomChunk](Han=100);[2c:HanDomChunk](Han=100);[3c:HanDomChunk](Han=100);[4c:HanDomChunk](Han=100);[5c:HanDomChunk](Han=100)" -o admixpops_test_data/han_chunks5000 --inChunks --outFormat VCF --zip')

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics_pruned -d admixpops_test_data/han -d admixpops_test_data/han_chunks5000 -f "HanDom,Han,HanDomChunk" -n han_merged -o admixpops_test_data/han_merged')

# mds
nd("admixpops_test_data/han_mds")
s('~/software/plink --bfile admixpops_test_data/han_merged/han_merged --genome --out admixpops_test_data/han_mds/pairwise_stats')
s('~/software/plink --bfile admixpops_test_data/han_merged/han_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/han_mds/pairwise_stats.genome --out admixpops_test_data/han_mds/mds')

# plot
mds_raw <- read_mds("admixpops_test_data/han_mds/mds.mds")

p <- mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

ggsave(
  "admixpops_test_data/plots/han_mds.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

#### one small population test ####

s('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a "[1:BantuSADom](BantuSA=100);[2:BantuSADom](BantuSA=100);[3:BantuSADom](BantuSA=100);[4:BantuSADom](BantuSA=100);[5:BantuSADom](BantuSA=100)" -o admixpops_test_data/BantuSA')
s('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a "[1c:BantuSADomChunk](BantuSA=100);[2c:BantuSADomChunk](BantuSA=100);[3c:BantuSADomChunk](BantuSA=100);[4c:BantuSADomChunk](BantuSA=100);[5c:BantuSADomChunk](BantuSA=100)" -o admixpops_test_data/BantuSA_chunks5000 --inChunks')

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics_pruned -d admixpops_test_data/BantuSA -d admixpops_test_data/BantuSA_chunks5000 -f "BantuSA,BantuSADom,BantuSADomChunk" -n BantuSA_merged -o admixpops_test_data/BantuSA_merged')

# mds
nd("admixpops_test_data/BantuSA_mds")
s('~/software/plink --bfile admixpops_test_data/BantuSA_merged/BantuSA_merged --genome --out admixpops_test_data/BantuSA_mds/pairwise_stats')
s('~/software/plink --bfile admixpops_test_data/BantuSA_merged/BantuSA_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/BantuSA_mds/pairwise_stats.genome --out admixpops_test_data/BantuSA_mds/mds')

# plot
mds_raw <- read_mds("admixpops_test_data/BantuSA_mds/mds.mds")

p <- mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

ggsave(
  "admixpops_test_data/plots/BantuSA_mds.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

#### two populations test ####

combinations <- partitions::compositions(n = 10, m = 2, include.zero = T) |>
  {\(x) x*10}() |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
  dplyr::mutate(
    id = paste(V1, V2, sep = "|"),
    unit = dplyr::case_when(
      V1 > V2 ~ "HanDom",
      V2 > V1 ~ "FrenchDom",
      TRUE ~ "Center"
    )
  )

combinations_chunks <- combinations |> {\(x) {
  x |>
    dplyr::mutate(unit = paste0(unit, "Chunk"))
  }}()

ind_admixpops <- combinations |> {\(x) { 
    purrr::pmap_chr(
      list(x$id, x$unit, x$V1, x$V2),
      \(a,b,c,d) {
        paste0("[",a,":",b,"]","(Han=",c,"+French=",d,")")
      }
    )
  }}() |>
  {\(x) paste(x, collapse = ";")}()

ind_admixpops_chunks <- combinations_chunks |> {\(x) { 
  purrr::pmap_chr(
    list(x$id, x$unit, x$V1, x$V2),
    \(a,b,c,d) {
      paste0("[",a,"c:",b,"]","(Han=",c,"+French=",d,")")
    }
  )
  }}() |>
  {\(x) paste(x, collapse = ";")}()

# run admixpops
s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a \"', ind_admixpops, '\" -o admixpops_test_data/hanfrench'))
s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a \"', ind_admixpops_chunks, '\" -o admixpops_test_data/hanfrench_chunks5000 --inChunks'))

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics_pruned -d admixpops_test_data/hanfrench -d admixpops_test_data/hanfrench_chunks5000 -f "Han,French,HanDom,FrenchDom,Center,HanDomChunk,FrenchDomChunk,CenterChunk" -n hanfrench_merged -o admixpops_test_data/hanfrench_merged')

# mds
nd("admixpops_test_data/mds")
s('~/software/plink --bfile admixpops_test_data/hanfrench_merged/hanfrench_merged --genome --out admixpops_test_data/mds/pairwise_stats')
s('~/software/plink --bfile admixpops_test_data/hanfrench_merged/hanfrench_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/mds/pairwise_stats.genome --out admixpops_test_data/mds/mds')

# plot
mds_raw <- read_mds("admixpops_test_data/mds/mds.mds")

p <- mds_raw |>
  dplyr::mutate(
    label = ifelse(grepl("\\|", IID), IID, NA)
  ) |>
  ggplot(aes(x = C1, y = C2, colour = FID, label = label)) +
  geom_point() +
  ggrepel::geom_label_repel()

ggsave(
  "admixpops_test_data/plots/mds_mds.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

#### three populations test ####

ind_admixpops2_table <- partitions::compositions(n = 10, m = 3, include.zero = T) |>
  {\(x) x*10}() |>
  as.matrix() |>
  t() |>
  as.data.frame() |>
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

ind_admixpops2_table_chunks <- ind_admixpops2_table |> {\(x) {
  x |>
    dplyr::mutate(unit = paste0(unit, "Chunk"))
  }}()

ind_admixpops2 <- ind_admixpops2_table |> {\(x) { 
  purrr::pmap_chr(
      list(x$id, x$unit, x$Mbuti, x$Han, x$French),
      \(a,b,c,d,e) {
        paste0("[",a,":",b,"]","(Mbuti=",c,"+Han=",d,"+French=",e,")")
      }
    )
  }}() |>
  {\(x) paste(x, collapse = ";")}()

ind_admixpops2_chunks <- ind_admixpops2_table_chunks |> {\(x) { 
  purrr::pmap_chr(
    list(x$id, x$unit, x$Mbuti, x$Han, x$French),
    \(a,b,c,d,e) {
      paste0("[",a,"c:",b,"]","(Mbuti=",c,"+Han=",d,"+French=",e,")")
    }
  )
}}() |>
  {\(x) paste(x, collapse = ";")}()

# run admixpops
s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a \"', ind_admixpops2, '\" -o admixpops_test_data/mbutihanfrench'))
s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a \"', ind_admixpops2_chunks, '\" -o admixpops_test_data/mbutihanfrench_chunks5000 --inChunks'))

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics_pruned -d admixpops_test_data/mbutihanfrench -d admixpops_test_data/mbutihanfrench_chunks5000 -f "Mbuti,Han,French,MbutiDom,HanDom,FrenchDom,Center,MbutiDomChunk,HanDomChunk,FrenchDomChunk,CenterChunk" -n mbutihanfrench_merged -o admixpops_test_data/mbutihanfrench_merged')

# mds
nd("admixpops_test_data/mbutihanfrench_mds")
s('~/software/plink --bfile admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged --genome --out admixpops_test_data/mbutihanfrench_mds/pairwise_stats')
s('~/software/plink --bfile admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/mbutihanfrench_mds/pairwise_stats.genome --out admixpops_test_data/mbutihanfrench_mds/mds')

# plot
mds_raw <- read_mds("admixpops_test_data/mbutihanfrench_mds/mds.mds")

p <- mds_raw |>
  dplyr::mutate(
    method = dplyr::case_when(
      grepl("Chunk", FID) ~ "in Chunks",
      (!grepl("Dom", FID) & FID != "Center") ~ "input",
      TRUE ~ "per SNP"
    )
  ) |>
  #dplyr::filter(method %in% c("per SNP", "input")) |>
  dplyr::filter(method %in% c("in Chunks", "input")) |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

ggsave(
  "admixpops_test_data/plots/mbutihanfrench_mds_chunks.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

#### three pops per SNP: now with --marginalizeMissing ####

# run admixpops
s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a \"', ind_admixpops2, '\" --marginalizeMissing -o admixpops_test_data/mbutihanfrench_mm'))

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics_pruned -d admixpops_test_data/mbutihanfrench_mm -f "MbutiDom,HanDom,FrenchDom,Center,Mbuti,Han,French" -n mbutihanfrench_mm_merged -o admixpops_test_data/mbutihanfrench_mm_merged')

# mds
nd("admixpops_test_data/mbutihanfrench_mm_mds")
s('~/software/plink --bfile admixpops_test_data/mbutihanfrench_mm_merged/mbutihanfrench_mm_merged --genome --out admixpops_test_data/mbutihanfrench_mm_mds/pairwise_stats')
s('~/software/plink --bfile admixpops_test_data/mbutihanfrench_mm_merged/mbutihanfrench_mm_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/mbutihanfrench_mm_mds/pairwise_stats.genome --out admixpops_test_data/mbutihanfrench_mm_mds/mds')

# plot
mds_raw <- read_mds("admixpops_test_data/mbutihanfrench_mm_mds/mds.mds")

p <- mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

ggsave(
  "admixpops_test_data/plots/mbutihanfrench_mm_mds.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

#### three pops: admixture analysis ####

# create .pop file for supervised admixture
fam <- readr::read_tsv(
  "admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged.fam",
  col_names = F
)

writeLines(
  ifelse(
    fam$X1 %in% c("French", "Mbuti", "Han"), fam$X1, "-"
  ),
  "admixpops_test_data/mbutihanfrench_merged/mbutihanfrench_merged.pop"
)

pw <- "..."
u <- "..."
h <- "daghead1.eva.mpg.de"

# upload data to cluster
eva.cluster::cluster_up(
  "~/agora/paagen/admixpops_test_data/mbutihanfrench_merged/" ~
    "/mnt/archgen/users/schmid/paagen/admixpops_test_data/mbutihanfrench_merged/",
  user = u, host = h, pw = pw
)

# run on cluster: qsub runAdmixture.sh

eva.cluster::cluster_down(
  "/mnt/archgen/users/schmid/paagen/admixpops_test_data/admixture_test/" ~
    "~/agora/paagen/admixpops_test_data/admixture_test/",
  user = u, host = h, pw = pw
)

# read all results and bring them together
merged_admixture_results_wide <- list.files(
  "~/agora/paagen/admixpops_test_data/admixture_test/3",
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

# transform data to ggtern
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

# plot
library(ggtern)
p <- ggtern() +
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

ggsave(
  "admixpops_test_data/plots/mbutihanfrench_admix.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

#### independent runs, but with other pops in MDS ####

s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a "[6:FrenchDom](French=100);[7:FrenchDom](French=100);[8:FrenchDom](French=100);[9:FrenchDom](French=100);[10:FrenchDom](French=100)" -o admixpops_test_data/french'))

s(paste0('xerxes admixpops -d admixpops_test_data/2012_PattersonGenetics_pruned -a "[11:MbutiDom](Mbuti=100);[12:MbutiDom](Mbuti=100);[13:MbutiDom](Mbuti=100);[14:MbutiDom](Mbuti=100);[15:MbutiDom](Mbuti=100)" -o admixpops_test_data/mbuti'))

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics_pruned -p admixpops_test_data/han/han.bed -p admixpops_test_data/french/french.bed -p admixpops_test_data/mbuti/mbuti.bed -f "MbutiDom,FrenchDom,HanDom,Han,French,Mbuti" -n independenthanfrenchmbuti_merged -o admixpops_test_data/independenthanfrenchmbuti_merged')

# mds
nd("admixpops_test_data/independenthanfrenchmbuti_mds")
s('~/software/plink --bfile admixpops_test_data/independenthanfrenchmbuti_merged/independenthanfrenchmbuti_merged --genome --out admixpops_test_data/independenthanfrenchmbuti_mds/pairwise_stats')
s('~/software/plink --bfile admixpops_test_data/independenthanfrenchmbuti_merged/independenthanfrenchmbuti_merged --cluster --mds-plot 2 --read-genome admixpops_test_data/independenthanfrenchmbuti_mds/pairwise_stats.genome --out admixpops_test_data/independenthanfrenchmbuti_mds/mds')

# plot
mds_raw <- read_mds("admixpops_test_data/independenthanfrenchmbuti_mds/mds.mds")

p <- mds_raw |>
  ggplot() +
  geom_point(aes(x = C1, y = C2, colour = FID))

ggsave(
  "admixpops_test_data/plots/independenthanfrenchmbuti_mds.jpeg",
  plot = p,
  device = "jpeg",
  width = 10,
  height = 6,
  scale = 0.8
)

