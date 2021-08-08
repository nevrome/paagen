# helper functions
s <- function(x, o = T, e = F) {
    redir <- if (e) { "2>&1" } else { "" }
    res <- system(paste(x, redir), intern=TRUE, ignore.stdout = !o)
    if (length(res) > 1) { cat(res, sep='\n') }
}
dd <- function(x) { unlink(x, recursive = T) }
nd <- function(x) {
    unlink(x, recursive = T)
    dir.create(x)
}

library(magrittr)

nd("spacetime_test_data")
s('trident fetch -d spacetime_test_data -f "*2012_PattersonGenetics*"')

janno_raw <- poseidonR::read_janno("spacetime_test_data/2012_PattersonGenetics")

# filtering to sampels with spati
janno_filtered <- janno_raw %>% dplyr::filter(
    !is.na(Latitude) & !is.na(Longitude)
)

# Nr_autosomal_SNPs: should be >= 20000 SNPs
janno_QC <- janno_filtered %>% dplyr::filter(
  Nr_autosomal_SNPs >= 20000
)
# Xcontam: if male, then should not be higher then 10%
janno_QC <- janno_QC %>% dplyr::filter(
  is.na(Xcontam) | Genetic_Sex == "F" | (Genetic_Sex == "M" & Xcontam < 0.1)
)
# Genetic_Sex: Individuals with unknown genetic sex should be removed
janno_QC <- janno_QC %>% dplyr::filter(Genetic_Sex != "U")
# Indicated as contaminated: Individuals which are indicated as potentially contaminated
# in their ID should be removed
janno_QC <- janno_QC %>% dplyr::filter(
  !grepl("cont|excluded|Ignore", x = Individual_ID, ignore.case = T) &
    !grepl("cont|excluded|Ignore", x = Group_Name, ignore.case = T)
)

janno_final <- janno_QC

save(janno_final, file = "spacetime_test_data/janno_final.RData")

load("spacetime_test_data/janno_final.RData")

# store ind list for poseidon extraction
tibble::tibble(
  #pop = sapply(janno_filtered_final$Group_Name, function(x) { x[[1]] }),
  ind = paste0("<", sort(janno_final$Individual_ID), ">")
) %>% 
  readr::write_delim(
    file = "spacetime_test_data/ind_list.txt",
    delim = " ",
    col_names = FALSE
  )

dd("spacetime_test_data/pat")
s('trident forge --forgeFile spacetime_test_data/ind_list.txt -d spacetime_test_data/2012_PattersonGenetics -n pat -o spacetime_test_data/pat')

#manual_pois <- tibble::tribble(
#    ~time, ~lat, ~lon,
#    2000, 46, 3,
#    2000, 46, 3,
#    2000, 46, 3,
#    2000, 55, 37    
#)

world <- spData::world 

poi_grid <- world %>% 
    dplyr::filter(continent != "Antarctica") %>%
    sf::st_make_grid(cellsize = 20, what = "centers") %>%
    sf::st_sf() %>%
    sf::st_intersection(world) %>%
    dplyr::mutate(
      lon = sf::st_coordinates(.)[,1],
      lat = sf::st_coordinates(.)[,2]
    ) %>%
    sf::st_drop_geometry() %>%
    dplyr::transmute(
      ind = paste0("poi", 1:(dplyr::n())),
      group = gsub(" ", "_", continent),
      time = 2000,
      lon = round(lon, 3),
      lat = round(lat, 3)
    )

poi_string <- purrr::pmap(poi_grid, function(ind, group, time, lat, lon) {
    paste0("[", ind, ":", group, "](", paste(time, lat, lon, sep = ","), ")")
}) %>% paste(collapse = ";")

poi_string

library(ggplot2)

janno_grouped <- janno_final %>%
  dplyr::group_by(Group_Name, Latitude, Longitude) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

ggplot() +
  geom_sf(data = world) +
  geom_point(data = janno_grouped, aes(x = Longitude, y = Latitude, size = n)) +
  ggrepel::geom_label_repel(
    data = janno_grouped, 
    aes(x = Longitude, y = Latitude, label = Group_Name), 
    color = "grey", size = 3, max.overlaps = 100
  ) +
  geom_text(data = poi_grid, aes(x = lon, y = lat, color = group, label = ind))

nd("spacetime_test_data/poi")
s(paste0('paagen spacetime -d spacetime_test_data/pat -p "', poi_string, '" --neighbors 100 -o spacetime_test_data/poi --outFormat EIGENSTRAT'))

dd("spacetime_test_data/poi_poseidon")
s('trident init --inFormat EIGENSTRAT --genoFile spacetime_test_data/poi/poi.geno --snpFile spacetime_test_data/poi/poi.snp --indFile spacetime_test_data/poi/poi.ind -o spacetime_test_data/poi_poseidon -n poi')

dd("spacetime_test_data/merged")
s('trident forge -d spacetime_test_data/pat -d spacetime_test_data/poi_poseidon -f "*pat*,*poi*" -o spacetime_test_data/merged -n merged')

# pruning
nd("spacetime_test_data/merged_pruned")
s('plink1.9 --bfile spacetime_test_data/merged/merged --exclude spacetime_test_data/myrange.txt --range --maf --make-bed --out spacetime_test_data/merged_pruned/merged.pruned')

# generate general pairwise stats
nd("spacetime_test_data/merge_pruned_distances")
s('plink1.9 --bfile spacetime_test_data/merged_pruned/merged.pruned --genome --out spacetime_test_data/merge_pruned_distances/merged.pruned')

# create mds table
nd("spacetime_test_data/mds")
s('plink1.9 --bfile spacetime_test_data/merged_pruned/merged.pruned --cluster --mds-plot 2 --read-genome spacetime_test_data/merge_pruned_distances/merged.pruned.genome --out spacetime_test_data/mds/mds')

mds_raw <- readr::read_delim(
        "spacetime_test_data/mds/mds.mds", " ", trim_ws = T,
        col_types = "ccddd_"
    )

load("spacetime_test_data/janno_final.RData")

input_spatpos <- janno_final %>% dplyr::transmute(
    ind = Individual_ID,
    group = sapply(Group_Name, function(x){x[1]}),
    time = 2000,
    lon = Longitude,
    lat = Latitude
)  %>% 
dplyr::left_join(
    mds_raw, by = c("ind" = "IID")
)

input_spatpos_grouped <- input_spatpos %>%
dplyr::group_by(FID) %>%
dplyr::summarise(
    C1 = mean(C1),
    C2 = mean(C2)
)

input_grid <- poi_grid %>% 
dplyr::left_join(
    mds_raw, by = c("ind" = "IID")
)

head(input_spatpos_grouped)

head(input_grid)

library(ggplot2)
options(repr.plot.width = 20, repr.plot.height = 7, repr.plot.res = 300)

ggplot() +
ggpointgrid::geom_textgrid(data = input_grid, aes(x = C1, y = C2, color = group, label = ind), size = 7) +
ggpointgrid::geom_textgrid(data = input_spatpos_grouped, aes(x = C1, y = C2, label = FID))
