source("technical_helpers.R")

# prepare data
nd("admixpops_test_data")
s('trident fetch -d admixpops_test_data -f "*2012_PattersonGenetics*"')

# run admixpops
s('paagen admixpops -d admixpops_test_data/2012_PattersonGenetics -a "[1:A](French=90+Han=10);[2:A](French=80+Han=20);[3:A](French=70+Han=30);[4:A](French=60+Han=40);[5:C](French=50+Han=50);[6:B](French=40+Han=60);[7:B](French=30+Han=70);[8:B](French=20+Han=80);[9:B](French=10+Han=90)" -o admixpops_test_data/hanfrench')

# create data subset
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/hanfrench -f "A,B,C,French,Han" -n hanfrench_merged -o admixpops_test_data/hanfrench_merged')

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

