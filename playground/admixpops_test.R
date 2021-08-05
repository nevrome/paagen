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

# prepare data
nd("admixpops_test_data")
s('trident fetch -d admixpops_test_data -f "*2012_PattersonGenetics*"')

# run admixpops
s('paagen admixpops -d admixpops_test_data/2012_PattersonGenetics -a "[A1:A](French=60+Han=40);[A2:A](French=70+Han=30);[A3:A](French=80+Han=20);[B1:B](French=40+Han=60);[B2:B](French=30+Han=70);[B3:B](French=20+Han=80)" -o admixpops_test_data/hanfrench')

# create data subset
s('trident init --inFormat EIGENSTRAT --snpSet HumanOrigins --genoFile admixpops_test_data/hanfrench/res.geno --indFile admixpops_test_data/hanfrench/res.ind --snpFile admixpops_test_data/hanfrench/res.snp -o admixpops_test_data/hanfrench_poseidon -n hanfrench_poseidon')
s('trident forge -d admixpops_test_data/2012_PattersonGenetics -d admixpops_test_data/hanfrench_poseidon -f "A,B,French,Han" -n hanfrench_merged -o admixpops_test_data/hanfrench_merged')

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

