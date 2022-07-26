# read mds
read_mds <- function(x) {
  readr::read_fwf(
    file = x, 
    col_positions = readr::fwf_empty(
      x,
      skip = 1,
      col_names = c("FID", "IID", "SOL", "C1", "C2"),
      n = 1000
    ),
    trim_ws = T,
    col_types = "ccddd_",
    skip = 1
  )
}

# run system command RStudio
s <- function(x) {
  termId <- rstudioapi::terminalCreate()
  rstudioapi::terminalSend(
    termId,
    paste0(
      x,
      "\n"
    )
  )
}

# run system command
sb <- function(x, o = T, e = T) {
  redir <- if (e) { "2>&1" } else { "" }
  res <- system(paste(x, redir), intern=TRUE, ignore.stdout = !o)
  if (length(res) > 1) { cat(res, sep='\n') }
}

# delete directory
dd <- function(x) { unlink(x, recursive = T) }

# delete directory and then create it again
nd <- function(x) {
  unlink(x, recursive = T)
  dir.create(x)
}
