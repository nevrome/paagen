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
