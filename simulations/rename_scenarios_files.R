# Cartella contenente i file
path <- "percorso/della/cartella"

# 20 repliche da conservare
# sel <- sample(1:40, 20)

stopifnot(
  length(sel) == 20,
  length(unique(sel)) == 20,
  all(sel %in% 1:40)
)

# Pattern valido per Y, Delta e Lambda
file_pattern <- "_(Y|Delta|Lambda)_m[0-9]+_H[0-9]+_T[0-9]+_rep[0-9]+\\.csv$"

# Elenco di tutti i file interessati
files <- list.files(
  path       = path,
  pattern    = file_pattern,
  full.names = TRUE
)

# Estrazione del numero rr
rr_file <- as.integer(
  sub(".*_rep([0-9]+)\\.csv$", "\\1", basename(files))
)

# ------------------------------------------------------------------
# 1. Elimina tutti i file Y, Delta e Lambda con rr non incluso in sel
# ------------------------------------------------------------------

files_to_delete <- files[!rr_file %in% sel]

if (length(files_to_delete) > 0) {
  deleted <- file.remove(files_to_delete)

  if (!all(deleted)) {
    warning("Non è stato possibile eliminare alcuni file.")
  }
}

# ------------------------------------------------------------------
# 2. Aggiunge "_temp" a tutti i file rimasti
# ------------------------------------------------------------------

files <- list.files(
  path       = path,
  pattern    = file_pattern,
  full.names = TRUE
)

temp_files <- file.path(
  dirname(files),
  sub("\\.csv$", "_temp.csv", basename(files))
)

renamed_to_temp <- file.rename(files, temp_files)

if (!all(renamed_to_temp)) {
  stop("Errore durante la rinomina temporanea di alcuni file.")
}

# ------------------------------------------------------------------
# 3. Rinumera le repliche:
#    sel[1] -> rep1
#    sel[2] -> rep2
#    ...
#    sel[20] -> rep20
# ------------------------------------------------------------------

for (i in seq_along(sel)) {

  old_rr <- sel[i]

  # Trova contemporaneamente Y, Delta e Lambda con questa replica
  current_files <- list.files(
    path       = path,
    pattern    = paste0(
      "_(Y|Delta|Lambda)_m[0-9]+_H[0-9]+_T[0-9]+_rep",
      old_rr,
      "_temp\\.csv$"
    ),
    full.names = TRUE
  )

  if (length(current_files) == 0) {
    warning("Nessun file trovato per rr = ", old_rr)
    next
  }

  new_files <- file.path(
    dirname(current_files),
    sub(
      pattern     = paste0("_rep", old_rr, "_temp\\.csv$"),
      replacement = paste0("_rep", i, ".csv"),
      x           = basename(current_files)
    )
  )

  renamed_final <- file.rename(current_files, new_files)

  if (!all(renamed_final)) {
    warning("Errore nella rinomina di alcuni file con rr = ", old_rr)
  }
}