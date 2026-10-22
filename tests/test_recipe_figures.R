source("scripts/render_recipe_figures.R")
output_dir <- tempfile()
render_recipe_figures(output_dir)
expected <- paste0(tools::file_path_sans_ext(recipe_files), ".png")
files <- sort(list.files(output_dir, pattern = "[.]png$"))
stopifnot(identical(files, sort(expected)), length(files) == 14L,
          all(file.info(file.path(output_dir, files))$size > 4000L))
headers <- vapply(file.path(output_dir, files), function(path) {
  paste(readBin(path, "raw", n = 8L), collapse = " ")
}, character(1))
stopifnot(all(headers == "89 50 4e 47 0d 0a 1a 0a"))
committed_previews <- c("workflows/rnaseq/example/pca.png",
                        "workflows/rnaseq/example/volcano.png",
                        "workflows/rnaseq/real_fastq/read_lengths.png")
stopifnot(all(file.exists(committed_previews)),
          all(file.info(committed_previews)$size > 4000L))
unlink(output_dir, recursive = TRUE)
cat("All recipe figures rendered successfully.\n")
