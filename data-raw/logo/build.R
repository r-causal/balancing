# Run from the package root: Rscript data-raw/logo/build.R
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
source_dir <- dirname(normalizePath(script))
root <- normalizePath(file.path(source_dir, "../.."))
svg <- paste(readLines(file.path(source_dir, "sticker.svg")), collapse = "\n")
families <- regmatches(svg, gregexpr('font-family="[^"]+"', svg))[[1]]
families <- unique(sub('font-family="([^"]+)"', "\\1", families))
missing <- setdiff(families, systemfonts::system_fonts()$family)
if (length(missing)) {
  stop(
    "Install the source fonts before rebuilding: ",
    paste(missing, collapse = ", ")
  )
}
art <- base64enc::dataURI(
  file = file.path(source_dir, "artwork.png"),
  mime = "image/png"
)
svg <- sub("artwork.png", art, svg, fixed = TRUE)
exports <- file.path(source_dir, "exports")
dir.create(exports, showWarnings = FALSE)
print_png <- file.path(exports, "logo.png")
raster <- rsvg::rsvg_png(charToRaw(svg), width = 1044, height = 1200)
magick::image_write(
  magick::image_read(raster),
  path = print_png,
  format = "png",
  density = "600x600"
)

# Outline the text separately to avoid resampling the embedded artwork.
text_only <- gsub("<image[^>]*/>", "", svg)
text_only <- gsub("<rect[^>]*/>", "", text_only)
text_only <- gsub('<polygon id="sticker-border"[^>]*/>', "", text_only)
print_svg <- file.path(exports, "logo.svg")
rsvg::rsvg_svg(charToRaw(text_only), print_svg, width = 1740, height = 2000)
outlined <- paste(readLines(print_svg), collapse = "\n")
inner <- sub("(?s).*?<svg[^>]*>", "", outlined, perl = TRUE)
inner <- sub("</svg>\\s*$", "", inner, perl = TRUE)
without_text <- gsub("(?s)<text[^>]*>.*?</text>", "", svg, perl = TRUE)
combined <- sub("</svg>", paste0(inner, "</svg>"), without_text, fixed = TRUE)
stopifnot(!grepl("<text[ >]", combined))
writeLines(combined, print_svg)

usethis::proj_set(root)
options(usethis.overwrite = TRUE)
usethis::use_logo(print_png)
