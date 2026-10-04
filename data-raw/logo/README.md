# Package logo

The text-free artwork was generated with the built-in image-generation tool.
`prompt.txt` records the artwork prompt. `sticker.svg` adds editable typography
and a border clipped to the 1.74 by 2 inch hexagon.

Run `Rscript data-raw/logo/build.R` from the package root to export a 600 dpi
print PNG and an SVG with outlined text to `exports/`, then install the README
logo with `usethis::use_logo()`.

The build requires the R packages base64enc, magick, rsvg, systemfonts, and
usethis, plus Baskerville installed locally.
The artwork, fonts, and print tools are not package runtime dependencies.
