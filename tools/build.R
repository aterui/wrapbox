usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)

lintr::lint_package()

## packagedown
usethis::use_pkgdown()
pkgdown::build_site()
pkgdown::preview_site()
