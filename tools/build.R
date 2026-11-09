usethis::use_mit_license(copyright_holder = "Akira Terui")
usethis::use_roxygen_md()
devtools::document()
devtools::load_all()
devtools::check(vignettes=FALSE)

lintr::lint_package()

## packagedown
# usethis::use_pkgdown()
# usethis::create_github_token()
# gitcreds::gitcreds_set()
# usethis::use_pkgdown_github_pages()

pkgdown::build_site()
pkgdown::preview_site()
