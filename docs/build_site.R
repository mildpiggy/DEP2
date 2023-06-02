
usethis::use_pkgdown()

pkgdown::build_site(install = F)

pkgdown::build_articles()
pkgdown::build_article(name = "Data_import")

pkgdown::build_site_github_pages(install = F)
