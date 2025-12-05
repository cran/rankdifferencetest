# Version number can't have leading zeros
# 1. Update Description
# 2. Update NEWS
# 3. Update README

devtools::document()
devtools::check(remote = TRUE)

# Check pdf manual
devtools::build_manual()

# Check package coverage
covr::package_coverage()
covr::report()

# goodpractice
goodpractice::gp()

# Check spelling
spelling::spell_check_package()

# Check URLs are correct
# install.packages('urlchecker', repos = 'https://r-lib.r-universe.dev')
urlchecker::url_check()

# updates any URLs which are permanent (301) redirects.
#urlchecker::url_update()

# Check on CRAN
# _win devel CRAN
#devtools::check_win_devel()
# _win release CRAN
#devtools::check_win_release()
# _macos CRAN
# Need to follow the URL proposed to see the results
#devtools::check_mac_release()

devtools::revdep()

# Release
#devtools::release()
