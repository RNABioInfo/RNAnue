option(CPM_USE_LOCAL_PACKAGES
       "Try `find_package` before downloading dependencies" ON)

CPMAddPackage(
    NAME matplotplusplus
    GITHUB_REPOSITORY alandefreitas/matplotplusplus
    GIT_TAG origin/master # or whatever tag you want
)
