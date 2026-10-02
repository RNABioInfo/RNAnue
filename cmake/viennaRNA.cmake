include(ExternalProject)

set(VIENNA_RNA_PREFIX ${CMAKE_BINARY_DIR}/submodules/viennaRNA-prefix)
set(VIENNA_RNA_INSTALL ${CMAKE_BINARY_DIR}/submodules/viennaRNA-install)
set(VIENNA_RNA_LIBRARY ${VIENNA_RNA_INSTALL}/lib/libRNA.a)

ExternalProject_Add(
    viennaRNA
    PREFIX ${VIENNA_RNA_PREFIX}
    DOWNLOAD_EXTRACT_TIMESTAMP true
    URL
        "https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.6.4/ViennaRNA-2.6.4.tar.gz"
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND
        ./configure CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
        --prefix=${VIENNA_RNA_INSTALL} --includedir=${VIENNA_RNA_INSTALL}/include
        --libdir=${VIENNA_RNA_INSTALL}/lib --disable-openmp --without-lto
        --without-gsl --without-perl --without-python --without-doc
        --without-kinfold --without-forester --without-kinwalker --without-cla
        --without-cla-pdf --without-doc-html --without-doc-pdf
        --without-tutorial --without-tutorial-html --without-tutorial-pdf
        --without-svm --without-swig --without-rnaxplorer --without-man
        --without-rnalocmin --disable-unittests --disable-check-executables
        --disable-check-perl --disable-check-python --disable-check-python2
    BUILD_COMMAND make -j 2 CXX=${CMAKE_CXX_COMPILER} CC=${CMAKE_C_COMPILER}
    INSTALL_COMMAND ""
)

# Use a custom install step with byproducts for CMake 3.24 compatibility.
# INSTALL_BYPRODUCTS requires CMake 3.26 or newer.
ExternalProject_Add_Step(viennaRNA install_library
    COMMAND make -j 2 install prefix=${VIENNA_RNA_INSTALL}
    WORKING_DIRECTORY <SOURCE_DIR>
    DEPENDEES build
    DEPENDERS install
    BYPRODUCTS "${VIENNA_RNA_LIBRARY}"
)

set(VIENNA_RNA_LIBRARIES ${VIENNA_RNA_LIBRARY})
set(VIENNA_RNA_INCLUDE_DIR ${VIENNA_RNA_INSTALL}/include)
