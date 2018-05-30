include(ExternalProject)

set(AngpowTag v0.4)

# Downloads and compiles Angpow
ExternalProject_Add(Angpow
        PREFIX ANGPOW
        GIT_REPOSITORY https://github.com/LSSTDESC/Angpow4CCL.git
        GIT_TAG ${AngpowTag}
        DOWNLOAD_NO_PROGRESS 1
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern)
set(ANGPOW_LIBRARY_DIRS ${CMAKE_BINARY_DIR}/extern/lib/ )
set(ANGPOW_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/extern/include/Angpow)
set(ANGPOW_LIBRARIES -langpow -lstdc++)
