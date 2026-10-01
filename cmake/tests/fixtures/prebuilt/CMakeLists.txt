######################################################################
## fixtures/prebuilt/CMakeLists.txt
##
## Builds the "already installed" packages the *found* cases in
## run_dependency_tests.cmake resolve via find_package()/<Name>_DIR: a real
## compiled static library for FakeSrcDep (its Config.cmake's IMPORTED_LOCATION
## names a fixed, pre-pinned output path rather than a generator expression,
## so the same Config.cmake content is valid for every configuration on a
## multi-config generator), plus a battery of found-package shapes that only
## exist to be consumed by
## gismo_add_dependency() - multi-element legacy *_INCLUDE_DIRS/*_LIBRARIES,
## a target with no INTERFACE_INCLUDE_DIRECTORIES at all, and four
## generator-expression shapes in a target's INTERFACE_INCLUDE_DIRECTORIES.
## Configured and built exactly once, into a scratch PREFIX_DIR, before any
## of the found-based test cases configure. Nothing here touches the source
## tree: PREFIX_DIR, FIXTURES_DIR and DIRA/DIRB are all driver-supplied
## scratch paths.
######################################################################

cmake_minimum_required(VERSION 3.1...3.10)
project(gsdeptest_prebuilt CXX)

if(NOT DEFINED PREFIX_DIR)
  message(FATAL_ERROR "PREFIX_DIR must be set (-DPREFIX_DIR=...)")
endif()
if(NOT DEFINED FIXTURES_DIR)
  message(FATAL_ERROR "FIXTURES_DIR must be set (-DFIXTURES_DIR=...)")
endif()
if(NOT DEFINED DIRA OR NOT DEFINED DIRB)
  message(FATAL_ERROR "DIRA/DIRB must be set (-DDIRA=... -DDIRB=...)")
endif()

file(MAKE_DIRECTORY "${DIRA}")
file(MAKE_DIRECTORY "${DIRB}")

## --- FakeHdrDep: header-only package, target defined by its own config ----
## Instantiated from the shared PkgConfig.cmake.in template (the header-only
## variant every plain found-package case with no generator-expression
## include dir can reuse).
file(COPY "${FIXTURES_DIR}/header_only/include" DESTINATION "${PREFIX_DIR}/FakeHdrDep")
set(GSDEP_NAME FakeHdrDep)
set(GSDEP_INCLUDE_DIR "${PREFIX_DIR}/FakeHdrDep/include")
file(MAKE_DIRECTORY "${PREFIX_DIR}/FakeHdrDep/lib/cmake/FakeHdrDep")
configure_file("${FIXTURES_DIR}/PkgConfig.cmake.in"
  "${PREFIX_DIR}/FakeHdrDep/lib/cmake/FakeHdrDep/FakeHdrDepConfig.cmake" @ONLY)

## --- FakeSrcDep: real static library, target defined by its own config ----
## The archive output directory is pinned to the same path for every
## configuration (both the plain, single-config property and every
## per-config <CONFIG> variant), so the Config.cmake below can name the
## archive with a literal path built from CMAKE_STATIC_LIBRARY_PREFIX/SUFFIX
## instead of a $<TARGET_FILE:...> generator expression - file(GENERATE)
## evaluates genex per-configuration, and on a multi-config generator
## (Visual Studio, Ninja Multi-Config) writing the same output file once per
## configuration with different resolved content is a hard CMake error
## ("Evaluation file to be written multiple times with different content").
file(COPY "${FIXTURES_DIR}/sources/include" DESTINATION "${PREFIX_DIR}/FakeSrcDep")
add_library(gsdeptest_prebuilt STATIC "${FIXTURES_DIR}/sources/src/dep_lib.cpp")
target_include_directories(gsdeptest_prebuilt PUBLIC "${FIXTURES_DIR}/sources/include")
set(_gsdep_archive_dir "${PREFIX_DIR}/FakeSrcDep/lib")
set_target_properties(gsdeptest_prebuilt PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  ARCHIVE_OUTPUT_DIRECTORY "${_gsdep_archive_dir}")
foreach(_gsdep_cfg Debug Release RelWithDebInfo MinSizeRel)
  string(TOUPPER "${_gsdep_cfg}" _gsdep_cfg_upper)
  set_target_properties(gsdeptest_prebuilt PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY_${_gsdep_cfg_upper} "${_gsdep_archive_dir}")
endforeach()
set(_gsdep_archive_file
  "${_gsdep_archive_dir}/${CMAKE_STATIC_LIBRARY_PREFIX}gsdeptest_prebuilt${CMAKE_STATIC_LIBRARY_SUFFIX}")
file(WRITE "${PREFIX_DIR}/FakeSrcDep/lib/cmake/FakeSrcDep/FakeSrcDepConfig.cmake"
"add_library(GsDepTest::FakeSrcDep STATIC IMPORTED)
set_target_properties(GsDepTest::FakeSrcDep PROPERTIES
  IMPORTED_LOCATION \"${_gsdep_archive_file}\"
  INTERFACE_INCLUDE_DIRECTORIES \"${PREFIX_DIR}/FakeSrcDep/include\")
set(FakeSrcDep_FOUND TRUE)
")

## --- MultiTwiceDep: legacy <PKG>_INCLUDE_DIRS/_LIBRARIES, two elements
## each, no config-defined target - exercises per-element dedup of a
## multi-path found result across re-configures.
file(MAKE_DIRECTORY "${PREFIX_DIR}/MultiTwiceDep/include_a")
file(MAKE_DIRECTORY "${PREFIX_DIR}/MultiTwiceDep/include_b")
file(MAKE_DIRECTORY "${PREFIX_DIR}/MultiTwiceDep/lib")
file(WRITE "${PREFIX_DIR}/MultiTwiceDep/lib/liba.fake" "")
file(WRITE "${PREFIX_DIR}/MultiTwiceDep/lib/libb.fake" "")
file(WRITE "${PREFIX_DIR}/MultiTwiceDep/lib/cmake/MultiTwiceDep/MultiTwiceDepConfig.cmake"
"set(MULTITWICEDEP_INCLUDE_DIRS \"${PREFIX_DIR}/MultiTwiceDep/include_a\" \"${PREFIX_DIR}/MultiTwiceDep/include_b\")
set(MULTITWICEDEP_LIBRARIES \"${PREFIX_DIR}/MultiTwiceDep/lib/liba.fake\" \"${PREFIX_DIR}/MultiTwiceDep/lib/libb.fake\")
set(MultiTwiceDep_FOUND TRUE)
")

## --- NoIncludeDep: config-defined target with no INTERFACE_INCLUDE_DIRECTORIES
## at all - exercises the empty-include-dir guard on the found branch.
file(WRITE "${PREFIX_DIR}/NoIncludeDep/lib/cmake/NoIncludeDep/NoIncludeDepConfig.cmake"
"add_library(GsDepTest::NoIncludeDep INTERFACE IMPORTED)
set(NoIncludeDep_FOUND TRUE)
")

## --- GenexBuildMultiDep: one BUILD_INTERFACE genex wrapping two paths -
## CMake list-splits the property on ';' before gismo_add_dependency() ever
## sees it, so this exercises rejoining a genex CMake has already split
## across two raw property elements.
file(WRITE "${PREFIX_DIR}/GenexBuildMultiDep/lib/cmake/GenexBuildMultiDep/GenexBuildMultiDepConfig.cmake"
"add_library(GsDepTest::GenexBuildMultiDep INTERFACE IMPORTED)
set_target_properties(GsDepTest::GenexBuildMultiDep PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES \"$<BUILD_INTERFACE:${DIRA};${DIRB}>\")
set(GenexBuildMultiDep_FOUND TRUE)
")

## --- GenexConcatDep: BUILD_INTERFACE and INSTALL_INTERFACE concatenated in
## one element with no separator, as target_include_directories() called once
## with both arguments produces - exercises locating each genex token's own
## span within a single property element rather than one anchored match per
## element.
file(WRITE "${PREFIX_DIR}/GenexConcatDep/lib/cmake/GenexConcatDep/GenexConcatDepConfig.cmake"
"add_library(GsDepTest::GenexConcatDep INTERFACE IMPORTED)
set_target_properties(GsDepTest::GenexConcatDep PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES \"$<BUILD_INTERFACE:${DIRA}>$<INSTALL_INTERFACE:include>\")
set(GenexConcatDep_FOUND TRUE)
")

## --- GenexInstallOnlyDep: only an INSTALL_INTERFACE genex - nothing survives
## normalisation, so resolution must fail with "no include directory".
file(WRITE "${PREFIX_DIR}/GenexInstallOnlyDep/lib/cmake/GenexInstallOnlyDep/GenexInstallOnlyDepConfig.cmake"
"add_library(GsDepTest::GenexInstallOnlyDep INTERFACE IMPORTED)
set_target_properties(GsDepTest::GenexInstallOnlyDep PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES \"$<INSTALL_INTERFACE:include>\")
set(GenexInstallOnlyDep_FOUND TRUE)
")

## --- GenexUnknownDep: an unrecognised $<...> element alongside a real,
## plain path - the unrecognised element must survive, exempt from EXISTS,
## while the plain sibling is still checked normally.
file(WRITE "${PREFIX_DIR}/GenexUnknownDep/lib/cmake/GenexUnknownDep/GenexUnknownDepConfig.cmake"
"add_library(GsDepTest::GenexUnknownDep INTERFACE IMPORTED)
set_target_properties(GsDepTest::GenexUnknownDep PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES \"$<$<CONFIG:Debug>:${DIRB}>;${DIRA}\")
set(GenexUnknownDep_FOUND TRUE)
")
