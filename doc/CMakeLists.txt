######################################################################
## CMakeLists.txt ---
## This file is part of the G+Smo library.
## 
## Author: Angelos Mantzaflaris 
## Copyright (C) 2012-2015 - RICAM-Linz.
######################################################################

project(GismoDoc)

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

find_package(Doxygen)
if(DOXYGEN_FOUND)

project(GismoDoc)

# Detect GraphViz dot
find_program(DOT_EXECUTABLE dot)
if (DOT_EXECUTABLE)
    #message("Detected GraphViz dot")
    SET(HAVE_DOT YES)
else ()
    SET(HAVE_DOT NO)
endif (DOT_EXECUTABLE)

if(CMAKE_COMPILER_IS_GNUCXX)
  if(CMAKE_SYSTEM_NAME MATCHES Linux)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g1")
  endif(CMAKE_SYSTEM_NAME MATCHES Linux)
endif(CMAKE_COMPILER_IS_GNUCXX)

# Set some Doxygen flags
set(GISMO_DOXY_PROJECT_NAME             "G+Smo")
set(GISMO_DOXY_OUTPUT_DIRECTORY_SUFFIX  "")
set(GISMO_DOXY_INPUT                    "\"${gismo_SOURCE_DIR}/src\" \"${gismo_SOURCE_DIR}/doc\" \"${gismo_SOURCE_DIR}/plugins\" \"${gismo_SOURCE_DIR}/extensions\"")
set(IMAGE_PATH        "\"${gismo_SOURCE_DIR}/doc/figs\"")
set(GISMO_DOXY_HTML_COLORSTYLE_HUE      "120") #220
set(GISMO_DOXY_HTML_COLORSTYLE_SAT      "100") #100
set(GISMO_DOXY_HTML_COLORSTYLE_GAMMA    "100")  # 80
set(GISMO_DOXY_TAGFILES                 "")
set(SHOW_INCLUDE_FILES                  "yes")
set(CLASS_FLAGS                         "yes")
set(COLLABORATION_GRAPH                 "yes")
set(SHOW_USED_FILES                     "no")
set(INCLUDED_BY_GRAPH                   "no")
set(GROUP_GRAPHS                        "no")
set(UML_LOOK                            "no")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

If(GISMO_PLAINDOX)
 configure_file(
   ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_header_empty.html.in
   ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_header.html )
 configure_file(
   ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_footer_empty.html.in
   ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_footer.html )
else()
 configure_file(
   ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_header.html.in
   ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_header.html )
 configure_file(
   ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_footer.html.in
   ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_footer.html )
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_layout.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_layout.xml )

#add_definitions("-DGISMO_MAKING_DOCS")

set(snippets_targets "")
add_subdirectory(snippets)

add_custom_target(
  doc-gismo-prerequisites
  COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/navtree_hacks.js           ${CMAKE_CURRENT_BINARY_DIR}/html/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

if ( NOT NOSNIPPETS) # drop dependancy on snippets
add_dependencies(doc-gismo-prerequisites all_snippets)
endif()

add_custom_target(doc #ALL
  COMMAND doxygen
  WORKING_DIRECTORY ${gismo_BINARY_DIR}/doc
  COMMENT "Generating API documentation with Doxygen" VERBATIM 
)

add_custom_target(doctar #ALL
  COMMAND ${CMAKE_COMMAND} -E rename html gismo-doc
  COMMAND ${CMAKE_COMMAND} -E remove gismo-doc/gismo-doc.tgz
  COMMAND ${CMAKE_COMMAND} -E tar cvfz gismo-doc/gismo-doc.tgz gismo-doc
  COMMAND ${CMAKE_COMMAND} -E rename gismo-doc html
  WORKING_DIRECTORY ${gismo_BINARY_DIR}/doc
  COMMENT "Generating API documentation with Doxygen plus tar file" VERBATIM 
  DEPENDS doc
)

add_dependencies(doc doc-gismo-prerequisites)


endif(DOXYGEN_FOUND)
