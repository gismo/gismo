######################################################################
## CMakeLists.txt --- doc
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project(GismoDoc)

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

find_package(Doxygen)
if(DOXYGEN_FOUND)

  # Detect GraphViz dot
  find_program(DOT_EXECUTABLE dot)
  if (DOT_EXECUTABLE)
    #message("Detected GraphViz dot")
    SET(HAVE_DOT YES)
  else ()
    SET(HAVE_DOT NO)
  endif (DOT_EXECUTABLE)

#  if(CMAKE_COMPILER_IS_GNUCXX)
#    if(CMAKE_SYSTEM_NAME MATCHES Linux)
#      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -g1")
#    endif(CMAKE_SYSTEM_NAME MATCHES Linux)
#  endif(CMAKE_COMPILER_IS_GNUCXX)

  # Set some Doxygen flags
  set(GISMO_DOXY_PROJECT_NAME             "G+Smo")
  set(GISMO_DOXY_OUTPUT_DIRECTORY_SUFFIX  "")
  set(GISMO_DOXY_INPUT                    "\"${gismo_SOURCE_DIR}/src\" \"${gismo_SOURCE_DIR}/doc\" \"${gismo_SOURCE_DIR}/extensions\" \"${gismo_BINARY_DIR}/gsCore\"") #\"${gismo_SOURCE_DIR}/plugins\" \"${gismo_SOURCE_DIR}/examples\"
#  set(GISMO_DOXY_INPUT                    \"${gismo_SOURCE_DIR}/extensions\")
  # Uncomment the line below to have doxygen process only
  # what is inside the doc folder when typing 'make doc'
  # This is ment to improve doc generation time. However, a very good way to
  # tweak doxygen temporarily is to go to CMAKE_BINARY_DIR/doc and edit
  # the Doxyfile directly. Then running doxygen inside that folder will
  # re-generate the doc html pages.
  #set(GISMO_DOXY_INPUT                    " \"${gismo_SOURCE_DIR}/doc\" \"${gismo_SOURCE_DIR}/extensions\" ")
  set(IMAGE_PATH        "\"${gismo_SOURCE_DIR}/doc/figs\"")
  set(GISMO_DOXY_EXCLUDE                    "\"${gismo_SOURCE_DIR}/external\" ${IMAGE_PATH} ")
  set(GISMO_DOXY_HTML_COLORSTYLE_HUE      "120") #220
  set(GISMO_DOXY_HTML_COLORSTYLE_SAT      "100") #100
  set(GISMO_DOXY_HTML_COLORSTYLE_GAMMA    "100")  # 80
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Eigen.doxytags")
    set(GISMO_DOXY_TAGFILES                 "${CMAKE_CURRENT_SOURCE_DIR}/Eigen.doxytags=http://eigen.tuxfamily.org/dox/")
  else()
    set(GISMO_DOXY_TAGFILES                 "")
  endif()
  set(SHOW_INCLUDE_FILES                  "yes")
  set(CLASS_FLAGS                         "yes")
  set(COLLABORATION_GRAPH                 "yes")
  set(SHOW_USED_FILES                     "no")
  set(INCLUDED_BY_GRAPH                   "no")
  set(GROUP_GRAPHS                        "no")
  set(UML_LOOK                            "no")

  # get the date in variable TODAY
  if(CMAKE_MAJOR_VERSION EQUAL 2 AND
      CMAKE_MINOR_VERSION LESS  9 AND
      CMAKE_PATCH_VERSION LESS 11 )
    if (WIN32)
      execute_process(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE TODAY)
      string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" TODAY ${TODAY})
    elseif(UNIX)
      execute_process(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE TODAY)
      string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" TODAY ${TODAY})
    else(WIN32)
      message(WARNING "date not implemented")
      set(TODAY 000000)
    endif(WIN32)
  else()
    string(TIMESTAMP TODAY "%d/%m/%Y")
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    )

  If(GISMO_PLAINDOX)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_header_empty.html.in
      ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_header.html COPYONLY)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_footer_empty.html.in
      ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_footer.html COPYONLY)
  else()
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_header.html.in
      ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_header.html COPYONLY)
    configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_footer.html.in
      ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_footer.html COPYONLY)
  endif()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/gismodoxy_layout.xml.in
    ${CMAKE_CURRENT_BINARY_DIR}/gismodoxy_layout.xml COPYONLY)

  #add_definitions("-DGISMO_MAKING_DOCS")

  add_custom_target(
    doc-gismo-prerequisites
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/html/
    #COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/navtree_hacks.js           ${CMAKE_CURRENT_BINARY_DIR}/html/
    DEPENDS ${snippets_out}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  if (NOT NOSNIPPETS)
    add_dependencies(doc-gismo-prerequisites doc-snippets)
  endif()

  add_custom_target(doc #ALL
    COMMAND doxygen
    WORKING_DIRECTORY ${gismo_BINARY_DIR}/doc
    COMMENT "Generating doc: ${gismo_BINARY_DIR}/doc/html/index.html" VERBATIM
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

else (DOXYGEN_FOUND)

  add_custom_target(doc
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Doxygen was not found! Please install doxygen and run cmake again to build the doc.")

endif(DOXYGEN_FOUND)
