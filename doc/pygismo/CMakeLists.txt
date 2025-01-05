######################################################################
## CMakeLists.txt --- pyGismoDoc
## This file is part of the G+Smo library.
##
## Author: Angelos Mantzaflaris
######################################################################

project(doc-pygismo)

set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL TRUE)

find_program(SPHINX_BUILD NAMES sphinx-build sphinx-build.exe)
# unchecked: sphinx_rtd_theme is installed

if(SPHINX_BUILD AND GISMO_WITH_PYBIND11)

    configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
    ${CMAKE_CURRENT_BINARY_DIR}/conf.py
    )

  add_custom_target(
    doc-pygismo
    COMMAND ${SPHINX_BUILD} -c ${CMAKE_CURRENT_BINARY_DIR} -b html -d ${CMAKE_CURRENT_BINARY_DIR}/doctrees . ${gismo_BINARY_DIR}/doc/html/pygismo
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMENT "Generating pygismo-doc: ${gismo_BINARY_DIR}/doc/html/pygismo/index.html" VERBATIM
    )

  add_dependencies(doc doc-pygismo)
  add_dependencies(doc-pygismo pygismo)

else (SPHINX_BUILD)

  add_custom_target(doc-pygismo
    ${CMAKE_COMMAND} -E cmake_echo_color --cyan "Documentation for pygismo cannot be compiled! Please install sphinx, sphinx_rtd_theme -- http://sphinx-doc.org , or python3-sphinx in package manager, or sphinx in pip -- and enable GISMO_WITH_PYBIND11 to build the pygismo  documentation."
    )

endif()
