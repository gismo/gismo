#####################################################################
## CMakeLists.txt --- cppyy
## This file is part of the G+Smo library.
##
## Author: Felix Scholz
######################################################################

set(GISMO_LIB_DIR  "${PROJECT_BINARY_DIR}/${LIB_INSTALL_DIR}")
set(PYTHON_PACKAGE_DIR ${PROJECT_BINARY_DIR}/cppyy)

set(PKG_NAME "gismo_cppyy")

configure_file(__init__.py.in ${PYTHON_PACKAGE_DIR}/${PKG_NAME}/__init__.py)
configure_file(pythonize.py.in ${PYTHON_PACKAGE_DIR}/${PKG_NAME}/pythonize.py)
configure_file(pyproject.toml.in ${PYTHON_PACKAGE_DIR}/pyproject.toml)
configure_file(setup.cfg.in ${PYTHON_PACKAGE_DIR}/setup.cfg)

set(PY_WHEEL_FILE "${PYTHON_PACKAGE_DIR}/dist/${PKG_NAME}-${gismo_VERSION}-py3-none-any.whl")

message("-- Python package was generated in ${PYTHON_PACKAGE_DIR}.
  To generate the pip wheel, run target \"wheel\".
  To install the bindings, run target \"install-bindings\".")

# Check if python-build is installed
execute_process(COMMAND python -c "import build"
        RESULT_VARIABLE EXIT_CODE
        OUTPUT_QUIET)
if(NOT ${EXIT_CODE} EQUAL 0)
    message("Python package \"build\" not found. If you want to generate the wheel, please install it with \"pip install build\"")
endif ()

add_custom_command(OUTPUT ${PY_WHEEL_FILE}
        COMMAND python -m build --wheel
        WORKING_DIRECTORY ${PYTHON_PACKAGE_DIR}
        COMMENT "  Generating wheel for the Python bindings...")
add_custom_target(wheel ALL DEPENDS ${PY_WHEEL_FILE})

add_custom_target(install-bindings COMMAND pip install ${PY_WHEEL_FILE}
        COMMENT "Installing the cppyy Python bindings using pip"
        DEPENDS wheel gismo)

install(CODE "execute_process(COMMAND pip install ${PY_WHEEL_FILE}")
