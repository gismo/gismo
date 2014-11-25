cmake_minimum_required(VERSION 2.8)
project(CPPLOT)

find_package(GLUT)
find_package(OpenGL)

find_package( Boost COMPONENTS thread )
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories ( ${Boost_INCLUDE_DIRS} )

#set(CMAKE_CXX_FLAGS "-ggdb -Wall -Wextra -fPIC")

## Build the library
set(SRC_DIR ".")
add_library(cpplot STATIC ${SRC_DIR}/cpplot.cpp ${SRC_DIR}/math.cpp ${SRC_DIR}/glut.cpp ${SRC_DIR}/figure.cpp ${SRC_DIR}/layer.cpp ${SRC_DIR}/axes.cpp ${SRC_DIR}/color.cpp ${SRC_DIR}/drawing.cpp ${SRC_DIR}/line.cpp ${SRC_DIR}/patch.cpp ${SRC_DIR}/surface.cpp ${SRC_DIR}/text.cpp ${SRC_DIR}/gl2ps.cpp)
target_link_libraries(cpplot
    ${GLUT_LIBRARY}
    ${OPENGL_LIBRARY}
    ${Boost_LIBRARIES}
    #boostpart   # Using local copy of part of boost libraries
    )
