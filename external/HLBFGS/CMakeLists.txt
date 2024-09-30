file(GLOB SOURCES *.cpp)
message(Calling LBFGS CMAKE)
add_library(hlbfgs ${SOURCES})
set_target_properties(hlbfgs PROPERTIES SOVERSION 1.2)
