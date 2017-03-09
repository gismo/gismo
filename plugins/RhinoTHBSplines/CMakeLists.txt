cmake_minimum_required(VERSION 2.8)

if (WIN32)
	project(THBSplines CXX)

	if(NOT "x${CMAKE_GENERATOR}" STREQUAL "xVisual Studio 14 2015 Win64")
	message(ERROR "Generator is NOT Visual Studio 14 2015 Win64")
	endif()

	if(NOT ("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "Debug") )
	message(ERROR "Build type is neither Debug nor Release")
	endif()
	
	#set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/plRhino)
	
	set(VERSION_MAJOR 1)
	set(VERSION_MINOR 0)
	set(VERSION_BUGFIX 0)

	if (${CMAKE_GENERATOR} MATCHES 2010)
		set (PLATFORM_TOOLSET v100)
	elseif (${CMAKE_GENERATOR} MATCHES 2015)
		set (PLATFORM_TOOLSET v140)
	endif()

	configure_file(Platform.targets.in ${CMAKE_BINARY_DIR}/Platform.targets)

	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/GismoSupport)
	
	add_subdirectory(GismoSupport)
	add_subdirectory(GismoExport)
	add_subdirectory(GismoImport)
	add_subdirectory(THBSplines)

endif(WIN32)
