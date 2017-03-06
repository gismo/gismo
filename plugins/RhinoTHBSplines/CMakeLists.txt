cmake_minimum_required(VERSION 2.8)

if (WIN32)
	project(THBSplines CXX)
	
	set(VERSION_MAJOR 1)
	set(VERSION_MINOR 0)
	set(VERSION_BUGFIX 0)

	if (${CMAKE_GENERATOR} MATCHES 2010)
		set (PLATFORM_TOOLSET v100)
	elseif (${CMAKE_GENERATOR} MATCHES 2015)
		set (PLATFORM_TOOLSET v140)
	endif()

	configure_file(Platform.targets.in ${CMAKE_BINARY_DIR}/Platform.targets)

	add_subdirectory(GismoSupport)
	add_subdirectory(GismoExport)
	add_subdirectory(GismoImport)
	add_subdirectory(THBSplines)

endif(WIN32)