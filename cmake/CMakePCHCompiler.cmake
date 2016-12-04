# CMakePCH
# --------
#
# Modified for G+Smo (Angelos Mantzaflaris)
#
# Author: Adam Strzelecki <ono@java.pl>
# Copyright (c) 2014-2015 Adam Strzelecki. All rights reserved.
# This code is licensed under the MIT License.

include(CMakeParseArguments)

# Create precompiled header
function(add_precompiled_header pch_target header)
       set(lang ${CMAKE_PCH_COMPILER_LANGUAGE})
       if(NOT MSVC AND
		NOT CMAKE_COMPILER_IS_GNU${lang} AND
		NOT CMAKE_${lang}_COMPILER_ID STREQUAL "GNU" AND
		NOT CMAKE_${lang}_COMPILER_ID STREQUAL "Clang" AND
		NOT CMAKE_${lang}_COMPILER_ID STREQUAL "AppleClang"
		)
		message(WARNING
			"Precompiled headers not supported for ${CMAKE_${lang}_COMPILER_ID}"
			)
		return()
	endif()

	if(ARGS_TYPE)
		set(header_type ${ARGS_TYPE})
	elseif(lang STREQUAL CXX)
		set(header_type "c++-header")
	elseif(lang STREQUAL C)
		set(header_type "c-header")
	else()
		message(WARNING "Unknown header type for language ${lang}")
		set(header_type "c++-header")
	endif()

	add_library(${pch_target} OBJECT ${header})
	set_target_properties(${pch_target} PROPERTIES
        COMPILE_DEFINITIONS gismo_EXPORTS
	POSITION_INDEPENDENT_CODE ON # bug: no effect, need fPIC later
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${pch_target}.dir")
	
	if(MSVC)
		# ensure pdb goes to the same location, otherwise we get C2859
		file(TO_NATIVE_PATH
			"${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${pch_target}.dir"
			pdb_dir
			)
		# /Yc - create precompiled header
		# /Fd - specify directory for pdb output
		set(flags "/Yc /Fd${pdb_dir}\\")
	else()
		set(flags "-x ${header_type}")
		set_target_properties(${pch_target} PROPERTIES
			COMPILE_FLAGS "-fPIC" )
	endif()

        set_source_files_properties(
		${header}
		PROPERTIES
		LANGUAGE ${lang}PCH
		COMPILE_FLAGS ${flags}
		)
endfunction()

# Use precompiled header
function(target_precompiled_header target pch_target)
       set(lang ${CMAKE_PCH_COMPILER_LANGUAGE})
       if(NOT MSVC AND
		NOT CMAKE_COMPILER_IS_GNU${lang} AND
		NOT CMAKE_${lang}_COMPILER_ID STREQUAL "GNU" AND
		NOT CMAKE_${lang}_COMPILER_ID STREQUAL "Clang" AND
		NOT CMAKE_${lang}_COMPILER_ID STREQUAL "AppleClang"
		)
		message(WARNING
			"Precompiled headers not supported for ${CMAKE_${lang}_COMPILER_ID}"
			)
		return()
	endif()
        add_dependencies(${target} ${pch_target})
	get_target_property(target_hdr ${pch_target} SOURCES)
	get_target_property(target_dir ${pch_target} LIBRARY_OUTPUT_DIRECTORY)
	if(MSVC)
		get_filename_component(win_header "${target_hdr}" NAME)
		file(TO_NATIVE_PATH "${target_dir}/${pch_header}" win_pch)
		# /Yu - use given include as precompiled header
		# /Fp - exact location for precompiled header
		# /FI - force include of precompiled header
		set(flags "/Yu${win_header} /Fp${win_pch} /FI${win_header}")
	else()
		#Note: ${target_dir}/${target_hdr} does not exist, so PCH is used
		set(flags "-Winvalid-pch -include ${target_dir}/${target_hdr}") # -H
	endif()
	set_target_properties(${target} PROPERTIES COMPILE_FLAGS "${flags}")
endfunction()

################################################################################
# PRIVATE MACROS
################################################################################

macro(__define_pch_compiler lang)
	if(NOT CMAKE_PCH_COMPILER_LANGUAGE)
		set(CMAKE_PCH_COMPILER_LANGUAGE ${lang})
	endif()

	# copy compiler settings from existing compiler
	set(CMAKE_${lang}PCH_COMPILE_OBJECT ${CMAKE_${lang}_COMPILE_OBJECT})
	set(CMAKE_INCLUDE_FLAG_${lang}PCH ${CMAKE_INCLUDE_FLAG_${lang}})
	set(CMAKE_INCLUDE_FLAG_SEP_${lang}PCH ${CMAKE_INCLUDE_FLAG_SEP_${lang}})

	if(CMAKE_COMPILER_IS_GNU${lang} OR
		CMAKE_${lang}_COMPILER_ID STREQUAL "GNU"
		)
		set(CMAKE_${lang}PCH_OUTPUT_EXTENSION .gch)
	else()
		set(CMAKE_${lang}PCH_OUTPUT_EXTENSION .pch)
	endif()

	# setup compiler & platform specific flags same way C/CXX does
	if(CMAKE_${lang}_COMPILER_ID)
		include(Platform/${CMAKE_SYSTEM_NAME}-${CMAKE_${lang}_COMPILER_ID}-${lang}PCH
			OPTIONAL
			)
	endif()

	# just use all settings from C/CXX compiler
	string(REPLACE "${lang}PCH" "${lang}"
		CMAKE_${lang}PCH_COMPILE_OBJECT
		${CMAKE_${lang}PCH_COMPILE_OBJECT}
		)

	if(MSVC)
		# redirect object file to NUL and just create precompiled header
		# /FoNUL - do not write output object file file
		# /Fp - specify location for precompiled header
		string(REPLACE " /Fo" " /FoNUL /Fp"
			CMAKE_${lang}PCH_COMPILE_OBJECT
			${CMAKE_${lang}PCH_COMPILE_OBJECT}
			)
		# disable pdb, we point to later to different location
		string(REPLACE " /Fd<TARGET_COMPILE_PDB>" ""
			CMAKE_${lang}PCH_COMPILE_OBJECT
			${CMAKE_${lang}PCH_COMPILE_OBJECT}
			)
	endif()

	# copy all initial settings for C/CXXPCH from C/CXX & watch them
	set(CMAKE_${lang}PCH_FLAGS "${CMAKE_${lang}_FLAGS_INIT}"
		CACHE STRING
		"Flags used by the compiler during all build types."
		)
	variable_watch(CMAKE_${lang}_FLAGS __watch_pch_variable)

	if(NOT CMAKE_NOT_USING_CONFIG_FLAGS)
		set(CMAKE_${lang}PCH_FLAGS_DEBUG "${CMAKE_${lang}_FLAGS_DEBUG_INIT}"
			CACHE STRING
			"Flags used by the compiler during debug builds."
			)
		set(CMAKE_${lang}PCH_FLAGS_MINSIZEREL "${CMAKE_${lang}_FLAGS_MINSIZEREL_INIT}"
			CACHE STRING
			"Flags used by the compiler during release builds for minimum size."
			)
		set(CMAKE_${lang}PCH_FLAGS_RELEASE "${CMAKE_${lang}_FLAGS_RELEASE_INIT}"
			CACHE STRING
			"Flags used by the compiler during release builds."
			)
		set(CMAKE_${lang}PCH_FLAGS_RELWITHDEBINFO "${CMAKE_${lang}_FLAGS_RELWITHDEBINFO_INIT}"
			CACHE STRING
			"Flags used by the compiler during release builds with debug info."
			)
		variable_watch(CMAKE_${lang}_FLAGS_DEBUG          __watch_pch_variable)
		variable_watch(CMAKE_${lang}_FLAGS_MINSIZEREL     __watch_pch_variable)
		variable_watch(CMAKE_${lang}_FLAGS_RELEASE        __watch_pch_variable)
		variable_watch(CMAKE_${lang}_FLAGS_RELWITHDEBINFO __watch_pch_variable)
	endif()
endmacro()


# copies all custom compiler settings to PCH compiler
macro(__watch_pch_variable variable access value)
	string(REPLACE _C_ _CPCH_ pchvariable ${variable})
	string(REPLACE _CXX_ _CXXPCH_ pchvariable ${pchvariable})
	set(${pchvariable} ${${variable}}) # because ${value} expands backslashes
endmacro()

macro(__configure_pch_compiler lang)
	set(CMAKE_${lang}PCH_COMPILER_ENV_VAR "${lang}PCH")
	set(CMAKE_${lang}PCH_COMPILER ${CMAKE_${lang}_COMPILER})

	if(SET_MSVC_${lang}PCH_ARCHITECTURE_ID)
		string(REPLACE _${lang}_ _${lang}PCH_
			${SET_MSVC_${lang}_ARCHITECTURE_ID}
			SET_MSVC_${lang}PCH_ARCHITECTURE_ID
			)
	endif()
	if(CMAKE_${lang}_SYSROOT_FLAG_CODE)
		string(REPLACE _${lang}_ _${lang}PCH_
			${CMAKE_${lang}_SYSROOT_FLAG_CODE}
			CMAKE_${lang}PCH_SYSROOT_FLAG_CODE
			)
	endif()
	if(CMAKE_${lang}_OSX_DEPLOYMENT_TARGET_FLAG_CODE)
		string(REPLACE _${lang}_ _${lang}PCH_
			${CMAKE_${lang}_OSX_DEPLOYMENT_TARGET_FLAG_CODE}
			CMAKE_${lang}PCH_OSX_DEPLOYMENT_TARGET_FLAG_CODE
			)
	endif()

	configure_file(
		${CMAKE_CURRENT_LIST_DIR}/CMake${lang}PCHCompiler.cmake.in
		${CMAKE_PLATFORM_INFO_DIR}/CMake${lang}PCHCompiler.cmake
		)
endmacro()


# Defines compiler command templates based on existing C++ compiler.
__define_pch_compiler(CXX)
