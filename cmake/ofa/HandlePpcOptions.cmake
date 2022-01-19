#=============================================================================
# Handling of PPC / PPC64 options
#
# This is a two-step process:
#
# 1. Generate a list of compiler flags for the specific CPU
#
# 2. Special compiler-specific treatment of "native" flag
#
# 3. Disabling of "broken" features based on OFA_xxx_INTRINSICS_BROKEN options
#
# 4. Set compiler-specific flags
#=============================================================================

include(ofa/AddCompilerFlag)
include(ofa/CommonMacros)
include(CheckIncludeFileCXX)

macro(OFA_HandlePpcOptions)
  set(_march_flag_list)
  set(_available_vector_units_list)
    
  # Define macros for PowerPC64
  macro(_power3)
    list(APPEND _march_flag_list "power3")
  endmacro()
  macro(_power4)
    list(APPEND _march_flag_list "power4")
    _power3()
  endmacro()
  macro(_power5)
    list(APPEND _march_flag_list "power5")
    _power4()
  endmacro()
  macro(_power5plus)
    list(APPEND _march_flag_list "power5+")
    _power5()
  endmacro()
  macro(_power6)
    list(APPEND _march_flag_list "power6")
    _power5()
  endmacro()
  macro(_power6x)
    list(APPEND _march_flag_list "power6x")
    _power6()
  endmacro()
  macro(_power7)
    list(APPEND _march_flag_list "power7")
    _power6()
  endmacro()
  macro(_power8)
    list(APPEND _march_flag_list "power8")
    list(APPEND _march_flag_list "pwr8")
    _power7()
  endmacro()
  macro(_power9)
    list(APPEND _march_flag_list "power9")
    list(APPEND _march_flag_list "pwr9")
    _power8()
  endmacro()
  macro(_power10)
    list(APPEND _march_flag_list "power10")
    list(APPEND _march_flag_list "pwr10")
    _power9()
  endmacro()
  
  # PowerPC64
  if(TARGET_ARCHITECTURE STREQUAL "power3")
    _power3()
  elseif(TARGET_ARCHITECTURE STREQUAL "power4")
    _power4()
  elseif(TARGET_ARCHITECTURE STREQUAL "power5")
    _power5()
  elseif(TARGET_ARCHITECTURE STREQUAL "power5+")
    _power5plus()
  elseif(TARGET_ARCHITECTURE STREQUAL "power6")
    _power6()
  elseif(TARGET_ARCHITECTURE STREQUAL "power6x")
    _power6x()
  elseif(TARGET_ARCHITECTURE STREQUAL "power7")
    _power7()
  elseif(TARGET_ARCHITECTURE STREQUAL "power8")
    _power8()
  elseif(TARGET_ARCHITECTURE STREQUAL "power9")
    _power9()
  elseif(TARGET_ARCHITECTURE STREQUAL "power10")
    _power10()

  # Others
  elseif(TARGET_ARCHITECTURE STREQUAL "generic")
    list(APPEND _march_flag_list "generic")
  elseif(TARGET_ARCHITECTURE STREQUAL "native")
    list(APPEND _march_flag_list "native")
  elseif(TARGET_ARCHITECTURE STREQUAL "none")
    # add this clause to remove it from the else clause

  else()
    message(FATAL_ERROR "Unknown target architecture: \"${TARGET_ARCHITECTURE}\". Please set TARGET_ARCHITECTURE to a supported value.")
  endif()

  # Special treatment for "native"
  if(TARGET_ARCHITECTURE STREQUAL "native")

  # Apply architecture flags
  elseif(NOT TARGET_ARCHITECTURE STREQUAL "none")

    # Disable "broken" features based on OFA_xxx_INTRINSICS_BROKEN options
    set(_disable_vector_unit_list)
    set(_enable_vector_unit_list)

    # Enable/disable macro
    macro(_enable_or_disable _name _flag _documentation _broken)
      if(_broken)
        set(_found false)
      else()
        _ofa_find(_available_vector_units_list "${_flag}" _found)
      endif()
      set(USE_${_name} ${_found} CACHE BOOL "${documentation}" ${_force})
      mark_as_advanced(USE_${_name})
      if(USE_${_name})
        list(APPEND _enable_vector_unit_list "${_flag}")
      else()
        list(APPEND _disable_vector_unit_list "${_flag}")
      endif()
    endmacro()

    # Enable/disable features
    
    # Add compiler flags
    if(CMAKE_CXX_COMPILER_ID MATCHES "NVHPC")

    elseif(CMAKE_CXX_COMPILER_ID MATCHES "XL")

    else()
      # Others: GNU, Clang and variants
      
      
    endif()
  endif()
endmacro(OFA_HandlePpcOptions)
