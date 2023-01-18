# Determine the host CPU feature set and determine the best set of compiler
# flags to enable all supported SIMD relevant features. Alternatively, the
# target CPU can be explicitly selected (for generating more generic binaries
# or for targeting a different system).
# Compilers provide e.g. the -march=native flag to achieve a similar result.
# This fails to address the need for building for a different microarchitecture
# than the current host.
# The script tries to deduce all settings from the model and family numbers of
# the CPU instead of reading the CPUID flags from e.g. /proc/cpuinfo. This makes
# the detection more independent from the CPUID code in the kernel (e.g. avx2 is
# not listed on older kernels).
#
# Usage:
# OptimizeForArchitecture()
#
# Optional inputs:
# TARGET_ARCHITECTURE=<name> specifies the target architecture (default=auto)
# TARGET_PROFILER=<name>     specifies the target profiler     (default=none)
# OFA_VERBOSE=<bool>         prints verbose output             (default=off)
#
# If any of the <feature>_broken flags are defined and set to true,
# the OptimizeForArchitecture macro will consequently disable the
# relevant features via compiler flags.
#
# Output:
# ARCHITECTURE_CXX_FLAGS compiler flags optimized for the target architecture
#
# Internal variables:
# USE_<feature>          boolean variable holding the status of <feature>
# HAVE_<feature>         boolean variable holding the compiler;s capability

#=============================================================================
# Copyright 2010-2016 Matthias Kretz <kretz@kde.org>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the names of contributing organizations nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

#=============================================================================
# Extension of the original version by Matthias Moller <m.moller@tudelft.nl>
#
# Changelog:
# - Update of CPUIDs for latest Intel and AMD processors
# - Added support for PPC64 (Clang, GCC, IBM XLC)
# - Added Support for ARM (Clang, GCC, ARM Clang, Cray, Fujitsu)
# - Restructuring and splitting into multiple files
#=============================================================================

#=============================================================================
# Autodetection of CPU
#=============================================================================

macro(OFA_AutodetectHostArchitecture)
  set(TARGET_ARCHITECTURE "generic")
  set(ARCHITECTURE_CXX_FLAGS CACHE STRING "CPU architecture compiler flags")
  
  if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
    include(ofa/AutodetectX86)
    OFA_AutodetectX86()
  elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "^(arm.*|ARM.*|aarch64.*|AARCH64.*)")
    include(ofa/AutodetectArm)
    OFA_AutodetectArm()
  elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "^(powerpc|ppc)64.*")
    include(ofa/AutodetectPpc)
    OFA_AutodetectPpc()
  else()
    message(WARNING "The CMAKE_SYSTEM_PROCESSOR '${CMAKE_SYSTEM_PROCESSOR}' is not supported by OptimizeForArchitecture")
  endif()
endmacro(OFA_AutodetectHostArchitecture)

#=============================================================================
# Handling of CPU options
#=============================================================================

macro(OptimizeForArchitecture)
  if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
    set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for. Using an incorrect setting here can result in crashes of the resulting binary because of invalid instructions used. Setting the value to \"auto\" will try to optimize for the architecture where cmake is called. Setting the value to \"native\" bypasses all checks and uses \"-march=native\" or the compiler equivalent flag. Other supported values are: \"none\", \"generic\", \"core\", \"core2\", \"merom\" (65nm Core2), \"penryn\" (45nm Core2), \"nehalem\", \"westmere\", \"sandybridge\", \"ivybridge\", \"haswell\", \"broadwell\", \"skylake\", \"skylake-xeon\", \"kabylake\", \"cannonlake\", \"cascadelake\", \"cooperlake\", \"icelake\", \"icelake-xeon\", \"tigerlake\", \"alderlake\", \"sapphirerapids\", \"bonnell\", \"silvermont\", \"goldmont\", \"goldmont-plus\", \"tremont\", \"knl\" (Knights Landing), \"knm\" (Knights Mill), \"atom\", \"k8\", \"k8-sse3\", \"barcelona\", \"istanbul\", \"magny-cours\", \"bulldozer\", \"interlagos\", \"piledriver\", \"steamroller\", \"excavator\", \"amd14h\", \"amd16h\", \"zen\", \"zen2\", \"zen3\"." )
  elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "^(arm.*|ARM.*|aarch64.*|AARCH64.*)")
    set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for. Using an incorrect setting here can result in crashes of the resulting binary because of invalid instructions used. Setting the value to \"auto\" will try to optimize for the architecture where cmake is called. Setting the value to \"native\" bypasses all checks and uses \"-march=native\" or the compiler equivalent flag. Other supported values are: \"none\", \"generic\", \"a64fx\", \"apple-a6\", \"apple-a7\", \"apple-a8\", \"apple-a9\", \"apple-a10\", \"apple-a11\", \"apple-a12\", \"apple-a13\", \"apple-a14\", \"apple-a15\", \"apple-a16\", \"apple-m1\", \"apple-m2\", \"arm1020e\", \"arm1020t\", \"arm1022e\", \"arm1026ej-s\", \"arm10e\", \"arm10tdmi\", \"arm1136j-s\", \"arm1136jf-s\", \"arm1156t2-s\", \"arm1156t2f-s\", \"arm1176jz-s\", \"arm1176jzf-s\", \"arm710t\", \"arm720t\", \"arm740t\", \"arm7tdmi-s\", \"arm7tdmi\", \"arm810\", \"arm8\", \"arm920\", \"arm920t\", \"arm922t\", \"arm926ej-s\", \"arm940t\", \"arm946e-s\", \"arm966e-s\", \"arm968e-s\", \"arm9\", \"arm9e\", \"arm9tdmi\", \"brahma-b15\", \"brahma-b53\", \"carmel\", \"cortex-a7\", \"cortex-a8\", \"cortex-a9\", \"cortex-a12\", \"cortex-a15.cortex-a7\", \"cortex-a15\", \"cortex-a17.cortex-a7\", \"cortex-a17\", \"cortex-a32\", \"cortex-a34\", \"cortex-a35\", \"cortex-a53\", \"cortex-a55\", \"cortex-a57.cortext-a53\", \"cortex-a57\", \"cortex-a5\", \"cortex-a72.cortext-a53\", \"cortex-a72\", \"cortex-a73.cortext-a35\", \"cortex-a73.cortext-a53\", \"cortex-a73\", \"cortex-a75.cortext-a55\", \"cortex-a75\", \"cortex-a76.cortext-a55\", \"cortex-a76\", \"cortex-a76ae\", \"cortex-a77\", \"cortex-a78\", \"cortex-a78ae\", \"cortex-a76c\", \"cortex-a510\", \"cortex-a710\", \"cortex-m0\", \"cortex-m0plus\", \"cortex-m1\", \"cortex-m23\", \"cortex-m33\", \"cortex-m35p\", \"cortex-m3\", \"cortex-m4\", \"cortex-m55\", \"cortex-m7\", \"cortex-r4\", \"cortex-r4f\", \"cortex-r52\", \"cortex-r5\", \"cortex-r7\", \"cortex-r8\", \"cortex-x1\", \"cortex-x2\", \"denver2\", \"denver\", \"exynos-m1\", \"fa526\", \"fa606te\", \"fa626\", \"fa626te\", \"fa726te\", \"falkor\", \"fmp626\", \"generic-armv7-a\", \"i80200\", \"i80321-400-b0\", \"i80321-400\", \"i80321-600-b0\", \"i80321-600\", \"ipx1200\", \"ipx425-266\", \"ipx425-400\", \"ipx425-533\", \"iwmmxt2\", \"iwmmxt\", \"krait\", \"kryo2\", \"kryo\", \"marvell-f\", \"marvell-pj4\", \"mpcore\", \"neoverse-e1\", \"neoverse-n1\", \"neoverse-n2\", \"neoverse-v1\", \"pxa210a\", \"pxa210b\", \"pxa210c\", \"pxa250a\", \"pxa250b\", \"pxa250c\", \"pxa27x\", \"pxa30x\", \"pxa31x\", \"pxa32x\", \"pxa930\", \"sa1110\", \"saphira\", \"scorpion\", \"strongarm1100\", \"strongarm110\", \"strongarm\", \"thunderx2\", \"thunderx2t99\", \"thunderx\", \"thunderxt81\", \"thunderxt83\", \"thunderxt88\", \"tsv110\", \"xgene1\", \"xscale\".")
  elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "^(powerpc|ppc)64.*")
    set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for. Using an incorrect setting here can result in crashes of the resulting binary because of invalid instructions used. Setting the value to \"auto\" will try to optimize for the architecture where cmake is called. Other supported values are: \"none\", \"generic\", \"power8\", \"power9\", \"power10\".")
  else()
    message(WARNING "The CMAKE_SYSTEM_PROCESSOR '${CMAKE_SYSTEM_PROCESSOR}' is not supported by OptimizeForArchitecture")
  endif()

  if(NOT OFA_VERBOSE)
    set(CMAKE_REQUIRED_QUIET true)
  endif()
  
  set(_force)
  if(NOT _last_target_arch STREQUAL "${TARGET_ARCHITECTURE}")
    message(STATUS "Target architecture changed from \"${_last_target_arch}\" to \"${TARGET_ARCHITECTURE}\"")
    set(_force FORCE)
  endif()
  set(_last_target_arch "${TARGET_ARCHITECTURE}" CACHE STRING "" FORCE)
  mark_as_advanced(_last_target_arch)
  string(TOLOWER "${TARGET_ARCHITECTURE}" TARGET_ARCHITECTURE)

  if(TARGET_ARCHITECTURE STREQUAL "auto")
    OFA_AutodetectHostArchitecture()
    message(STATUS "Detected Host CPU: ${TARGET_ARCHITECTURE}")
  endif()

  message(STATUS "Checking Host CPU features. This can take some time ...")
  if("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "i686.*|i386.*|x86.*|amd64.*|x86_64.*|AMD64.*")
    include(ofa/HandleX86Options)
    OFA_HandleX86Options()
  elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "^(arm.*|ARM.*|aarch64.*|AARCH64.*)")
    include(ofa/HandleArmOptions)
    OFA_HandleArmOptions()
  elseif("${CMAKE_SYSTEM_PROCESSOR}" MATCHES "^(powerpc|ppc)64.*")
    include(ofa/HandlePpcOptions)
    OFA_HandlePpcOptions()
  endif()
endmacro(OptimizeForArchitecture)
