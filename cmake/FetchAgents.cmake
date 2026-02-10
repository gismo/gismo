######################################################################
## FetchAgents.cmake
## This file is part of the gsAgent compiler.
##
## Creates symlinks from the project root to agent instruction files
## in output/<provider>/agents/. Falls back to downloading from GitHub.
##
## Supported providers: claude, gemini, opencode, cursor, copilot,
##                      windsurf, cline
######################################################################

# Provider link registry: PROVIDER|ROOT_PATH|AGENTS_PATH
# ROOT_PATH  is relative to CMAKE_SOURCE_DIR (the project root).
# AGENTS_PATH is relative to output/.

# set(_AGENT_LINKS
#     "claude|.claude/commands|claude/commands"
#     "claude|.claude/agents|claude/agents"
#     "claude|.claude/rules|claude/rules"
#     "gemini|GEMINI.md|gemini/agents/GEMINI.md"
#     "opencode|.opencode/agents|opencode/agents"
#     "copilot|.github/agents|github/agents"
# )

set (_AGENT_LINKS
    "claude/agents|.claude/agents"
    "claude/commands|.claude/commands"
    "claude/skills|.claude/skills"
    "claude/rules|.claude/rules"
    "opencode/agents|.opencode/agents"
    "opencode/commands|.opencode/commands"
    "opencode/skills|.opencode/skills"
    "opencode/rules|.opencode/rules"
    "copilot/agents|.github/agents"
    "copilot/commands|.github/commands"
    "copilot/skills|.github/skills"
)

function(_get_agent_dir_name PROVIDER OUT_VAR)
    foreach(_entry ${_AGENT_LINKS})
        string(REPLACE "|" ";" _parts "${_entry}")
        list(GET _parts 0 _prov)
        list(GET _parts 1 _agents_path)
        if(_prov STREQUAL "${PROVIDER}")
            set(${OUT_VAR} "${_agents_path}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
    set(${OUT_VAR} "" PARENT_SCOPE)
endfunction()

## When the script is run via `cmake -P`, invoke the fetcher.
#if(NOT DEFINED CMAKE_RUNNING_FROM_PARENT)
#    # Expect GISMO_AGENT_PROVIDERS to be set by the caller (CMake cache or -D)
#    if(NOT GISMO_AGENT_PROVIDERS)
#        set(GISMO_AGENT_PROVIDERS claude gemini opencode)
#    endif()
#    fetch_gismo_agents()
#endif()

function(fetch_gismo_agents)
    set(_dir "${CMAKE_SOURCE_DIR}")

    message(STATUS "gsAgent Fetcher: Using local gsAgent clone")

    # --- Remote fallback: download from GitHub ---
    message(STATUS "gsAgent Fetcher: Downloading from GitHub")

    if(NOT GISMO_AGENT_TYPES)
        set(GISMO_AGENT_TYPES agents commands rules skills)
    endif()

    foreach(_type ${GISMO_AGENT_TYPES})
        foreach(_prov ${GISMO_AGENT_PROVIDERS})

            _get_agent_dir_name("${_prov}/${_type}" _subdirname)
            if(NOT _subdirname)
                continue()
            endif()

            # message(STATUS "Processing provider ${_prov} (registry path: ${_subdirname})")

            set(_base "https://raw.githubusercontent.com/gismo/gsAgents2/refs/heads/agents-${_prov}/${_subdirname}")
            file(DOWNLOAD "${_base}/manifest.txt"
                    "${CMAKE_BINARY_DIR}/${_type}/manifest_${_prov}.txt" STATUS _st)
            list(GET _st 0 _code)
            if(NOT _code EQUAL 0)
                # message(STATUS "  No manifest for ${_type}/${_prov}, skipping.")
                continue()
            endif()

            file(READ "${CMAKE_BINARY_DIR}/${_type}/manifest_${_prov}.txt" _raw)
            string(STRIP "${_raw}" _raw)
            separate_arguments(_items NATIVE_COMMAND "${_raw}")

            foreach(_item ${_items})
                set(_filename "${_item}.md")
                set(_dest "${CMAKE_SOURCE_DIR}/${_subdirname}/${_filename}")

                # Skip if the agent path exists
                if(EXISTS "${_dest}")
                    message(STATUS "-> ${_dest} exists, skipping")
                    continue()
                endif()

                get_filename_component(_destdir "${_dest}" DIRECTORY)
                if(_destdir AND NOT IS_DIRECTORY "${_destdir}")
                    file(MAKE_DIRECTORY "${_destdir}")
                endif()

                file(DOWNLOAD "${_base}/${_item}.md" "${_dest}" STATUS _dl)
                message("${_base}/${_item}.md")
                list(GET _dl 0 _dlc)
                message(STATUS "-> [${_prov}/${_item}]: Downloaded ${_filename}")
            endforeach()
        endforeach()
    endforeach()
endfunction()
