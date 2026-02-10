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

    if(GISMO_AGENTS_FORCE_DOWNLOAD)
        message(STATUS "gsAgent Fetcher: Force-download enabled; will download/overwrite agent files from remote")
    else()
        message(STATUS "gsAgent Fetcher: Using local gsAgent clone when available")
    endif()

    # --- Get AGENT.md file (generic for all agents) ---
    if (GISMO_AGENTS_FORCE_DOWNLOAD OR NOT EXISTS "${CMAKE_SOURCE_DIR}/AGENTS.md")
        set(_base "https://raw.githubusercontent.com/gismo/gsAgents2/refs/heads/main/AGENTS.md")
        file(DOWNLOAD "${_base}" "${CMAKE_SOURCE_DIR}/AGENTS.md" STATUS _st)
        list(GET _st 0 _code)
        if(NOT _code EQUAL 0)
            message(STATUS "  No AGENTS.md found or download failed, skipping AGENTS.md.")
        else()
            message(STATUS "  Downloaded AGENTS.md")
        endif()
    else()
        message(STATUS "  AGENTS.md already exists, skipping download.")
    endif()

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


            set(_base "https://raw.githubusercontent.com/gismo/gsAgents2/refs/heads/output-${_prov}/${_subdirname}")
            file(DOWNLOAD "${_base}/manifest.txt"
                    "${CMAKE_BINARY_DIR}/${_prov}/${_type}/manifest.txt" STATUS _st)
            list(GET _st 0 _code)
            if(NOT _code EQUAL 0)
                # message(STATUS "  No manifest for ${_type}/${_prov}, skipping.")
                continue()
            endif()

            file(READ "${CMAKE_BINARY_DIR}/${_prov}/${_type}/manifest.txt" _raw)
            string(STRIP "${_raw}" _raw)
            separate_arguments(_items NATIVE_COMMAND "${_raw}")

            foreach(_item ${_items})
                if(_type STREQUAL "skills")
                    # Skills manifest can list either a skill name ("my-skill") which implies
                    # "my-skill/SKILL.md" or explicit file paths within the skill dir
                    # ("my-skill/examples.md", "my-skill/scripts/run.sh").
                    string(FIND "${_item}" "/" _has_slash)
                    if(_has_slash EQUAL -1)
                        set(_filename "${_item}/SKILL.md")
                    else()
                        set(_filename "${_item}")
                    endif()

                    set(_dest "${CMAKE_SOURCE_DIR}/${_subdirname}/${_filename}")

                    # Skip if the file exists (unless forcing)
                    if(EXISTS "${_dest}" AND NOT GISMO_AGENTS_FORCE_DOWNLOAD)
                        message(STATUS "-> ${_dest} exists, skipping")
                        continue()
                    endif()

                    get_filename_component(_destdir "${_dest}" DIRECTORY)
                    if(_destdir AND NOT IS_DIRECTORY "${_destdir}")
                        file(MAKE_DIRECTORY "${_destdir}")
                    endif()

                    file(DOWNLOAD "${_base}/${_filename}" "${_dest}" STATUS _dl)
                    list(GET _dl 0 _dlc)
                    if(_dlc EQUAL 0)
                        if(GISMO_AGENTS_FORCE_DOWNLOAD)
                            message(STATUS "-> [${_prov}/${_item}]: Downloaded/overwritten ${_filename}")
                        else()
                            message(STATUS "-> [${_prov}/${_item}]: Downloaded ${_filename}")
                        endif()
                    else()
                        message(STATUS "-> [${_prov}/${_item}]: Failed to download ${_filename} (code ${_dlc})")
                    endif()
                else()
                    set(_filename "${_item}.md")
                    set(_dest "${CMAKE_SOURCE_DIR}/${_subdirname}/${_filename}")

                    # Skip if the agent path exists (unless force-download is enabled)
                    if(EXISTS "${_dest}" AND NOT GISMO_AGENTS_FORCE_DOWNLOAD)
                        message(STATUS "-> ${_dest} exists, skipping")
                        continue()
                    endif()

                    get_filename_component(_destdir "${_dest}" DIRECTORY)
                    if(_destdir AND NOT IS_DIRECTORY "${_destdir}")
                        file(MAKE_DIRECTORY "${_destdir}")
                    endif()

                    file(DOWNLOAD "${_base}/${_item}.md" "${_dest}" STATUS _dl)
                    list(GET _dl 0 _dlc)
                    if(_dlc EQUAL 0)
                        if(GISMO_AGENTS_FORCE_DOWNLOAD)
                            message(STATUS "-> [${_prov}/${_item}]: Downloaded/overwritten ${_filename}")
                        else()
                            message(STATUS "-> [${_prov}/${_item}]: Downloaded ${_filename}")
                        endif()
                    else()
                        message(STATUS "-> [${_prov}/${_item}]: Failed to download ${_filename} (code ${_dlc})")
                    endif()
                endif()
            endforeach()
        endforeach()
    endforeach()
endfunction()
