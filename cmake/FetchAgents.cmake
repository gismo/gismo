function(fetch_gismo_agents)
    set(GITHUB_REPO "gismo/gsAgents")

    foreach(PROVIDER ${GISMO_AGENT_PROVIDERS})
        # Use the /refs/heads/ path as you found
        set(BRANCH "agents-${PROVIDER}")
        set(BASE_URL "https://raw.githubusercontent.com/${GITHUB_REPO}/refs/heads/${BRANCH}")
        
        message(STATUS "G+Smo Agent Fetcher: Checking ${PROVIDER}...")

        # 1. Download manifest
        file(DOWNLOAD "${BASE_URL}/manifest.txt" "${CMAKE_BINARY_DIR}/manifest_${PROVIDER}.txt" STATUS DOWNLOAD_STATUS)
        
        list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
        if(NOT STATUS_CODE EQUAL 0)
            message(STATUS "  -> No manifest found for ${PROVIDER}. Skipping.")
            continue()
        endif()

        file(READ "${CMAKE_BINARY_DIR}/manifest_${PROVIDER}.txt" RAW_LIST)
        string(STRIP "${RAW_LIST}" RAW_LIST)
        separate_arguments(AGENT_LIST NATIVE_COMMAND "${RAW_LIST}")

        foreach(AGENT ${AGENT_LIST})
            # 2. Lookup the filename (Removed QUIET to fix the error)
            file(DOWNLOAD "${BASE_URL}/${AGENT}/.filename" "${CMAKE_BINARY_DIR}/${AGENT}_${PROVIDER}.name")
            
            if(EXISTS "${CMAKE_BINARY_DIR}/${AGENT}_${PROVIDER}.name")
                file(READ "${CMAKE_BINARY_DIR}/${AGENT}_${PROVIDER}.name" REAL_FILENAME)
                string(STRIP "${REAL_FILENAME}" REAL_FILENAME)
                
                # 3. Download the actual agent file
                message(STATUS "  -> Agent [${AGENT}]: Downloading ${REAL_FILENAME}")
                file(DOWNLOAD "${BASE_URL}/${AGENT}/${REAL_FILENAME}" 
                              "${CMAKE_SOURCE_DIR}/.agents/${AGENT}/${REAL_FILENAME}")
            endif()
        endforeach()
    endforeach()
endfunction()