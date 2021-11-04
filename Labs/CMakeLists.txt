cmake_minimum_required (VERSION 3.8)

find_package (MPI REQUIRED)

file (GLOB TASK_DIRECTORIES "*")

foreach (DIR ${TASK_DIRECTORIES})
    if (IS_DIRECTORY ${DIR})
        get_filename_component (CURRENT_TARGET_NAME ${DIR} NAME)

        file (GLOB CURRENT_TARGET_SOURCE
                "${CURRENT_TARGET_NAME}/*.h"
                "${CURRENT_TARGET_NAME}/*.c"
                "${CURRENT_TARGET_NAME}/*.cpp")

        add_executable (${CURRENT_TARGET_NAME} ${CURRENT_TARGET_SOURCE})
        target_include_directories (${CURRENT_TARGET_NAME} PRIVATE ${MPI_CXX_INCLUDE_DIRS})
        target_link_libraries (${CURRENT_TARGET_NAME} PRIVATE ${MPI_CXX_LIBRARIES})

        install (TARGETS ${CURRENT_TARGET_NAME}
                RUNTIME DESTINATION "${CMAKE_SOURCE_DIR}/bin/"
                )
    endif ()
endforeach ()
