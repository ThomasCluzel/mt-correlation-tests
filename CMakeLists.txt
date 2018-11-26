cmake_minimum_required(VERSION 3.0)

set(DEBUG "yes")

# Le projet est étiqueté hello
project("mt-correlation-tests"
        VERSION 0.1
        LANGUAGES C)

# add warning flags to the compiler
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
if(${DEBUG} EQUAL "yes")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# Include the library
add_subdirectory(libmt)

# Add the tests
add_subdirectory(tests)

# Source files
file(GLOB SRC *.c)

# Compile the executable => prog
add_executable(prog ${SRC})

# Linking of the MT library
target_link_libraries(prog mt19937ar)
