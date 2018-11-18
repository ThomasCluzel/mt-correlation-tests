cmake_minimum_required(VERSION 3.0)

# Le projet est étiqueté hello
project("mt-correlation-tests"
        VERSION 0.1
        LANGUAGES C)

# Include the library
add_subdirectory(libmt)

file(GLOB SRC *.c)

# add warning flags to the compiler
set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")

# Compile the executable => prog
add_executable(prog ${SRC})

# Linking of the MT library
target_link_libraries(prog mt19937ar)
