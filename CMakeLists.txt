cmake_minimum_required(VERSION 3.0)

set(DEBUG YES)

# Le projet est étiqueté hello
project("mt-correlation-tests"
        VERSION 0.1
        LANGUAGES C)

# add warning flags to the compiler
if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Debug or release
if(${DEBUG})
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
