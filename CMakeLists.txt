cmake_minimum_required(VERSION 3.0)

# The name of the project
project(mt-correlation-tests
        VERSION 1.0
        LANGUAGES C CXX
        DESCRIPTION "Project to test Mersenne Twister with C. Ismay's tests")

# Build type (release by default for efficiency)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
endif()

# Include MT and Ismay's tests
add_subdirectory(externals/libmt)
add_subdirectory(externals/libmultivartest)

# Add warnings
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra)
endif()

# Add the unit tests
add_subdirectory(tests)

# Compile the executable => prog
add_executable(prog src/main.cpp
                    src/allexperiments.cpp src/allexperiments.hpp
                    src/experiment.cpp src/experiment.hpp)

# Linking of MT and Ismay's tests
target_link_libraries(prog PRIVATE mt19937ar)
target_link_libraries(prog PRIVATE multivartest)
