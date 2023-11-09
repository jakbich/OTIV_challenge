cmake_minimum_required(VERSION 3.5)
project(cpp_challenge)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

add_executable(main src/main.cpp src/railmap_helpers.cpp)
