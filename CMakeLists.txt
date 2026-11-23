cmake_minimum_required(VERSION 3.16)
project(Mousquito LANGUAGES CXX VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Link the MSVC runtime statically so the exe runs without the VC++ redistributable.
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(MOUSQUITO_WERROR "Treat compiler warnings as errors (used by CI)" OFF)

# Platform-independent core: the movement engine. No OS calls in here.
add_library(mousquito_core STATIC
  src/trajectory.cpp
  src/version.cpp
)
target_include_directories(mousquito_core PUBLIC src)
target_compile_definitions(mousquito_core PRIVATE MQ_VERSION="${PROJECT_VERSION}")

if(MSVC)
  target_compile_options(mousquito_core PUBLIC /W4 /utf-8)
  if(MOUSQUITO_WERROR)
    target_compile_options(mousquito_core PUBLIC /WX)
  endif()
else()
  target_compile_options(mousquito_core PUBLIC -Wall -Wextra)
  if(MOUSQUITO_WERROR)
    target_compile_options(mousquito_core PUBLIC -Werror)
  endif()
endif()

add_executable(mousquito src/main.cpp)
target_link_libraries(mousquito PRIVATE mousquito_core)

enable_testing()
add_executable(mousquito_tests tests/core_tests.cpp)
target_link_libraries(mousquito_tests PRIVATE mousquito_core)
add_test(NAME core COMMAND mousquito_tests)
