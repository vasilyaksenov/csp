cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("csp")

set (CMAKE_CXX_STANDARD 17)

# Build OR-tools dependencies.
set(BUILD_DEPS ON)

set(BUILD_SAMPLES OFF)
set(BUILD_EXAMPLES OFF)

# Download or-tools and it's dependecies
include(FetchContent)
FetchContent_Declare(
  or-tools
  GIT_REPOSITORY https://github.com/google/or-tools.git
  GIT_TAG        v9.6
)
FetchContent_MakeAvailable(or-tools)

FetchContent_Declare(Argumentum
   GIT_REPOSITORY https://github.com/mmahnic/argumentum.git
   GIT_TAG master
)
FetchContent_MakeAvailable(Argumentum)

add_executable (csp "csp.cpp" "csp.h")

target_link_libraries(csp ortools::ortools Argumentum::headers)
