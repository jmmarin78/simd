
cmake_minimum_required (VERSION 3.6)

set(FRAMEWORK_NAME "simd")                                # <== Set to your framework's name

include_directories("../pre/include")

set(USER_FRAMEWORK_OWNER_PATH "${CMAKE_CURRENT_LIST_DIR}")
include ("../../cmake/build_framework.cmake.txt")
