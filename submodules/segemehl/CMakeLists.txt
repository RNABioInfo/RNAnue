cmake_minimum_required(VERSION 3.10)
project(segemehl)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -DSORTEDUNMAPPED -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGNFO -DSHOWALIGN -DDBGLEVEL=0 -DPROGNFO -Wno-unused-parameter")
# HTSlib include paths come from find_package and target_include_directories
# below. Shell substitutions in compiler flags produce invalid Ninja files.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm -lpthread -lz -lncurses -lform -lmenu")

include_directories(include)

file(GLOB sources "libs/*.c")

find_package(HTSlib)
add_library(segemehl STATIC ${sources})

target_link_libraries(segemehl PRIVATE ${HTSlib_LIBRARIES})
target_include_directories(segemehl PUBLIC include ${HTSlib_INCLUDE_DIRS})
