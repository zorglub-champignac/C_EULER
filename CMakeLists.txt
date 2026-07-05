cmake_minimum_required(VERSION 4.2)
project(C_EULER C)

set(CMAKE_C_STANDARD 11)

add_executable(fix_final fix_final.c)
add_executable(solve-iter solve-iter.c)
add_executable(solve solve.c)
