cmake_minimum_required (VERSION 3.0)
project (FTCHEB)

if (APPLE)

    include_directories(/opt/local/include)
    find_library(macblas
        NAMES ACCELERATE
        PATHS ${CMAKE_OSX_SYSROOT}/System/Library
        PATH_SUFFIXES Frameworks
        NO_DEFAULT_PATH)

    set(blas_lapack ${macblas}/Accelerate)
    set(compiler_flags "-std=c99 -O2 -Wall -Weverything -Wno-shorten-64-to-32 -Wno-padded -Wno-missing-prototypes -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-missing-noreturn -Wno-reserved-id-macro -Wno-bad-function-cast -Wno-double-promotion")

    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /usr/local/lib)
else ()
    set(compiler_flags "-g -O2 -ftree-vectorize -Wall -Wextra -std=c99 -fno-omit-frame-pointer")

    set(CMAKE_C_COMPILER "gcc")

    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    include_directories(/usr/include)
    set(blas_lapack "-lblas -llapack")

endif()

find_library(c3lib c3 PATHS “/usr/local/lib”)
include_directories(/usr/local/include)
#set_target_properties(c3 PROPERTIES #    IMPORTED_LOCATION "/Users/msoley/Downloads/c3/build/src/libc3.dylib" #    INTERFACE_INCLUDE_DIRECTORIES "/Users/msoley/Downloads/c3/build/include")

include_directories(${BLAS_INCLUDE_DIR}, ${LAPACK_INCLUDE_DIR})

set(CMAKE_C_FLAGS "${compiler_flags}")

add_executable(ft_cheb 2dpro.c auxiliary.c)
target_link_libraries(ft_cheb
    -g 
    ${blas_lapack} 
    -lm
    ${c3lib} )


