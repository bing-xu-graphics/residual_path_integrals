cmake_minimum_required(VERSION 3.20)
project(residual_path_integrals LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(Threads REQUIRED)
set(EMBREE_ROOT "" CACHE PATH "Embree 3 installation prefix")
find_path(EMBREE_INCLUDE_DIR embree3/rtcore.h
    HINTS "${EMBREE_ROOT}" PATH_SUFFIXES include)
find_library(EMBREE_LIBRARY NAMES embree3
    HINTS "${EMBREE_ROOT}"
    PATH_SUFFIXES lib lib64 lib-win32 lib-linux lib-macos)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Embree REQUIRED_VARS EMBREE_INCLUDE_DIR EMBREE_LIBRARY)

add_library(glacier_lib STATIC
    src/3rdparty/miniz.c src/3rdparty/pugixml.cpp
    src/parsers/load_serialized.cpp src/parsers/parse_obj.cpp
    src/parsers/parse_ply.cpp src/parsers/parse_scene.cpp
    src/camera.cpp src/filter.cpp src/image.cpp src/intersection.cpp
    src/light.cpp src/material.cpp src/medium.cpp src/parallel.cpp
    src/phase_function.cpp src/render.cpp src/scene.cpp src/shape.cpp
    src/table_dist.cpp src/transform.cpp src/volume.cpp)
target_include_directories(glacier_lib PUBLIC src "${EMBREE_INCLUDE_DIR}")
target_link_libraries(glacier_lib PUBLIC "${EMBREE_LIBRARY}" Threads::Threads)
if(MSVC)
    target_compile_definitions(glacier_lib PUBLIC _USE_MATH_DEFINES)
    target_compile_options(glacier_lib PUBLIC $<$<COMPILE_LANGUAGE:CXX>:/EHsc>)
endif()
add_executable(glacier src/main.cpp)
target_link_libraries(glacier PRIVATE glacier_lib)
add_executable(image_compare tools/image_compare.cpp)
target_link_libraries(image_compare PRIVATE glacier_lib)
if(WIN32)
    get_filename_component(embree_library_dir "${EMBREE_LIBRARY}" DIRECTORY)
    set(EMBREE_RUNTIME_DIR "${embree_library_dir}/../bin" CACHE PATH "Embree and TBB DLL directory")
    file(GLOB embree_runtime_dlls "${EMBREE_RUNTIME_DIR}/*.dll")
    foreach(target glacier image_compare)
        if(embree_runtime_dlls)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${embree_runtime_dlls} "$<TARGET_FILE_DIR:${target}>")
        endif()
    endforeach()
endif()

include(CTest)
if(BUILD_TESTING)
    foreach(test_name filter frame image intersection materials matrix mipmap)
        add_executable(test_${test_name} src/tests/${test_name}.cpp)
        target_link_libraries(test_${test_name} PRIVATE glacier_lib)
        add_test(NAME ${test_name} COMMAND test_${test_name})
        set_tests_properties(${test_name} PROPERTIES PASS_REGULAR_EXPRESSION "SUCCESS")
    endforeach()
endif()
