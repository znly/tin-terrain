cmake_minimum_required(VERSION 3.0)

project(tin-terrain)

option(TNTN_TEST "include test targets in the buildsystem" OFF)
option(TNTN_DOWNLOAD_DEPS "download dependencies during cmake configure" ON)
#option(TNTN_USE_ADDONS "" OFF)

set(CMAKE_CXX_STANDARD 14)


# get the current working branch
execute_process(
  COMMAND git describe --all --dirty
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE TNTN_GIT_DESCRIPTION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

# last commit hash
execute_process(
  COMMAND git log -1 --format=%H
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE TNTN_GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )
string(TIMESTAMP TNTN_TIMESTAMP "%Y%m%dT%H%MZ" UTC)

find_package(GDAL)
if(NOT GDAL_FOUND)
  message(STATUS "GDAL not found, cannot build command line tool")
else()
  if(NOT GDAL_CONFIG)
    message(FATAL_ERROR "gdal-config command not found (not in PATH?), cannot proceed")
  endif()

  execute_process(
    COMMAND ${GDAL_CONFIG} --version
    OUTPUT_VARIABLE SYSTEM_GDAL_VERSION
    )

  if(SYSTEM_GDAL_VERSION VERSION_LESS "2.2")
    message(FATAL_ERROR "GDAL version \"${SYSTEM_GDAL_VERSION}\" is too old, at least 2.2 is required")
  endif()
  set(TNTN_HAVE_GDAL 1)
endif()

# Boost
# set(Boost_USE_STATIC_LIBS ON CACHE BOOL "" FORCE)
# set(Boost_USE_STATIC_LIBS ON)
# find_package(Boost REQUIRED COMPONENTS
#     program_options
#     filesystem
#     system
# )
# if(TNTN_DOWNLOAD_DEPS)
#   include("download-deps.cmake")
#   find_path(TNTN_LIBGLM_SOURCE_DIR NAMES "glm/glm.hpp" HINTS "${PROJECT_SOURCE_DIR}/3rdparty/glm-0.9.9.0/")
#   find_path(TNTN_LIBFMT_SOURCE_DIR NAMES "include/fmt/format.h" HINTS "${PROJECT_SOURCE_DIR}/3rdparty/fmt-5.1.0/")
#   find_path(TNTN_CATCH2_SOURCE_DIR NAMES "catch.hpp" HINTS "${PROJECT_SOURCE_DIR}/3rdparty/Catch2-2.3.0/")
#   if(NOT TNTN_LIBGLM_SOURCE_DIR)
#     message(SEND_ERROR "GLM math library source dir not found, please download and set TNTN_LIBGLM_SOURCE_DIR")
#   endif()
#   if(NOT TNTN_LIBFMT_SOURCE_DIR)
#     message(SEND_ERROR "libfmt string formatting library source dir not found, please download and set TNTN_LIBFMT_SOURCE_DIR")
#   endif()
#   #add_subdirectory("${TNTN_LIBFMT_SOURCE_DIR}")
# endif()

set(TNTN_COMMON_SOURCE_FILES
  include/tntn/File.h
  include/tntn/logging.h
  include/tntn/util.h
  src/File.cpp
  src/logging.cpp
  src/util.cpp
  )

set(TNTN_MESHING_SOURCE_FILES
  include/tntn/DelaunayMesh.h
  include/tntn/DelaunayTriangle.h
  include/tntn/Mesh.h
  include/tntn/Mesh2Raster.h
  include/tntn/MeshIO.h
  include/tntn/MeshMode.h
  include/tntn/OFFReader.h
  include/tntn/QuadEdge.h
  include/tntn/SuperTriangle.h
  include/tntn/SurfacePoints.h
  include/tntn/TerraMesh.h
  include/tntn/TerraUtils.h
  include/tntn/ZemlyaMesh.h
  include/tntn/geometrix.h
  include/tntn/raster_tools.h
  include/tntn/simple_meshing.h
  src/DelaunayMesh.cpp
  src/DelaunayTriangle.cpp
  src/Delaunator.cpp
  src/Mesh.cpp
  src/Mesh2Raster.cpp
  src/MeshIO.cpp
  src/OFFReader.cpp
  src/QuadEdge.cpp
  src/SuperTriangle.cpp
  src/SurfacePoints.cpp
  src/TerraMesh.cpp
  src/TerraUtils.cpp
  src/ZemlyaMesh.cpp
  src/geometrix.cpp
  src/raster_tools.cpp
  src/simple_meshing.cpp
  )

set(TNTN_SOURCE_FILES
  "${CMAKE_CURRENT_BINARY_DIR}/version_info.cpp"
  include/delaunator_cpp/Delaunator.h
  include/tntn/BinaryIO.h
  include/tntn/FileFormat.h
  include/tntn/MercatorProjection.h
  include/tntn/MeshWriter.h
  include/tntn/ObjPool.h
  include/tntn/Points2Mesh.h
  include/tntn/QuantizedMeshIO.h
  include/tntn/Raster.h
  include/tntn/RasterIO.h
  include/tntn/RasterOverviews.h
  include/tntn/TileMaker.h
  include/tntn/ZoomRange.h
  include/tntn/benchmark_workflow.h
  include/tntn/dem2tintiles_workflow.h
  include/tntn/endianness.h
  include/tntn/gdal_init.h
  include/tntn/println.h
  include/tntn/terra_meshing.h
  include/tntn/tntn_assert.h
  include/tntn/version_info.h
  include/tntn/zemlya_meshing.h
  src/BinaryIO.cpp
  src/Delaunator.cpp
  src/MercatorProjection.cpp
  src/MeshWriter.cpp
  src/Points2Mesh.cpp
  src/QuantizedMeshIO.cpp
  src/RasterIO.cpp
  src/RasterOverviews.cpp
  src/TileMaker.cpp
  src/benchmark_workflow.cpp
  src/dem2tintiles_workflow.cpp
  src/gdal_init.cpp
  src/println.cpp
  src/terra_meshing.cpp
  src/version_info.cpp.in
  src/zemlya_meshing.cpp
  )

configure_file(src/version_info.cpp.in version_info.cpp)

if(TNTN_USE_ADDONS)
  list(APPEND TNTN_SOURCE_FILES
    include/tntn/Raster2Mesh.h
    src/Raster2Mesh.cpp
    )
endif()

add_library(tntn_meshing STATIC
  ${TNTN_MESHING_SOURCE_FILES}
  ${TNTN_COMMON_SOURCE_FILES}
  )

target_include_directories(tntn_meshing
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include/
  )

target_link_libraries(tntn_meshing
  PUBLIC
  Boost::boost
  Boost::filesystem
  fmt::fmt
  glm
  )

if (TNTN_HAVE_GDAL)

  add_library(tntn STATIC
    ${TNTN_SOURCE_FILES}
    ${TNTN_MESHING_SOURCE_FILES}
    ${TNTN_COMMON_SOURCE_FILES}
    )

  target_compile_definitions(tntn PUBLIC $<$<CONFIG:Debug>:TNTN_DEBUG>)

  target_compile_definitions(tntn PUBLIC
    TNTN_HAVE_GDAL)

  #target_compile_definitions(tntn PUBLIC GLM_FORCE_SWIZZLE GLM_ENABLE_EXPERIMENTAL)

  if(TNTN_DOWNLOAD_ADDONS)
    target_compile_definitions(tntn PUBLIC TNTN_USE_ADDONS=1)
  endif()


  target_include_directories(tntn
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/
    PRIVATE
    ${GDAL_INCLUDE_DIR}
    )
  if(TNTN_USE_ADDONS)
    target_include_directories(tntn
      PUBLIC
      3rdparty/addons/include/
      )
  endif()

  target_link_libraries(tntn
    PUBLIC
    Boost::boost
    Boost::filesystem
    fmt::fmt
    glm
    PRIVATE
    ${GDAL_LIBRARY}
    )

  add_executable(tin-terrain
    src/cmd.cpp
    )
  target_link_libraries(tin-terrain
    PRIVATE
    Boost::program_options
    tntn
    )

  if(TNTN_TEST)
    add_subdirectory(test)
  endif()

endif()
