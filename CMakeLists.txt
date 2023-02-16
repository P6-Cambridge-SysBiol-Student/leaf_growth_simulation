# Basic GLFW application, FJN 13.11.2021
set(CMAKE_CXX_STANDARD 11)

project(iLEAF)
cmake_minimum_required(VERSION 3.21)

set(TARGET iLEAF)

include_directories("deps")
find_package(OpenGL REQUIRED)

set(GLAD_GL "deps/glad/gl.h")

add_executable(${TARGET} WIN32 MACOSX_BUNDLE main.cc ${ICON} ${GLAD_GL})

target_link_libraries(${TARGET} "${PROJECT_SOURCE_DIR}/deps/libglfw3.a")
target_link_libraries(${TARGET} OpenGL::GL)

set_target_properties(${TARGET} PROPERTIES C_STANDARD 99)

if (APPLE)
    set_target_properties(${TARGET} PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME "iLEAF")
    target_link_libraries(${TARGET}
    "-framework Cocoa"
    "-framework OpenGL"
    "-framework IOKit"
    )
    set_source_files_properties(glfw.icns PROPERTIES
                                MACOSX_PACKAGE_LOCATION "Resources")
    set_target_properties(${GUI_ONLY_BINARIES} PROPERTIES
                          MACOSX_BUNDLE_SHORT_VERSION_STRING ${GLFW_VERSION}
                          MACOSX_BUNDLE_LONG_VERSION_STRING ${GLFW_VERSION}
                          MACOSX_BUNDLE_ICON_FILE glfw.icns
                          MACOSX_BUNDLE_INFO_PLIST "${GLFW_SOURCE_DIR}/CMake/Info.plist.in")
endif()
