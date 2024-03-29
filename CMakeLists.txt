project(LEAFSIM)
cmake_minimum_required(VERSION 3.21)

set(TARGET leafsim)

include_directories("deps")
find_package(OpenGL REQUIRED)

set(GLAD_GL "deps/glad/gl.h" createTriangles.h polish.h vector.h hormone.h arrays.h sigmoid.h graphics.h springs.h writing.h fitness.h)

add_executable(${TARGET} WIN32 MACOSX_BUNDLE main.cc ${ICON} ${GLAD_GL})

target_link_libraries(${TARGET} "${PROJECT_SOURCE_DIR}/deps/libglfw3.a")
target_link_libraries(${TARGET} OpenGL::GL)

set_target_properties(${TARGET} PROPERTIES C_STANDARD 99)

if (APPLE)
    set(ICON deps/glfw.icns)
    set_target_properties(${TARGET} PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME "Leafsim")
    target_link_libraries(${TARGET}
    "-framework Cocoa"
    "-framework OpenGL"
    "-framework IOKit"
    )
else()
    set(ICON deps/glfw.rc)
endif()
