cmake_minimum_required(VERSION 3.16)

project(universe_sandbox_clone C OBJC)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Emit compile_commands.json for IDE tooling.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenGL REQUIRED)

# SDL2
find_package(SDL2 QUIET)
if(NOT SDL2_FOUND)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(SDL2 REQUIRED sdl2)
endif()

# GLEW
find_package(GLEW QUIET)
if(NOT GLEW_FOUND)
  if(NOT PKG_CONFIG_FOUND)
    find_package(PkgConfig REQUIRED)
  endif()
  pkg_check_modules(GLEW REQUIRED glew)
endif()

add_executable(universe
  src/main.c
  src/util.c
  src/ui_text.c
  src/ui_overlay.c
  src/ui_labels.c
  src/input_mouse.c
  src/mat4.c
  src/camera.c
  src/scene_fbo.c
  src/render_bg.c
  src/render_lines.c
  src/render_particles.c
  src/render_rings.c
  src/presets.c
  src/sim.c
  src/glutil.c
  src/textures.c
  src/metal_accel.m
)

target_include_directories(universe PRIVATE third_party)

target_compile_options(universe PRIVATE -Wall -Wextra)

if(SDL2_FOUND)
  if(TARGET SDL2::SDL2)
    target_link_libraries(universe PRIVATE SDL2::SDL2)
  else()
    target_include_directories(universe PRIVATE ${SDL2_INCLUDE_DIRS})
    target_link_libraries(universe PRIVATE ${SDL2_LIBRARIES})
  endif()
else()
  target_include_directories(universe PRIVATE ${SDL2_INCLUDE_DIRS})
  target_link_libraries(universe PRIVATE ${SDL2_LIBRARIES})
endif()

if(GLEW_FOUND)
  if(TARGET GLEW::GLEW)
    target_link_libraries(universe PRIVATE GLEW::GLEW)
  else()
    target_include_directories(universe PRIVATE ${GLEW_INCLUDE_DIRS})
    target_link_libraries(universe PRIVATE ${GLEW_LIBRARIES})
  endif()
else()
  target_include_directories(universe PRIVATE ${GLEW_INCLUDE_DIRS})
  target_link_libraries(universe PRIVATE ${GLEW_LIBRARIES})
endif()

target_link_libraries(universe PRIVATE OpenGL::GL)

if(APPLE)
  target_link_libraries(universe PRIVATE "-framework Foundation" "-framework Metal")
  set_source_files_properties(src/metal_accel.m PROPERTIES COMPILE_FLAGS "-fobjc-arc")
endif()
