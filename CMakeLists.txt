cmake_minimum_required(VERSION 3.16)

project(OriginalPinballCompilation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED IMPORTED_TARGET sdl2)
pkg_check_modules(SDL2IMAGE REQUIRED IMPORTED_TARGET SDL2_image)
pkg_check_modules(SDL2MIXER REQUIRED IMPORTED_TARGET SDL2_mixer)

add_executable(pinball engine.cpp)
target_link_libraries(pinball PRIVATE PkgConfig::SDL2 PkgConfig::SDL2IMAGE PkgConfig::SDL2MIXER)

