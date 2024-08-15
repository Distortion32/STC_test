cmake_minimum_required(VERSION 3.5)

project(STC_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(STC_test main.cpp
    header.h)

include(GNUInstallDirs)
install(TARGETS STC_test
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
