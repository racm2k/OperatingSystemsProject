cmake_minimum_required(VERSION 3.17)
project(Projeto C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(Projeto
        main_sockets_server.c)
