cmake_minimum_required(VERSION 3.05)
project(Coverage)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

    add_executable(Coverage main.cpp utilityfunction.h ParSSP.h RepSampling.h
            RoundRobinSieve.h RandomMultiGreedy.h)
