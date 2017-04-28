cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)
project (meep)

add_library(libmeep STATIC
  genotype.cc
  pileup.cc
  plpdata.cc
  reftype.cc
  samio.cc
  seqem.cc
)

target_link_libraries(libmeep hts)

add_executable(meep meep.cc)
target_link_libraries(meep hts)
target_link_libraries(meep libmeep)

