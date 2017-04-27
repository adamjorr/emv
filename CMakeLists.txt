cmake_minimum_required(VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)
project (meep)
add_executable(meep meep.cc)
target_link_libraries(meep hts)

add_library(libmeep STATIC
  em.cc
  genotype.cc
  pileup.cc
  plpdata.cc
  reftype.cc
  samio.cc
  seqem.cc
)

target_link_libraries(meep libmeep)

