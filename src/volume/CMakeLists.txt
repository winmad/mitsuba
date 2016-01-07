# Volumetric data sources
if (NOT MTS_VERSION)
  message(FATAL_ERROR "Use the top level configuration file")
endif()

macro(add_volume)
  add_mts_plugin(${ARGN} TYPE volume)
endmacro()

add_volume(shellmap2 shellmap2.cpp)
add_volume(instancedvol instancedvol.cpp)
add_volume(tiledvol    tiledvol.cpp)
add_volume(gridvolumeSimple gridvol_simple.cpp)

add_volume(constvolume constvolume.cpp)
add_volume(gridvolume  gridvolume.cpp)
add_volume(hgridvolume hgridvolume.cpp)
add_volume(volcache    volcache.cpp)
