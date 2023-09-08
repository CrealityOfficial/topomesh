__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)

__assert_target(trimesh2)
__assert_target(eigen)
__assert_target(clipper)

__files_group(${CMAKE_CURRENT_SOURCE_DIR}/topomesh/ SRCS)
__files_group(${CMAKE_CURRENT_SOURCE_DIR}/internal/ Inters)
		 
list(APPEND SRCS ${Inters})
set(ILIBS trimesh2)
set(LIBS clipper eigen)

set(TOPOMESH_STATIC 1)
__add_common_library(topomesh OPENMP)