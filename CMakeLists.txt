__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)
__cc_find(cxlog)

__files_group(${CMAKE_CURRENT_SOURCE_DIR}/topomesh/ SRCS)
		 
set(ILIBS trimesh2 cxnd cxutil mmesh eigen)
set(LIBS boost_system mmesh clipper)

set(TOPOMESH_STATIC 1)
__add_common_library(topomesh OPENMP)