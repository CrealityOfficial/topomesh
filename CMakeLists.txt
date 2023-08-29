__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)
__cc_find(VCG)
__cc_find(cgal)
__cc_find(cxlog)

__files_group(${CMAKE_CURRENT_SOURCE_DIR}/topomesh/ SRCS)
		 
set(ILIBS trimesh2 cxnd cxutil eigen)
set(LIBS vcg boost_system mmesh rapidjson clipper)

set(INCS)
if(CGAL_INCLUDE_DIRS)
	list(APPEND INCS ${CGAL_INCLUDE_DIRS})
endif()

set(TOPOMESH_STATIC 1)
__add_common_library(topomesh OPENMP)