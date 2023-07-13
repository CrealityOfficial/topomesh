__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)
__cc_find(VCG)
__cc_find(cgal)
__cc_find(cxlog)

__files_group(${CMAKE_CURRENT_SOURCE_DIR}/topomesh/ SRCS)
		 
set(ILIBS trimesh2 eigen mmesh clipper rapidjson)
set(LIBS vcg boost_system cxlog)

set(INCS)
if(CGAL_INCLUDE_DIRS)
	list(APPEND INCS ${CGAL_INCLUDE_DIRS})
endif()

__add_real_target(topomesh lib SOURCE ${SRCS}
							   ILIB ${ILIBS}
							   LIB ${LIBS}
							   INC ${INCS}
							   INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							   OPENMP
							   )