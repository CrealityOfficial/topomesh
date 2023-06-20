__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)
__cc_find(VCG)
__cc_find(cgal)
__cc_find(cxlog)

set(SRCS topomesh/data/mmesht.h
		 topomesh/data/mmesht.cpp
		 topomesh/data/mmeshFace.h
		 topomesh/data/mmeshFace.cpp
		 topomesh/data/mmeshVertex.h
		 topomesh/data/mmeshVertex.cpp		
		 topomesh/alg/letter.h
		 topomesh/alg/letter.cpp
		 topomesh/alg/convert.h
		 topomesh/alg/convert.cpp
		 topomesh/alg/cutting.h
		 topomesh/alg/cutting.cpp
		 topomesh/alg/fillhoneycombs.h
		 topomesh/alg/fillhoneycombs.cpp
		 )
		 
set(ILIBS trimesh2 eigen mmesh clipper)
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
							   )