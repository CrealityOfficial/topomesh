__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)

set(SRCS topomesh/data/mmesht.h
		 topomesh/data/mmesht.cpp
		 topomesh/alg/letter.h
		 topomesh/alg/letter.cpp
		 )
		 
set(LIBS trimesh2 eigen clipper)
__add_real_target(topomesh lib SOURCE ${SRCS}
							   LIB ${LIBS}
							   INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							   )