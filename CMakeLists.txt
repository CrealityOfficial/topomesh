__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)
__cc_find(VCG)

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
		 )
		 
set(ILIBS trimesh2 eigen mmesh clipper)
set(LIBS vcg)

__add_real_target(topomesh lib SOURCE ${SRCS}
							   ILIB ${ILIBS}
							   LIB ${LIBS}
							   INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							   )