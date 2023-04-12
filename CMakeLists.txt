__cc_find(TriMesh2)
__cc_find(Eigen)
__cc_find(Clipper)

set(SRCS topomesh/data/mmesht.h
		 topomesh/data/mmesht.cpp
		 topomesh/data/mmeshFace.h
		 topomesh/data/mmeshFace.cpp
		 topomesh/data/mmeshVertex.h
		 topomesh/data/mmeshVertex.cpp
		 topomesh/data/convert.h
		 topomesh/data/convert.cpp
		 topomesh/alg/letter.h
		 topomesh/alg/letter.cpp
		 )
		 
set(ILIBS trimesh2 eigen mmesh clipper)
__add_real_target(topomesh lib SOURCE ${SRCS}
							   ILIB ${ILIBS}
							   INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}
							   )