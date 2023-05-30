#include "cutting.h"

namespace topomesh
{
	bool ModleCutting(const std::vector<trimesh::TriMesh*>& inMesh, std::vector<trimesh::TriMesh*>& outMesh, const SimpleCamera& camera,
		const TriPolygon& paths, ccglobal::Tracer* tracer)
	{
		CameraParam cp;
		cp.lookAt = camera.center;
		cp.pos = camera.pos;
		cp.up = camera.up;
		cp.n = camera.n; cp.f = camera.f;
		cp.fov = camera.fov; cp.aspect = camera.aspect;
		loadCameraParam(cp);
		Eigen::Matrix4f viewMatrix;
		Eigen::Matrix4f projectionMatrix;
		getViewMatrixAndProjectionMatrix(cp, viewMatrix, projectionMatrix);
		std::vector<MMeshT> meshts;
		for (int i = 0; i < inMesh.size(); i++)
		{
			inMesh[i]->need_adjacentfaces();
			inMesh[i]->need_neighbors();
			TransformationMesh(inMesh[i], viewMatrix, projectionMatrix);
			if (JudgeMeshIsVaild(inMesh[i]))
			{				
				meshts.push_back(inMesh[i]);
			}
		}
		std::vector<std::vector<trimesh::vec2>> lines(1);
		for (int i = 0; i < paths.size(); i++)
		{
			lines[0].push_back(trimesh::vec2(paths[i].x, paths[i].y));
		}
		std::vector<MMeshT> newMeshs;
		for (int i = 0; i < meshts.size(); i++)
		{
			std::vector<int> faceindex;
			for (int fi = 0; fi < meshts[i].faces.size(); fi++)
			{
				if (!meshts[i].faces[fi].IsL())
					faceindex.emplace_back(fi);
			}
			int beforeVertexSize = meshts[i].vertices.size();
			embedingAndCutting(&meshts[i], lines,faceindex);			
			for (int vi = beforeVertexSize; vi < meshts[i].vertices.size(); vi++)
			{
				meshts[i].vertices[vi].SetB();
				/*for (int vf = 0; vf < meshts[i].vertices[vi].connected_face.size(); vf++)
					meshts[i].vertices[vi].connected_face[vf]->SetB();*/
			}
			splitMesh(&meshts[i], newMeshs);
			trimesh::TriMesh* trimesh1 = new trimesh::TriMesh();
			trimesh::TriMesh* trimesh2 = new trimesh::TriMesh();			
			newMeshs[0].mmesh2trimesh(trimesh1);
			newMeshs[1].mmesh2trimesh(trimesh2);			
			unTransformationMesh(trimesh1, viewMatrix, projectionMatrix);
			unTransformationMesh(trimesh2, viewMatrix, projectionMatrix);			
			trimesh1->write("trimesh1.ply");
			trimesh2->write("trimesh2.ply");	

			/*trimesh::TriMesh* trimesh1 = new trimesh::TriMesh();
			meshts[i].mmesh2trimesh(trimesh1);
			unTransformationMesh(trimesh1, viewMatrix, projectionMatrix);
			trimesh1->write("trimesh1.ply");*/
		}
		if (meshts.empty()) return false;

		return true;
	}

	//------nlogn---
	bool JudgeCloseOfPath(const TriPolygon& paths)
	{
		for (int i = 0; i < paths.size() - 1; i++)
		{
			trimesh::point l1b = trimesh::point(paths[i].x, paths[i].y, 0);
			trimesh::point l1e = trimesh::point(paths[i + 1].x, paths[i + 1].y, 0);
			trimesh::point l1 = l1e - l1b;
			for (int j = i + 1; j < paths.size() - 1; j++)
			{
				trimesh::point l2b = trimesh::point(paths[j].x, paths[j].y, 0);
				trimesh::point l2e = trimesh::point(paths[j+1].x, paths[j+1].y, 0);
				trimesh::point l2 = l2e - l2b;
				trimesh::point n1 = l1b - l2b;
				trimesh::point n2 = l1b - l2e;
				trimesh::point n3 = l1e - l2b;								
				if (((l1 % -n1) ^ (l1 % -n2)) >= 0 || ((l2 % n1) ^ (l2 % n3)) >= 0) continue;
				trimesh::point m = l1 % l2;
				trimesh::point n = n1 % l1;
				float t = (m.x != 0 && n.x != 0) ? n.x / m.x : (m.y != 0 && n.y != 0) ? n.y / m.y : (m.z != 0 && n.z != 0) ? n.z / m.z : -1;
				if (t > 0 && t < 1)						
				{
					return true;
				}
			}
		}
		return false;
	}

	bool JudgeMeshIsVaild(const trimesh::TriMesh* inMesh)
	{
		for (int i = 0; i < inMesh->vertices.size(); i++)
		{
			if (inMesh->vertices[i].x > -1 && inMesh->vertices[i].x < 1 && inMesh->vertices[i].y < 1 && inMesh->vertices[i].y>-1)
				return true;
		}
		return false;
	}


	void setMarkOfCorssPoint(std::vector<MMeshT*>& meshs, const TriPolygon& paths, const std::vector<std::pair<int, int>>& corssPoint)
	{
		for (int i = 0; i < meshs.size(); i++)
		{
			for (int fi = 0; fi < meshs[i]->faces.size(); fi++)
			{
				for (int j = 0; j < corssPoint.size(); j++)
				{
					if (meshs[i]->faces[fi].innerFace(paths[corssPoint[j].first]) || meshs[i]->faces[fi].innerFace(paths[corssPoint[j].second]))
					{
						meshs[i]->faces[fi].SetL(); break;
					}
				}
			}
		}
	}
	

	void splitMesh(MMeshT* mesh, std::vector<MMeshT>& outmesh)
	{				
		for (int i = 0; i < mesh->faces.size(); i++)
			mesh->faces[i].ClearS();
		std::queue<int> queue;
		int fn = 0; int it = 0;
		while (fn <mesh->FN())
		{
			int index = -1;
			for (int i = 0; i < mesh->faces.size(); i++)
				if (!mesh->faces[i].IsS()&&!mesh->faces[i].IsD())
				{
					index = i; break;
				}	
			if (index == -1)
				break;
			queue.push(index);		
			mesh->faces[index].SetS();
			std::vector<int> faceindex;
			faceindex.push_back(index);
			while (!queue.empty())
			{				
				std::vector<int> boundpointindex;
				for (int i = 0; i < mesh->faces[queue.front()].connect_vertex.size(); i++)
				{
					if (mesh->faces[queue.front()].connect_vertex[i]->IsB())
						boundpointindex.push_back(mesh->faces[queue.front()].connect_vertex[i]->index);
				}	
				std::sort(boundpointindex.begin(), boundpointindex.end());
				std::vector<std::pair<int, int>> seq;
				for (int i = 1; i < boundpointindex.size(); i++)
				{
					if (boundpointindex[i] - boundpointindex[i - 1] == 1)
						seq.push_back(std::make_pair(boundpointindex[i], boundpointindex[i - 1]));
				}			
				
				for (int i = 0; i < mesh->faces[queue.front()].connect_face.size(); i++)
				{
					bool pass = true;
					for (int k = 0; k < seq.size(); k++)
					{
						int n = 0;
						for (int j = 0; j < mesh->faces[queue.front()].connect_face[i]->connect_vertex.size(); j++)
						{
							if (mesh->faces[queue.front()].connect_face[i]->connect_vertex[j]->index == seq[k].first ||
								mesh->faces[queue.front()].connect_face[i]->connect_vertex[j]->index == seq[k].second)
								n++;
						}
						if (n == 2)
						{
							pass = false; break;
						}
					}
					if (!mesh->faces[queue.front()].connect_face[i]->IsS()&&pass)
					{							
						queue.push(mesh->faces[queue.front()].connect_face[i]->index);
						mesh->faces[queue.front()].connect_face[i]->SetS();
						faceindex.push_back(mesh->faces[queue.front()].connect_face[i]->index); fn++;
					}
				}
				
				queue.pop();
			}		
			
			/*for (int i = 0; i < mesh->vertices.size(); i++)
			{
				if(mesh->vertices[i].IsV()&&!mesh->vertices[i].IsB())
					for (int j = 0; j < mesh->vertices[i].connected_face.size(); j++)
					{
						if (!mesh->vertices[i].connected_face[j]->IsS())
						{
							mesh->vertices[i].connected_face[j]->SetS(); faceindex.push_back(mesh->vertices[i].connected_face[j]->index); fn++;
						}
					}				
			}*/
			MMeshT mt(mesh, faceindex);
			/*trimesh::TriMesh* trimesh1 = new trimesh::TriMesh();
			mt.mmesh2trimesh(trimesh1);
			trimesh1->write("trimesh1.ply");*/
			auto&& mtrr = std::move(mt);			
			outmesh.push_back(std::move(mtrr));	
			it++;
		}
		if (it == 1)
			outmesh.clear();
		
	}
}