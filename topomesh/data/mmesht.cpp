#include "mmesht.h"

namespace topomesh
{
	MMeshT::MMeshT(trimesh::TriMesh* currentMesh)
	{		
		if(currentMesh->faces.size()<4096)
			this->faces.reserve(8192);
		else
			this->faces.reserve((unsigned)(currentMesh->faces.size() * 1.5));
		if (currentMesh->vertices.size() < 4096)
			this->vertices.reserve(8192);
		else
			this->vertices.reserve((unsigned)(currentMesh->vertices.size() * 1.5));		
		
		int vn = 0;
		if (currentMesh->normals.size() > 0)
			this->VertexNormals = true;
		for (trimesh::point apoint : currentMesh->vertices)
		{
			this->vertices.push_back(apoint);
			this->vertices.back().index = vn;
			if (this->is_VertexNormals())
			{
				this->vertices.back().normal = currentMesh->normals.at(vn);
			}
			vn++;
		}
		this->vn = vn;
		for (int i = 0; i < this->vn; i++)
		{
			if (currentMesh->neighbors.size() > 0)
			{
				this->VVadjacent = true;
				for (int j = 0; j < currentMesh->neighbors[i].size(); j++)
					this->vertices[i].connected_vertex.push_back(&this->vertices[currentMesh->neighbors[i][j]]);
			}
		}

		int fn = 0;
		for (trimesh::TriMesh::Face f : currentMesh->faces)
		{
			this->faces.push_back(f);
			this->faces.back().index = fn;
			for (int i = 0; i < 3; i++)
				this->faces.back().connect_vertex.push_back(&this->vertices[f[i]]);
			fn++;
		}
		this->fn = fn;
		if (currentMesh->adjacentfaces.size() > 0)
		{
			this->FFadjacent = true;
			this->VFadjacent = true;
			for (int i = 0; i < this->vn; i++)
			{
				for (int j = 0; j < currentMesh->adjacentfaces[i].size(); j++)
					this->vertices[i].connected_face.push_back(&this->faces[currentMesh->adjacentfaces[i][j]]);
			}
			for (int n = 0; n < this->faces.size(); n++)
			{
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < this->faces[n].connect_vertex[i]->connected_face.size(); j++)
						this->faces[n].connect_vertex[i]->connected_face[j]->SetA();

				for (int i = 0; i < 3; i++)
					for (int j = 0; j < this->faces[n].connect_vertex[i]->connected_face.size(); j++)
						if (this->faces[n].connect_vertex[i]->connected_face[j]->IsA(2) && !this->faces[n].connect_vertex[i]->connected_face[j]->IsS())
						{
							this->faces[n].connect_vertex[i]->connected_face[j]->SetS();
							this->faces[n].connect_face.push_back(this->faces[n].connect_vertex[i]->connected_face[j]);
						}

				for (int i = 0; i < 3; i++)
					for (int j = 0; j < this->faces[n].connect_vertex[i]->connected_face.size(); j++)
					{
						this->faces[n].connect_vertex[i]->connected_face[j]->ClearA();
						this->faces[n].connect_vertex[i]->connected_face[j]->ClearS();
					}
			}
		}
	}


	void MMeshT::mmesh2trimesh(trimesh::TriMesh* currentMesh)
	{
		shrinkMesh();
		currentMesh->clear();
		for (int i = 0; i < this->vertices.size(); i++)
		{
			if (this->vertices[i].IsD()) continue;
			currentMesh->vertices.push_back(this->vertices[i].p);
			if(this->is_VertexNormals())
				currentMesh->normals.push_back(this->vertices[i].normal);
			if (this->is_VVadjacent())
			{
				std::vector<int> vvadj;
				for (MMeshVertex* v : this->vertices[i].connected_vertex)
					vvadj.push_back(v->index);
				currentMesh->neighbors.push_back(vvadj);
			}
			if (this->is_VFadjacent())
			{
				std::vector<int> vfadj;
				for (MMeshFace* f : this->vertices[i].connected_face)
					vfadj.push_back(f->index);
				currentMesh->adjacentfaces.push_back(vfadj);
			}
		}
		for (int i = 0; i < this->faces.size(); i++)
		{
			if (this->faces[i].IsD()) continue;
			currentMesh->faces.push_back(trimesh::TriMesh::Face(this->faces[i].connect_vertex[0]->index,
				this->faces[i].connect_vertex[1]->index,
				this->faces[i].connect_vertex[2]->index));
			if (this->is_FFadjacent())
			{
				//.....
			}
		}
	}

	void MMeshT::shrinkMesh()
	{
		int deleteVNum = 0;
		int deleteFNum = 0;
		for (MMeshVertex& v : this->vertices)
		{
			if (v.IsD())
			{
				deleteVNum++; continue;
			}
			v.index -= deleteVNum;
		}
		for (MMeshFace& f : this->faces)
		{
			if (f.IsD())
			{
				deleteFNum++; continue;
			}
			f.index -= deleteFNum;
		}
		if(this->is_VVadjacent())
			for (int i = 0; i < this->vertices.size(); i++)
			{
				for (int j = 0; j < this->vertices[i].connected_vertex.size(); j++)
					this->vertices[i].connected_vertex[j] = &this->vertices[this->vertices[i].connected_vertex[j]->index];
				for (int j = 0; j < this->vertices[i].connected_face.size(); j++)
					this->vertices[i].connected_face[j] = &this->faces[this->vertices[i].connected_face[j]->index];
			}
		if(this->is_VFadjacent())
			for (int i = 0; i < this->faces.size(); i++)
			{
				for (int j = 0; j < this->faces[i].connect_vertex.size(); j++)
					this->faces[i].connect_vertex[j] = &this->vertices[this->faces[i].connect_vertex[j]->index];
				for (int j = 0; j < this->faces[i].connect_face.size(); j++)
					this->faces[i].connect_face[j] = &this->faces[this->faces[i].connect_face[j]->index];
			}
		for (int i = 0; i < this->vertices.size(); i++)
			if (this->vertices[i].IsD())
			{
				this->vertices.erase(this->vertices.begin() + i);
				i--;
			}
		for (int i = 0; i < this->faces.size(); i++)
			if (this->faces[i].IsD())
			{
				this->faces.erase(this->faces.begin() + i);
				i--;
			}
	}

	void MMeshT::getMeshBoundary()
	{
		if (!this->is_VVadjacent()) return;
		for (MMeshVertex& v : this->vertices)
		{
			for (int i = 0; i < v.connected_face.size(); i++)
			{
				v.connected_face[i]->V1(&v)->SetA();
				v.connected_face[i]->V2(&v)->SetA();
			}
			for (int i = 0; i < v.connected_vertex.size(); i++)
			{
				if (v.connected_vertex[i]->IsA(1))
				{
					v.connected_vertex[i]->SetB();
					v.connected_vertex[i]->ClearA();
				}
			}
		}
		for (MMeshVertex& v : this->vertices)
			v.ClearA();
	}

	void MMeshT::getEdge(std::vector<trimesh::ivec2>& edge, bool is_select)
	{
		for (MMeshVertex& v : this->vertices)
		{
			if (is_select && !v.IsS()) continue;
			for (int i = 0; i < v.connected_vertex.size(); i++)
			{
				if (v.connected_vertex[i]->IsV()) continue;
				edge.push_back(trimesh::ivec2(v.index, v.connected_vertex[i]->index));
			}
			v.SetV();
		}
		for (MMeshVertex& v : this->vertices)
		{
			v.ClearV();
		}
	}

	void MMeshT::calculateCrossPoint(std::vector<trimesh::ivec2>& edge, std::pair<trimesh::point, trimesh::point>& line, std::vector<std::pair<float, trimesh::ivec2>>& tc)
	{		
		for (trimesh::ivec2& e : edge)
		{
			trimesh::point b1 = trimesh::point(this->vertices[e.x].p.x, this->vertices[e.x].p.y, 0);
			trimesh::point b2 = trimesh::point(this->vertices[e.y].p.x, this->vertices[e.y].p.y, 0);			
			trimesh::point a = line.first - line.second;//a2->a1
			trimesh::point b = b1 - b2;//b2->b1
			trimesh::point n1 = line.first - b1;//b1->a1
			trimesh::point n2 = line.first - b2;//b2->a1
			trimesh::point n3 = line.second - b1;//b1->a2								
			if (((-a % -n2) ^ (-a % -n1)) >= 0 || ((-b % n1) ^ (-b % n3)) >= 0) continue;
			trimesh::point m = a % b;
			trimesh::point n = n1 % a;
			float t = (m.x != 0 && n.x != 0) ? n.x / m.x : (m.y != 0 && n.y != 0) ? n.y / m.y : (m.z != 0 && n.z != 0) ? n.z / m.z : -1;
			if (t > 0 && t < 1)//b1-t*b							
			{
				tc.push_back(std::make_pair(t, trimesh::ivec2(e.x, e.y)));				
			}
		}		
	}

	void MMeshT::deleteVertex(MMeshVertex& v)
	{
		if (v.IsD()) return;
		v.SetD();
		if (this->is_VFadjacent() && this->is_VVadjacent())
		{
			std::vector<MMeshFace*> deleteface = v.connected_face;
			for (int i = 0; i < deleteface.size(); i++)
				deleteFace(*deleteface[i]);
		}
		this->vn--;
	}

	void MMeshT::deleteVertex(int i)
	{
		deleteVertex(this->vertices[i]);
	}

	void MMeshT::deleteFace(MMeshFace& f)
	{
		if (f.IsD()) return;
		f.SetD();
		if (this->is_VFadjacent())
		{
			for (int i = 0; i < f.connect_vertex.size(); i++)
			{
				for (int j = 0; j < f.connect_vertex[i]->connected_face.size(); j++)
				{
					if (f == f.connect_vertex[i]->connected_face[j])
					{
						f.connect_vertex[i]->connected_face.erase(f.connect_vertex[i]->connected_face.begin() + j);
						break;
					}
				}
			}
		}
		if (this->is_VVadjacent())
		{
			for (MMeshVertex* v : f.connect_vertex)
				v->connected_vertex.clear();
			for (MMeshVertex* v : f.connect_vertex)
			{
				for (MMeshFace* vf : v->connected_face)
				{
					v->connected_vertex.push_back(vf->V1(v));
					v->connected_vertex.push_back(vf->V2(v));
				}
			}
			/*for (int i = 0; i < f.connect_vertex.size(); i++)
			{
				if (f.V0(i)->IsB() && f.V1(i)->IsB())
				{
					for (int j = 0; j < f.V0(i)->connected_vertex.size(); j++)
					{
						if (f.V1(i) == f.V0(i)->connected_vertex[j])
						{
							f.V0(i)->connected_vertex.erase(f.V0(i)->connected_vertex.begin() + j);
							break;
						}
					}
				}
				if (f.V0(i)->IsB() && f.V2(i)->IsB())
				{
					for (int j = 0; j < f.V0(i)->connected_vertex.size(); j++)
					{
						if (f.V2(i) == f.V0(i)->connected_vertex[j])
						{
							f.V0(i)->connected_vertex.erase(f.V0(i)->connected_vertex.begin() + j);
							break;
						}
					}
				}
			}*/
		}
		if (this->is_FFadjacent())
		{
			for (int i = 0; i < f.connect_face.size(); i++)
			{
				for (int j = 0; j < f.connect_face[i]->connect_face.size(); j++)
				{
					if (f == f.connect_face[i]->connect_face[j])
					{
						f.connect_face[i]->connect_face.erase(f.connect_face[i]->connect_face.begin() + j);
						break;
					}
				}
			}
		}
		this->fn--;
	}

	void MMeshT::deleteFace(int i)
	{
		deleteFace(this->faces[i]);
	}

	void MMeshT::appendVertex(MMeshVertex& v)
	{
		this->vertices.push_back(v);
		this->vertices.back().index = this->vertices.size()-1;
		this->vn++;
	}

	void MMeshT::appendVertex(trimesh::point& v)
	{
		this->vertices.push_back(v);
		this->vertices.back().index = this->vertices.size()-1;
		this->vn++;
	}

	void MMeshT::appendFace(MMeshVertex& v0, MMeshVertex& v1, MMeshVertex& v2)
	{				
		bool invert = false;		
		this->faces.push_back(MMeshFace(&this->vertices[v0.index], &this->vertices[v1.index], &this->vertices[v2.index]));		
		this->faces.back().index = this->faces.size()-1;	
		
		if (this->is_FFadjacent())
		{
			for (int i = 0; i < v0.connected_face.size(); i++)
				v0.connected_face[i]->SetA();
			for (int i = 0; i < v1.connected_face.size(); i++)
				v1.connected_face[i]->SetA();
			for (int i = 0; i < v2.connected_face.size(); i++)
				v2.connected_face[i]->SetA();
			//---可能重复添加connect_face
			for (int i = 0; i < v0.connected_face.size(); i++)
				if (v0.connected_face[i]->IsA(2) && !v0.connected_face[i]->IsL())
				{
					v0.connected_face[i]->SetL();
					this->faces.back().connect_face.push_back(v0.connected_face[i]);
					v0.connected_face[i]->connect_face.push_back(&this->faces.back());
					if (v0.connected_face[i]->V1(&v0) == &v1 || v0.connected_face[i]->V2(&v0) == &v2)
						invert = true;
				}
			for (int i = 0; i < v1.connected_face.size(); i++)
				if (v1.connected_face[i]->IsA(2) && !v1.connected_face[i]->IsL())
				{
					v1.connected_face[i]->SetL();
					this->faces.back().connect_face.push_back(v1.connected_face[i]);
					v1.connected_face[i]->connect_face.push_back(&this->faces.back());
					if (v1.connected_face[i]->V1(&v1) == &v2 && v1.connected_face[i]->V2(&v1) == &v0)
						invert = true;
					
				}
			for (int i = 0; i < v2.connected_face.size(); i++)
				if (v2.connected_face[i]->IsA(2) && !v2.connected_face[i]->IsL())
				{
					v2.connected_face[i]->SetL();
					this->faces.back().connect_face.push_back(v2.connected_face[i]);
					v2.connected_face[i]->connect_face.push_back(&this->faces.back());
					if (v2.connected_face[i]->V1(&v2) != &v0 && v2.connected_face[i]->V2(&v2) != &v1)
						invert = true;
				}
			for (int i = 0; i < v0.connected_face.size(); i++) {
				v0.connected_face[i]->ClearL(); v0.connected_face[i]->ClearA();
			}
			for (int i = 0; i < v1.connected_face.size(); i++) {
				v1.connected_face[i]->ClearL(); v1.connected_face[i]->ClearA();
			}
			for (int i = 0; i < v2.connected_face.size(); i++) {
				v2.connected_face[i]->ClearL(); v2.connected_face[i]->ClearA();
			}
		}
		
		if (this->is_VFadjacent())
		{
			v0.connected_face.push_back(&this->faces.back());
			v1.connected_face.push_back(&this->faces.back());
			v2.connected_face.push_back(&this->faces.back());
		}
		
		if (this->is_VVadjacent())
		{
			bool edge1 = false, edge2 = false, edge3 = false;
			for (int i = 0; i < v0.connected_vertex.size(); i++)
			{
				if (v0.connected_vertex[i] == &v1)
					edge1 = true;
				if (v0.connected_vertex[i] == &v2)
					edge2 = true;
			}
			for (int i = 0; i < v1.connected_vertex.size(); i++)
				if (v1.connected_vertex[i] == &v2)
					edge3 = true;
			if (!edge1)
			{
				v0.connected_vertex.push_back(&v1); v1.connected_vertex.push_back(&v0);
			}
			if (!edge2)
			{
				v0.connected_vertex.push_back(&v2); v2.connected_vertex.push_back(&v0);
			}
			if (!edge3)
			{
				v1.connected_vertex.push_back(&v2); v2.connected_vertex.push_back(&v1);
			}
		}
		if (invert)
			std::swap(this->faces.back().connect_vertex[1], this->faces.back().connect_vertex[2]);
		
		this->fn++;
	}

	void MMeshT::appendFace(int i0, int i1, int i2)
	{
		appendFace(this->vertices[i0], this->vertices[i1], this->vertices[i2]);
	}

	double  getTotalArea(std::vector<trimesh::point>& inVertices)
	{
		double sum = 0.0f;
		for (int n = 1; n < inVertices.size() - 1; n++)
		{
			sum += MMeshT::det(inVertices[0], inVertices[n], inVertices[n + 1]);
		}
		return abs(sum);
	}

	void MMeshT::getFacesNormals()
	{
		for (MMeshFace& f : this->faces)
		{
			trimesh::point ave_normal;
			ave_normal += f.V0(0)->normal;
			ave_normal += f.V0(1)->normal;
			ave_normal += f.V0(2)->normal;
			ave_normal /= 3.0f;
			f.normal = ave_normal;
		}
		this->FacesNormals = true;
	}
}