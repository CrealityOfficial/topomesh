#include "mmesht.h"

namespace topomesh
{
	MMeshT::MMeshT(trimesh::TriMesh* currentMesh)
	{
		if (currentMesh->vertices.size() < 3000)	this->vertices.reserve(3000);
		if (currentMesh->faces.size() < 1024)	this->faces.reserve(1024);
		int vn = 0;
		for (trimesh::point apoint : currentMesh->vertices)
		{
			this->vertices.push_back(apoint);
			this->vertices.back().index = vn;
			if (currentMesh->normals.size() > 0)
				this->vertices.back().normal = currentMesh->normals.at(vn);
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
		for (int i = 0; i < this->vertices.size(); i++)
		{
			for (int j = 0; j < this->vertices[i].connected_vertex.size(); j++)
				this->vertices[i].connected_vertex[j] = &this->vertices[this->vertices[i].connected_vertex[j]->index];
			for (int j = 0; j < this->vertices[i].connected_face.size(); j++)
				this->vertices[i].connected_face[j] = &this->faces[this->vertices[i].connected_face[j]->index];
		}
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

	void MMeshT::calculateCrossPoint(std::vector<trimesh::ivec2>& edge, std::vector<trimesh::point>& line, std::vector<trimesh::vec3>& tc)
	{
		for (trimesh::ivec2& e : edge)
		{
			int size = line.size();
			for (int i = 0; i < line.size(); i++)
			{
				trimesh::point b1 = trimesh::point(this->vertices[e.x].p.x, this->vertices[e.x].p.y, 0);
				trimesh::point b2 = trimesh::point(this->vertices[e.y].p.x, this->vertices[e.y].p.y, 0);
				trimesh::point a = line[i] - line[(i + 1) % size];//a2->a1
				trimesh::point b = b1 - b2;//b2->b1
				trimesh::point n1 = line[i] - b1;//b1->a1
				trimesh::point n2 = line[i] - b2;//b2->a1
				trimesh::point n3 = line[(i + 1) % size] - b1;//b1->a2
				if (((-a % -n2) ^ (-a % -n1)) >= 0 || ((-b % n1) ^ (-b % n3)) >= 0) continue;
				trimesh::point m = a % b;
				trimesh::point n = n1 % a;
				float t = (m.x != 0 && n.x != 0) ? n.x / m.x : (m.y != 0 && n.y != 0) ? n.y / m.y : (m.z != 0 && n.z != 0) ? n.z / m.z : -1;
				if (t > 0 && t < 1)//b1-t*b			
					tc.push_back(trimesh::vec3(t, e.x, e.y));
			}

		}
	}

	void MMeshT::deleteVertex(MMeshVertex& v)
	{
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
			for (int i = 0; i < f.connect_vertex.size(); i++)
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
			}
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
		this->vertices.back().index = this->vn;
		this->vn++;
	}

	void MMeshT::appendVertex(trimesh::point& v)
	{
		this->vertices.push_back(v);
		this->vertices.back().index = this->vn;
		this->vn++;
	}

	void MMeshT::appendFace(MMeshVertex& v0, MMeshVertex& v1, MMeshVertex& v2)
	{
		bool invert = false;
		this->faces.push_back(MMeshFace(&this->vertices[v0.index], &this->vertices[v1.index], &this->vertices[v2.index]));
		this->faces.back().index = this->fn;

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
}