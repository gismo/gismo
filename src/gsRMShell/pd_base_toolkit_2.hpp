#pragma  once
//
//#include <vector>
//#include <set>
//#include <math.h>
//#include <map>
//#include <algorithm>
//#include <assert.h>
//#include <iostream>
//#include <fstream>
//#include <iterator>
////#include <xfunctional>
//#include "pd_base_toolkit.hpp"
//
//using namespace std;
//
//using namespace gismo;
//
//using namespace Eigen;
//
//namespace DLUT
//{
//	namespace SAE
//	{
//
//		namespace PERIDYNAMIC
//		{
//			//	计算四边形面积
//			template<typename Type>
//			double Calculate_Area_Square(const Type& p1, const Type& p2,
//				const Type& p3, const Type& p4)
//			{
//				double A1 = Calculate_Area_Triangle(p1, p2, p3);
//				double A2 = Calculate_Area_Triangle(p1, p3, p4);
//				return (A1 + A2);
//			}
//			//	计算四边形形心
//			template<typename Type>
//			Type Calculate_ShapeCenter_Square(const Type& p1, const Type& p2,
//				const Type& p3, const Type& p4)
//			{
//				Type center, center1, center2;
//				center1 = Calculate_ShapeCenter_Triangle(p1, p2, p3);
//				center2 = Calculate_ShapeCenter_Triangle(p1, p3, p4);
//
//				double area1, area2;
//				area1 = Calculate_Area_Triangle(p1, p2, p3);
//				area2 = Calculate_Area_Triangle(p1, p3, p4);
//
//				center.x() = (area1 * center1.x() + area2 * center2.x()) / (area1 + area2);
//				center.y() = (area1 * center1.y() + area2 * center2.y()) / (area1 + area2);
//				center.z() = (area1 * center1.z() + area2 * center2.z()) / (area1 + area2);
//
//				return center;
//			}
//			//  FA
//			template<typename Type>
//			double Calculate_Intersect_Area_FA(const Type& c1, double r1, const Type& c2, double r2)
//			{
//				double r = Min(r1, r2);
//				double R = Max(r1, r2);
//				double d = Distance_2pt(c1, c2);
//				if (d > R + 1E-3)
//				{
//					return 0;
//				}
//				else
//				{
//					return PI * r * r;
//				}
//
//			}
//			//  LAMMPS
//			template<typename Type>
//			double Calculate_Intersect_Area_LAMMPS(const Type& c1, double r1, const Type& c2, double r2, double meshsize)
//			{
//				double r = Min(r1, r2);
//				double R = Max(r1, r2);
//				double half_of_meshsize = 0.5 * meshsize;
//				double d = Distance_2pt(c1, c2);
//				if (d > R + 1E-3)
//				{
//					return 0;
//				}
//				else if ((d + half_of_meshsize) <= R + 1E-3)
//				{
//					return PI * r * r;
//				}
//				else
//				{
//					return (R + half_of_meshsize - d) / meshsize * PI * r * r;
//				}
//
//			}
//		}
//
//
//
//		namespace PERIDYNAMIC
//		{
//			//	Node
//			class TFemNode
//			{
//			public:
//				TFemNode() : m_p_local_id(NULL)
//				{
//					m_set_adj_element_ids.clear();
//					m_volume_of_node = 0.0;
//				}
//				TFemNode(const Vector3d& pt, int nid = -1)
//				{
//					m_p_local_id = NULL;
//					m_volume_of_node = 0.0;
//
//					m_tp_point = pt;
//					m_nid_global = nid;
//					m_set_adj_element_ids.clear();
//				}
//				void				Dispose()
//				{
//					if (m_p_local_id)
//					{
//						delete m_p_local_id;
//						m_p_local_id = NULL;
//					}
//					m_set_adj_element_ids.clear();
//				}
//			public:
//				Vector3d& Point() { return m_tp_point; }
//				Vector3d			Point() const { return m_tp_point; }
//
//				int& Id()
//				{
//					if (m_p_local_id == NULL) m_p_local_id = new int(-1);
//					return *m_p_local_id;
//				}
//				int					Id() const
//				{
//					if (m_p_local_id == NULL)
//						return -1;
//					else
//						return *m_p_local_id;
//				}
//
//				//int&				Id() { return m_nid_global; }
//				//int				Id() const { return m_nid_global; }
//				int& IdGlobal() { return m_nid_global; }
//				int					IdGlobal() const { return m_nid_global; }
//
//				double& VolumeOfNode() { return m_volume_of_node; }
//				double				VolumeOfNode() const { return m_volume_of_node; }
//
//				void				InsertAdjElement(int* eid) { m_set_adj_element_ids.insert(eid); }
//				void				DeleteAdjElement(int* eid) { m_set_adj_element_ids.erase(eid); }
//
//				set<int>			GetAdjElementIds() const
//				{
//					set<int> res;
//					for (int* eid : m_set_adj_element_ids)
//					{
//						res.insert(*eid);
//					}
//					return res;
//				}
//
//				vector<Vector3d>& PointsOfVolume() { return m_PointsVolume; }
//				vector<Vector3d>    PointsOfVolume() const { return m_PointsVolume; }
//
//			public:
//				Vector3d			m_tp_point;							//	Coordinate
//				set<int*>			m_set_adj_element_ids;				//	Adjoint element set
//				vector<Vector3d>    m_PointsVolume;                     //  Points to shape the volume
//				double				m_volume_of_node;					//	Volume of this node
//				int* m_p_local_id;						//	Id of this node
//				int					m_nid_global;						//	Global Id of this node
//			};
//
//			//	Element
//			class TFemElement
//			{
//			public:
//				TFemElement() : m_p_local_id(NULL) {}
//				void				Dispose()
//				{
//					if (m_p_local_id)
//					{
//						delete m_p_local_id;
//						m_p_local_id = NULL;
//					}
//				}
//			public:
//				void				InsertNode(int* nid) { m_vec_nids.push_back(nid); }
//				vector<int>			GetNodeIds() const
//				{
//					vector<int> res;
//					for (int* pnid : m_vec_nids)
//					{
//						res.push_back(*pnid);
//					}
//					return res;
//				}
//
//				int& Id()
//				{
//					if (m_p_local_id == NULL)
//						m_p_local_id = new int(-1);
//
//					return *m_p_local_id;
//				}
//				int					Id() const
//				{
//					if (m_p_local_id == NULL)
//						return -1;
//					else
//						return *m_p_local_id;
//				}
//				int& IdGlobal() { return m_eid_global; }
//				int					IdGlobal() const { return m_eid_global; }
//
//				int& PartId() { return m_part_id; }
//				int					PartId() const { return m_part_id; }
//
//				double& VolumeOfElement() { return m_volume_of_element; }
//				double				VolumeOfElement() const { return m_volume_of_element; }
//
//				int& IdPd() { return m_nid_pd; }
//				int					IdPd() const { return m_nid_pd; }
//			public:
//				vector<int*>		m_vec_nids;						//	Node Ids of this element
//				int					m_part_id;						//	Part Id of this element
//				int* m_p_local_id;					//	Id of this element
//				double				m_volume_of_element;			//	Volume of this element
//				int					m_eid_global;					//	Global Id of this element
//				int					m_nid_pd;						//	Id in PD model
//				Vector3d            m_shapeCenter;                  //  Shape center
//			};
//
//			//	FemMesh
//			class TFemMeshCoreData
//			{
//			public:
//				TFemMeshCoreData()
//				{}
//				~TFemMeshCoreData()
//				{
//					for (vector<TFemNode>::iterator iterNode = m_nodes.begin();
//						iterNode != m_nodes.end(); ++iterNode)
//					{
//						iterNode->Dispose();
//					}
//					m_nodes.clear();
//
//					for (vector<TFemElement>::iterator iterElem = m_elements.begin();
//						iterElem != m_elements.end(); ++iterElem)
//					{
//						iterElem->Dispose();
//					}
//					m_elements.clear();
//
//					m_vec_output_node_ids.clear();
//					m_vec_output_element_ids.clear();
//				}
//			public:
//				void Initialize()
//				{
//					m_nodes.clear();
//					m_elements.clear();
//				}
//			public:
//				void InsertNode(const Vector3d& pt, int nid)
//				{
//					TFemNode node(pt, nid);
//					//m_nodes.push_back(TFemNode(pt, nid));
//					m_nodes.push_back(node);
//					m_nodes.back().Id() = (int)(m_nodes.size() - 1);
//				}
//				int					NodeCount() const { return (int)(m_nodes.size()); }
//				TFemNode& FemNodeById(int nid_local)
//				{
//					return m_nodes[nid_local];
//				}
//				set<int>			GetFemNodeIds()
//				{
//					set<int> nids;
//
//					for (const TFemElement& element : m_elements)
//					{
//						const vector<int>& ele_nids = element.GetNodeIds();
//						for (int nid : ele_nids)
//						{
//							nids.insert(nid);
//						}
//					}
//
//					return nids;
//				}
//
//				void InsertElement(const vector<int>& nids, int id_global, int part_id)
//				{
//					m_elements.push_back(TFemElement());// TFemElement()是一个类，下面的Id(),IdGlobal()...都是其成员函数
//					TFemElement& element = m_elements.back(); // element 指向了 m_element 的末尾
//					element.Id() = (int)(m_elements.size() - 1);
//					element.IdGlobal() = id_global;
//					element.PartId() = part_id;
//					for (int nid : nids)
//					{
//						element.InsertNode(&(m_nodes[nid].Id()));
//
//						TFemNode& node = FemNodeById(nid); // return m_nodes[nid]
//						node.InsertAdjElement(&(element.Id()));
//					}
//				}
//				int	ElementCount() const { return (int)(m_elements.size()); }
//
//				TFemElement& FemElementById(int eid_local)
//				{
//					return m_elements[eid_local];
//				}
//				const TFemElement& FemElementById(int eid_local) const
//				{
//					return m_elements[eid_local];
//				}
//
//				set<int> GetFemElementIds()
//				{
//					set<int> eids;
//					for (const TFemElement& element : m_elements)
//					{
//						eids.insert(element.Id());
//					}
//					return eids;
//				}
//
//			public:
//				void RefreshInfo()
//				{
//					for (vector<TFemElement>::iterator iter = m_elements.begin();
//						iter != m_elements.end(); ++iter)
//					{
//						TFemElement& element = *iter;
//						const vector<int>& nids = element.GetNodeIds();
//						int nCount = (int)(nids.size());
//						switch (nCount)
//						{
//						case 4:
//						{
//							const Vector3d& n0 = m_nodes[nids[0]].Point();
//							const Vector3d& n1 = m_nodes[nids[1]].Point();
//							const Vector3d& n2 = m_nodes[nids[2]].Point();
//							const Vector3d& n3 = m_nodes[nids[3]].Point();
//							double a = Calculate_Area_Triangle(n0, n1, n2);
//							double b = Calculate_Area_Triangle(n2, n3, n0);
//							element.VolumeOfElement() = a + b;
//							break;
//						}
//						default:
//							break;
//						}
//					}
//
//					for (vector<TFemNode>::iterator iter = m_nodes.begin();
//						iter != m_nodes.end(); ++iter)
//					{
//
//					}
//				}
//			public:
//				set<int>			GetPdNodeIdsByPart(int pid)
//				{
//					set<int> res;
//					res.clear();
//					for (const TFemElement& element : m_elements)
//					{
//						if (element.PartId() == pid)
//						{
//							res.insert(element.IdPd());
//						}
//					}
//
//					return res;
//				}
//				void				AddOutputNode(int nid_local)
//				{
//					m_vec_output_node_ids.insert(nid_local);
//				}
//				void				AddOutputElement(int eid_local)
//				{
//					m_vec_output_element_ids.insert(eid_local);
//				}
//				const set<int>& GetOutputNodeIds() const { return m_vec_output_node_ids; }
//				const set<int>& GetOutputElementIds() const { return m_vec_output_element_ids; }
//			public:
//				vector<TFemNode>		m_nodes;
//				vector<TFemElement>		m_elements;
//			public:
//				set<int>				m_vec_output_node_ids;
//				set<int>				m_vec_output_element_ids;
//			};
//
//			//计算一个点到一个网格的最小距离
//			template<typename Type>
//			inline void Distance_Node2Mesh(Type& pt,
//				TFemMeshCoreData& mesh,
//				double& mini_dis,
//				int& id)
//			{
//				set<int> nodeIds = mesh.GetFemNodeIds();
//				mini_dis = 1E10;
//				id = -1;
//				for (set<int>::const_iterator it = nodeIds.begin(); it != nodeIds.end(); it++)
//				{
//					double disT = Distance_2pt(pt, (mesh.FemNodeById(*it)).Point());
//					if (mini_dis > disT)
//					{
//						mini_dis = disT;
//						id = *it;
//					}
//				}
//				return;
//			}
//
//		}
//
//	}
//}
