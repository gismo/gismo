/*All rights reserved

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
--Please append file description informations here --

== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
Date            Name                    Description of Change

$HISTORY$
== == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
*/

#ifndef DLUT_SAE_PERIDYNAMIC_PD_DATABASE_HXX_20171220
#define DLUT_SAE_PERIDYNAMIC_PD_DATABASE_HXX_20171220

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include <map>
#include <string>
//#include "Eigen/Dense"
//#include "Eigen/IterativeLinearSolvers"
//#include "Eigen/SparseLU"
//#include "Eigen/Eigenvalues"
//#include "Eigen"
#include "umfSolver.h"
#include "pd_base_toolkit.hpp"
#include "fem_database.hpp"
#include <ppl.h>

using namespace Eigen;
using namespace std;
using namespace concurrency;

extern double RATIO_OF_HORIZON_MESHSIZE;
extern double HORIZON;

namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			
			//static double RATIO_OF_HORIZON_MESHSIZE = 3.0;
			//static double HORIZON = 3.0;
			static bool USE_CONSTANT_HORIZON = false;

			/************************************************************************/
			/* ����Beam PDʱ���ڵ����6�����ɶ�,����BBPDʱ���ڵ����3�����ɶ�       */
			/************************************************************************/
			typedef Eigen::MatrixXd SingleStiffness;

			class TPdCalculateParas
			{
			public:
				TPdCalculateParas()
				{
					cax = 0;
					cby = 0;
					cbz = 0;
					ctor = 0;
					csy = 0;
					density = 0;
					b_facture = true;

					c = 0;
					s0 = 0;
				}
				~TPdCalculateParas() {}
			public:
				double	cax;
				double	cby;
				double	cbz;
				double	ctor;
				double	csy;
				double	density;
				bool	b_facture;

				double	c;
				double	s0;
			};

			//	PD Bond
			class TPdBond
			{
			public:
				TPdBond(double volume, TCoordinate center)
				{
					m_volume = volume;
					m_center = center;
				}
				TPdBond(double volume = 0)
				{
					m_volume = volume;
					m_center.setZero();
				}
				~TPdBond() {}
			public:
				Eigen::VectorXd& ForceOfBond() { return m_force; }
				const Eigen::VectorXd& ForceOfBond() const { return m_force; }

				double& Volume() { return m_volume; }
				double						Volume() const { return m_volume; }

				TCoordinate& ShapeCenter() { return m_center; }
				const TCoordinate& ShapeCenter() const { return m_center; }

				SingleStiffness& SK() { return m_single_stiffness; }
				const SingleStiffness& SK() const { return m_single_stiffness; }
			private:
				/************************************************************************/
				/* F=K*d                                   */
				/************************************************************************/
				double						m_volume;				//	Modified volume of NODE_j
				TCoordinate					m_center;				//	Shape center of NODE_j
				Eigen::MatrixXd				m_single_stiffness;		//	Single stiffness matrix of Bond_ij
				Eigen::VectorXd				m_force;				//	FORCE vector of Bond_ij
			};

			typedef map<int, TPdBond>		MAP_NJ_TPDBOND;

			typedef TNodeBase TPdNode;
			//	Element
			class TPdElement : public TElementBase
			{
			public:
				TPdElement(vector<TNodeBase>& vecNode) : TElementBase(vecNode) {}
			public:
				void					Dispose()
				{
					TElementBase::Dispose();
				}
				const TPdElement& operator=(const TPdElement& right)
				{
					TElementBase::operator=(right);
					m_map_family_bonds = right.m_map_family_bonds;
					m_pd_paras = right.m_pd_paras;
					m_d_init_volumes = right.m_d_init_volumes;

					return *this;
				}
			public:
				//	The family nodes informations and operations
				int						FamilyElementCount() const { return (int)m_map_family_bonds.size(); }
				MAP_NJ_TPDBOND& FamilyElementBonds() { return m_map_family_bonds; }
				const MAP_NJ_TPDBOND& FamilyElementBonds() const { return m_map_family_bonds; }
				TPdBond& FamilyElementBond(int eId) { return m_map_family_bonds[eId]; }
				void					InsertFamilyElement(int eId, double volume, TCoordinate center) { m_map_family_bonds.insert(pair<int, TPdBond>(eId, TPdBond(volume, center))); }
				void					DeleteFamilyElement(int eId) { m_map_family_bonds.erase(eId); }
				void					ClearFamilyElements() { m_map_family_bonds.clear(); }
			public:
				const TPdCalculateParas& CalParas() const { return m_pd_paras; }
				TPdCalculateParas& CalParas() { return m_pd_paras; }
			public:
				void					InitDamageIndex()
				{
					m_d_init_volumes = 0;
					const MAP_NJ_TPDBOND& familyBonds = FamilyElementBonds();
					for (const pair<int, TPdBond>& mnt : familyBonds)
					{
						m_d_init_volumes += mnt.second.Volume();
					}
				}
			public:
				double					DamageIndex() const
				{
					double current_volumes = 0;
					const MAP_NJ_TPDBOND& familyBonds = FamilyElementBonds();
					for (const pair<int, TPdBond>& mnt : familyBonds)
					{
						current_volumes += mnt.second.Volume();
					}
					return (1.0 - (current_volumes / m_d_init_volumes));
				}

			private:
				MAP_NJ_TPDBOND			m_map_family_bonds;					//	Bond informations
				TPdCalculateParas		m_pd_paras;							//	Calculate parameters of PD Element
				double					m_d_init_volumes;					//	Initial volumes for damage calculation
			};

			typedef TMeshCoreTemplate<TPdNode, TPdElement> TPdMeshCore;

			class TPdDataCollector
			{
			public:
				TPdDataCollector()
				{
					Initialize();
				}
				~TPdDataCollector() {}
			public:
				void						Initialize()
				{
					m_pd_meshcore.Initialize();

					m_vec_parts.clear();
					m_vec_materials.clear();
					m_vec_sections.clear();
					m_vec_curves.clear();

					m_vec_boundary_spc_nodes.clear();
					m_vec_initial_velocity_nodes.clear();
					m_vec_load_node_points.clear();
					m_vec_boundary_prescribed_motions.clear();
					m_vec_crevice.clear();
				}
			public:
				TPdMeshCore& PdMeshCore() { return m_pd_meshcore; }
				const TPdMeshCore& PdMeshCore() const { return m_pd_meshcore; }
			public:
				TPart& Part(int part_id)
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_parts.size()); ++loop)
					{
						if (part_id == m_vec_parts[loop].Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_parts[cur_loop];
				}
				const TPart& Part(int part_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_parts.size()); ++loop)
					{
						if (part_id == m_vec_parts[loop].Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_parts[cur_loop];
				}
				void						AddPart(const TPart& part)
				{
					bool b_exist = false;
					for (int loop = 0; loop < (int)(m_vec_parts.size()); ++loop)
					{
						TPart& cur_part = m_vec_parts[loop];
						//	�����ͬ��ID�ŵ�PART���������Ϣ����
						if (cur_part.Id() == part.Id())
						{
							m_vec_parts[loop] = part;
							b_exist = true;
							break;
						}
					}
					//	���������ͬ��ID��PART�����²���һ��
					if (!b_exist)
					{
						m_vec_parts.push_back(part);
					}
				}
				int							PartCounts() const { return (int)(m_vec_parts.size()); }
				int							PartId(int i_count) const { return m_vec_parts[i_count].Id(); }
				const vector<TPart>& Parts() { return m_vec_parts; }

				TMaterial& Material(int mat_id)
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_materials.size()); ++loop)
					{
						if (mat_id == m_vec_materials[loop].Id())
						{
							cur_loop = loop;
						}
					}
					return m_vec_materials[cur_loop];
				}
				const TMaterial& Material(int mat_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_materials.size()); ++loop)
					{
						if (mat_id == m_vec_materials[loop].Id())
						{
							cur_loop = loop;
						}
					}
					return m_vec_materials[cur_loop];
				}
				bool						MaterialExist(int mat_id)
				{
					bool res = false;
					for (int loop = 0; loop < (int)(m_vec_materials.size()); ++loop)
					{
						if (mat_id == m_vec_materials[loop].Id())
						{
							res = true;
							break;
						}
					}
					return res;
				}
				void						AddMaterial(const TMaterial& mat)
				{
					m_vec_materials.push_back(mat);
				}
				int							MaterialCounts() const { return (int)(m_vec_materials.size()); }
				int							MaterialId(int i_count) const { return m_vec_materials[i_count].Id(); }
				const vector<TMaterial>& Materials() { return m_vec_materials; }

				TSection& Section(int sec_id)
				{
					int cur_loop = -1;
					for (int loop = 0; loop < (int)(m_vec_sections.size()); ++loop)
					{
						if (sec_id == m_vec_sections[loop].Id())
						{
							cur_loop = loop;
						}
					}
					assert(cur_loop >= 0);
					return m_vec_sections[cur_loop];
				}
				const TSection& Section(int sec_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)(m_vec_sections.size()); ++loop)
					{
						if (sec_id == m_vec_sections[loop].Id())
						{
							cur_loop = loop;
						}
					}
					return m_vec_sections[cur_loop];
				}
				bool						SectionExist(int sec_id)
				{
					bool res = false;
					for (int loop = 0; loop < (int)(m_vec_sections.size()); ++loop)
					{
						if (sec_id == m_vec_sections[loop].Id())
						{
							res = true;
							break;
						}
					}
					return res;
				}
				void						AddSection(const TSection& sec)
				{
					m_vec_sections.push_back(sec);
				}
				int							SectionCounts() const { return (int)(m_vec_sections.size()); }
				int							SectionId(int i_count) const { return m_vec_sections[i_count].Id(); }
				const vector<TSection>& Sections() { return m_vec_sections; }

				TCurve& Curve(int cur_id)
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)m_vec_curves.size(); ++loop)
					{
						TCurve& curve = m_vec_curves[loop];
						if (cur_id == curve.Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_curves[cur_loop];
				}
				const TCurve& Curve(int cur_id) const
				{
					int cur_loop = 0;
					for (int loop = 0; loop < (int)m_vec_curves.size(); ++loop)
					{
						const TCurve& curve = m_vec_curves[loop];
						if (cur_id == curve.Id())
						{
							cur_loop = loop;
							break;
						}
					}
					return m_vec_curves[cur_loop];
				}
				bool						CurveExist(int cur_id)
				{
					bool res = false;
					for (int loop = 0; loop < (int)m_vec_curves.size(); ++loop)
					{
						const TCurve& curve = m_vec_curves[loop];
						if (cur_id == curve.Id())
						{
							res = true;
							break;
						}
					}
					return res;
				}
				void						AddCurve(const TCurve& cur)
				{
					m_vec_curves.push_back(cur);
				}
				int							CurveCounts() const { return (int)(m_vec_curves.size()); }
				int							CurveId(int i_count) { return m_vec_curves[i_count].Id(); }
				const vector<TCurve>& Curves() { return m_vec_curves; }
			public:
				//	*BOUNDARY_SPC_NODE informations into PD MODEL
				void						AddBounarySpcNode(const TBoundarySpcNode& boundary)
				{
					m_vec_boundary_spc_nodes.push_back(boundary);
				}
				const vector<TBoundarySpcNode>& BoundarySpcNodes() const { return m_vec_boundary_spc_nodes; }
				vector<TBoundarySpcNode>& BoundarySpcNodes() { return m_vec_boundary_spc_nodes; }

				//	*INITIAL_VELOCITY informations into PD MODEL
				void						AddInitialVelocityNode(const TInitialVelocityNode& initvel)
				{
					m_vec_initial_velocity_nodes.push_back(initvel);
				}
				const vector<TInitialVelocityNode>& InitialVelocityNodes() const { return m_vec_initial_velocity_nodes; }
				vector<TInitialVelocityNode>& InitialVelocityNodes() { return m_vec_initial_velocity_nodes; }

				//	*LOAD_NODE_POINT informations into PD MODEL
				void						AddLoadNodePoint(const TLoadNodePoint& loadnode)
				{
					m_vec_load_node_points.push_back(loadnode);
				}
				const vector<TLoadNodePoint>& LoadNodePoints() const { return m_vec_load_node_points; }
				vector<TLoadNodePoint>& LoadNodePoints() { return m_vec_load_node_points; }

				//	*BOUNDARY_PRESCRIBED_MOTION_NODE information into the PD MODEL 
				void						AddBoundaryPreMotionNode(const TBoundaryPrescribedMotion& motion)
				{
					m_vec_boundary_prescribed_motions.push_back(motion);
				}
				const vector<TBoundaryPrescribedMotion>& BoundaryPrescribedMotions() const { return m_vec_boundary_prescribed_motions; }
				vector<TBoundaryPrescribedMotion>& BoundaryPrescribedMotions() { return m_vec_boundary_prescribed_motions; }

				//	Add the initial crevice into the PD model, and the crevice is *ELEMENT_SEATBELT
				void						AddCrevice(const TCrevice& cre)
				{
					m_vec_crevice.push_back(cre);
				}
				const vector<TCrevice>& Crevice() const { return m_vec_crevice; }
			protected:
				TPdMeshCore							m_pd_meshcore;			//	PD MESH CORE

				vector<TPart>						m_vec_parts;			//	*PART
				vector<TMaterial>					m_vec_materials;		//	*MAT
				vector<TSection>					m_vec_sections;			//	*SECTION
				vector<TCurve>						m_vec_curves;			//	*CURVE

				vector<TBoundarySpcNode>			m_vec_boundary_spc_nodes;				//	*BOUNDARY_SPC_NODE
				vector<TInitialVelocityNode>		m_vec_initial_velocity_nodes;			//	*INITIAL_VELOCITY_NODE
				vector<TLoadNodePoint>				m_vec_load_node_points;					//	*LOAD_NODE_POINT
				vector<TBoundaryPrescribedMotion>	m_vec_boundary_prescribed_motions;		//	*BOUNDARY_PRESCRIBED_MOTION_NODE
				vector<TCrevice>					m_vec_crevice;							//	*ELEMENT_SEATBELT
			};

			//	PD instantiation model
			class TPdModel : public TPdDataCollector
			{
			public:
				//	Update the *PART informations into PD MODEL
				void				UpdatePartInfo()
				{
					double start = clock();
					if (PartCounts() == 0)
					{
						cout << "ERROR: Have no any part information in this LSDYNA file!" << endl;
						return;
					}

					//set<int> eids = m_pd_meshcore.GetElementIdsByAll();
					//for (int eid : eids)
					//{
					//	m_pd_meshcore.AddSeparateElement(eid);						
					//}

					for (int pid = 0; pid < PartCounts(); ++pid)
					{
						TPart& part = Part(PartId(pid));

						set<int> eids = m_pd_meshcore.GetElementIdsByPart(part.Id());
						part.AddElementId(eids);
						//	���õ�Ԫ��THICKNESS
						const TSection& section = Section(part.SectionId());
						string stype = section.Type();
						if (stype == "PLANE_STRESS" ||
							stype == "PLANE_STRAIN")
						{
							double thickness = section.GetSectionValue("THICKNESS");
							for (int eid : eids)
							{
								TPdElement& pd_element = m_pd_meshcore.Element(eid);
								pd_element.Thickness() = thickness;
							}
						}
					}

					//	���²�����Ϣ
					for (const TPart& part : Parts())
					{
						if (!MaterialExist(part.MaterialId()))
						{
							cout << "Material " << part.MaterialId() << " is not exist in this KEYWORD file." << endl;
							return;
						}
						const TMaterial& material = Material(part.MaterialId());
						if (material.Name() == string("MAT_RIGID"))
						{
							int cmo = (int)(material.GetMatValue("CMO"));
							if (cmo == 1)
							{
								int con1 = (int)(material.GetMatValue("CON1"));
								int con2 = (int)(material.GetMatValue("CON2"));
								set<int> nids = part.GetElementIds();
								for (int nid : nids)
								{
									switch (con1)
									{
									case 1:
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										break;
									case 2:
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										break;
									case 3:
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										break;
									case 4:
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										break;
									case 5:
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										break;
									case 6:
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										break;
									case 7:
										AddBounarySpcNode(TBoundarySpcNode(nid, 1));
										AddBounarySpcNode(TBoundarySpcNode(nid, 2));
										AddBounarySpcNode(TBoundarySpcNode(nid, 3));
										break;
									default:
										break;
									}
								}
							}
						}
					}

					double total_time = (clock() - start) / 1000;
					cout << "UpdatePartInfo():\t\t" << total_time << endl;
				}
				void				UpdateFamilyInParts()
				{
					double start = clock();
					//	First, delete all exist family information
					for (int loop = 0; loop < m_pd_meshcore.ElementCount(); ++loop)
					{
						m_pd_meshcore.Element(loop).ClearFamilyElements();
					}

					//	Then, rebuild the family information for every node
					for (const TPart& part : Parts())
					{
						int sid = part.SectionId();
						const TSection& section = Section(sid);
						string stype = section.Type();

						const set<int> eleIds = part.GetElementIds();
						// 
						parallel_for_each(eleIds.begin(), eleIds.end(), [&](int ei)
							{
								TPdElement& element_i = m_pd_meshcore.Element(ei);
								TCoordinate coor_i = element_i.CoordinateInElement(0, 0);

								double dx = element_i.SideLength();
								//	Extend the search range to include all nodes

								double dis_for_juge = 0;
								if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
								{
									dis_for_juge = HORIZON + dx + ERR_VALUE;
								}
								else
								{
									double m = RATIO_OF_HORIZON_MESHSIZE + 1.0;
									dis_for_juge = m * dx + ERR_VALUE;
								}

								for (int ej : eleIds)
								{
									if (ei == ej)
										continue;

									TPdElement& element_j = m_pd_meshcore.Element(ej);
									TCoordinate coor_j = element_j.CoordinateInElement(0, 0);

									if ((abs(coor_i.x() - coor_j.x()) < dis_for_juge) &&
										(abs(coor_i.y() - coor_j.y()) < dis_for_juge) &&
										(abs(coor_i.z() - coor_j.z()) < dis_for_juge))
									{
										double xi = Distance_2pt(coor_i, coor_j);
										/*输出ei的编号 swan 2021-01-22 */
										//std::cout << "elem i : " << ei+1 << "\t\t->nj : ";

										if (xi < dis_for_juge)
										{
											Vector3d cij;
											// 体积修正
											double intersect_volume = CalModifiedVolume(element_i, element_j, stype, cij);
											element_i.InsertFamilyElement(ej, intersect_volume, cij);

											//std::cout << "elem j: " << ej+1 << "\t-->volume: " << intersect_volume << endl;

											/*输出ni nj的坐标 swan 2021-01-22
											cout << ei + 1 << "\t (";
											cout << coor_i.x() << ", ";
											cout << coor_i.y() << ", ";
											cout << coor_i.z() << ")\t--> ";
											cout << ej + 1 << "\t (";
											cout << coor_j.x() << ", ";
											cout << coor_j.y() << ", ";
											cout << coor_j.z() << ")\n";*/
										}
									}
								}
							});
					}
					double total_time = (clock() - start) / 1000;
					cout << "UpdateFamilyInParts():\t\t" << total_time << endl;
				}
			public:
				//	Refresh the *BOUNDARY_SPC_NODE informations into PD MODEL
				void				RefreshBoundarySpcInfo()
				{
					for (const TBoundarySpcNode& tbsn : m_vec_boundary_spc_nodes)
					{
						int nid = tbsn.Id();
						int dof = tbsn.Dof() - 1;
						m_pd_meshcore.Node(nid).Displacement()(dof) = 0;
					}
				}
				//	Refresh the *INITIAL_VELOCITY informations into PD MODEL
				void				RefreshInitVelocityInfo()
				{
					for (const TInitialVelocityNode& tivn : m_vec_initial_velocity_nodes)
					{
						int nid = tivn.Id();
						m_pd_meshcore.Node(nid).Velocity() = tivn.Velocity();
					}
				}
				//	Refresh the *LOAD_NODE_POINT informations into PD MODEL
				void				RefreshLoadNodeInfo(double cur_time)
				{
					for (const TLoadNodePoint& tlnp : m_vec_load_node_points)
					{
						int nid = tlnp.Id();
						int curid = tlnp.Lcid();

						double value_force = 0;
						if (CurveExist(curid))
						{
							TCurve& curve = Curve(curid);
							value_force = curve.GetValueByX(cur_time) * tlnp.Sf();
						}
						else
						{
							//	û�����ߣ������������Ϊ��ֵ
							value_force = tlnp.Sf();
						}

						switch (tlnp.Dof())
						{
						case 1:
						{
							m_pd_meshcore.Node(nid).OuterForce().x() = value_force;
							break;
						}
						case 2:
						{
							m_pd_meshcore.Node(nid).OuterForce().y() = value_force;
							break;
						}
						case 3:
						{
							m_pd_meshcore.Node(nid).OuterForce().z() = value_force;
							break;
						}
						default:
							break;
						}
					}
				}
				//	Refresh the *BOUNDARY_PRESCRIBED_MOTION_NODE information into the PD MODEL 
				void				RefreshPreMotionInfo(double cur_time)
				{
					for (const TBoundaryPrescribedMotion& bpm : m_vec_boundary_prescribed_motions)
					{
						string type = bpm.MotionType();
						set<int> nids;
						nids.clear();
						if (type == "RIGID")
						{
							int partId = bpm.Id();
							TPart& part = Part(partId);
							nids = part.GetElementIds();
						}
						else if (type == "NODE")
						{
							nids.insert(bpm.Id());
						}
						for (int nid : nids)
						{
							int curid = bpm.Lcid();
							double value = 0;
							if (CurveExist(curid))
							{
								TCurve& curve = Curve(curid);
								value = curve.GetValueByX(cur_time) * bpm.Sf();
							}
							else
							{
								value = bpm.Sf();
							}

							switch (bpm.Dof())
							{
								//	DOF=1 -> X translation
							case 1:
							{
								switch (bpm.Vda())
								{
									//	VDA = 0 -> VELOCITY
								case 0:
								{
									m_pd_meshcore.Node(nid).Displacement().x() = value * cur_time;
									break;
								}
								//	VDA = 1 -> ACCELERATION
								case 1:
								{
									m_pd_meshcore.Node(nid).Displacement().x() = 0.5 * value * cur_time * cur_time;
								}
								//	VDA = 2 -> DISPLACEMENT
								case 2:
								{
									m_pd_meshcore.Node(nid).Displacement().x() = value;
								}
								default:
									break;
								}
								break;
							}
							//	DOF=2 -> Y translation
							case 2:
							{
								switch (bpm.Vda())
								{
								case 0:
								{
									m_pd_meshcore.Node(nid).Displacement().y() = value * cur_time;
									break;
								}
								case 1:
								{
									m_pd_meshcore.Node(nid).Displacement().y() = 0.5 * value * cur_time * cur_time;
								}
								case 2:
								{
									m_pd_meshcore.Node(nid).Displacement().y() = value;
								}
								default:
									break;
								}
								break;
							}
							//	DOF=3 -> Z translation
							case 3:
							{
								switch (bpm.Vda())
								{
								case 0:
								{
									m_pd_meshcore.Node(nid).Displacement().z() = value * cur_time;
									break;
								}
								case 1:
								{
									m_pd_meshcore.Node(nid).Acceleration().z() = 0.5 * value * cur_time * cur_time;
								}
								case 2:
								{
									m_pd_meshcore.Node(nid).Displacement().z() = value;
								}
								default:
									break;
								}
								break;
							}
							default:
								break;
							}
						}
					}
				}
			public:
				//	Generate initial crevice for PD model
				void				GenerateInitCrevice()
				{
					double start = clock();

					for (const TCrevice& crevice : m_vec_crevice)
					{
						for (const TPart& part : Parts())
						{
							const set<int> eids = part.GetElementIds();
							parallel_for_each(eids.begin(), eids.end(), [&](int ei) {
								TPdElement& element_i = m_pd_meshcore.Element(ei);
								const TCoordinate& ei_coord = element_i.CoordinateInElement(0, 0);
								MAP_NJ_TPDBOND& familyBonds = element_i.FamilyElementBonds();

								for (MAP_NJ_TPDBOND::iterator iter = familyBonds.begin();
									iter != familyBonds.end();)
								{
									int ej = iter->first;
									TPdElement& element_j = m_pd_meshcore.Element(ej);
									const TCoordinate& ej_coord = element_j.CoordinateInElement(0, 0);

									bool res = IsTowLineIntersect_xy_plane<Vector3d>(ei_coord, ej_coord, crevice.Start(), crevice.End());
									if (res)
									{
										++iter;
										element_i.DeleteFamilyElement(element_j.Id());
									}
									else
									{
										++iter;
									}
								}
								});
						}
					}

					double total_time = (clock() - start) / 1000;
					cout << "GenerateInitCrevice():\t\t" << total_time << endl;
				}
				//	Generate initial Damage Value for all PD nodes
				void				GenerateInitDamage()
				{
					double start = clock();
					//	Set FamilyNodeCount for Initial Damage Value
					set<int> eleIds = m_pd_meshcore.GetElementIdsByAll();
					for (int i : eleIds)
					{
						TPdElement& element = m_pd_meshcore.Element(i);
						element.InitDamageIndex();
					}

					double total_time = (clock() - start) / 1000;
					cout << "GenerateInitDamage():\t\t" << total_time << endl;
				}
			public:
				//	Influence index
				double				CalInfluenceIndex(double idist, double horizon)
				{
					return 1.0;
					/*return 1 - (idist / horizon);*/
					/*return pow(1 - pow( (idist / horizon), 2), 2);*/
				}
				double				CalModifiedVolume(const TPdElement& element_i, const TPdElement& element_j, string stype, Vector3d& center)
				{
					double volume_of_j = 0;

					double Ri = 0;
					if (DLUT::SAE::PERIDYNAMIC::USE_CONSTANT_HORIZON)
					{
						Ri = HORIZON;
					}
					else
					{
						Ri = element_i.SideLength() * RATIO_OF_HORIZON_MESHSIZE;
					}

					if (stype == "PLANE_STRESS" ||
						stype == "PLANE_STRAIN")
					{
						int partId = element_i.PartId();
						const TPart& part = Part(partId);
						int sectionId = part.SectionId();
						const TSection& section = Section(sectionId);
						double thickness = section.GetSectionValue("THICKNESS");

						vector<int> nids = element_j.NodeIds();
						int nCount = (int)(nids.size());

						Vector3d pt1 = m_pd_meshcore.Node(nids[0]).Coordinate();
						Vector3d pt2 = m_pd_meshcore.Node(nids[1]).Coordinate();
						Vector3d pt3 = m_pd_meshcore.Node(nids[2]).Coordinate();

						if (element_j.MeshElementType() == TRIANGLE_ELEMENT)
						{
							volume_of_j = Calculate_Intersect_Area_Circle_Triangle(element_i.CoordinateInElement(0, 0), Ri, pt1, pt2, pt3, center) * thickness;
						}
						else
						{
							Vector3d pt4 = m_pd_meshcore.Node(nids[3]).Coordinate();
							volume_of_j = Calculate_Intersect_Area_Circle_Quadrangle(element_i.CoordinateInElement(0, 0), Ri, pt1, pt2, pt3, pt4, center) * thickness;
						}
					}

					return volume_of_j;
				}
			};

			/************************************************************************/
			/* �� vector< map<int, double> > ת���� ϡ�����                        */
			/************************************************************************/
			void TransVecMap2SparseMatrix(const vector< map<int, double> >& vec_map_matrix, SparseMatrix<double>& sparse_matrix)
			{
				const double ERROR_FOR_SM = ERR_VALUE;
				int nCount = (int)(vec_map_matrix.size());
				sparse_matrix.resize(nCount, nCount);

				vector< Triplet<double> > tri;
				tri.clear();
				for (int i = 0; i < (int)(vec_map_matrix.size()); ++i)
				{
					for (const pair<int, double>& j_v : vec_map_matrix[i])
					{
						if (abs(j_v.second) > ERROR_FOR_SM)
						{
							tri.push_back(Triplet<double>(i, j_v.first, j_v.second));
						}
					}
				}
				sparse_matrix.setFromTriplets(tri.begin(), tri.end());
			}
			/************************************************************************/
			/* �� vector< map<int, double> > ת���� Matrix����                      */
			/************************************************************************/
			void TransVecMap2Matrix(const vector< map<int, double> >& vec_map_matrix, MatrixXd& matrix)
			{
				int nCount = (int)(vec_map_matrix.size());
				matrix.resize(nCount, nCount);
				matrix.setZero();
				for (int i = 0; i < (int)(vec_map_matrix.size()); ++i)
				{
					for (const pair<int, double>& j_v : vec_map_matrix[i])
					{
						matrix(i, j_v.first) = j_v.second;
					}
				}
			}
		} //	end of namespace PERIDYNAMIC
	} // end of namespace SAE
} // end of namespace DLUT

#endif