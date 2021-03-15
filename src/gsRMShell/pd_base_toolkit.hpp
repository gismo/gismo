/*==============================================================================

Copyright 2017 Dalian University of Technology .
All rights reserved

--------------------------------------------------------------------=
-- Please append file description informations here --
Some base functions for PD
--------------------------------------------------------------------=
Date            Name                    Description of Change
2017/12/20		Zheng Guojun			Create
2019/11/09		Zheng Guojun			Added Calculate_Intersect_Area_Circle_Quadrangle()
2019/11/10		Zheng Guojun			Modified Normalizer() to avoid division by zero
$HISTORY$
--------------------------------------------------------------------=
*/
#ifndef DLUT_SAE_PERIDYNAMIC_PD_BASE_TOOLKIT_HXX_20171220
#define DLUT_SAE_PERIDYNAMIC_PD_BASE_TOOLKIT_HXX_20171220

#include <vector>
#include <set>
#include <math.h>
#include <map>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <iterator>

using namespace std;
namespace DLUT
{
	namespace SAE
	{
		namespace PERIDYNAMIC
		{
			const static double ERR_VALUE = 1E-10;
			const static double PI = acos(-1);

			template<typename Type, typename VALUE>
			bool IsValueExistInMultiMap(multimap<Type, VALUE>& _mmap, pair<Type, VALUE> _value)
			{
				bool b_find_value = false;
				auto range = _mmap.equal_range(_value.first);
				for (auto i = range.first; i != range.second; ++i)
				{
					if (i->second == _value.second)
					{
						b_find_value = true;
						break;
					}
				}
				return b_find_value;
			}
			//	计算最大值
			template<typename Type>
			Type Max(const Type& lhs, const Type& rhs)
			{
				return lhs > rhs ? lhs : rhs;
			}
			//	计算最小值
			template<typename Type>
			Type Min(const Type& lhs, const Type& rhs)
			{
				return lhs < rhs ? lhs : rhs;
			}
			//	计算两个点之间的距离
			template<typename Type>
			inline double Distance_2pt(const Type& pt1, const Type& pt2)
			{
				return sqrt((pt1.x() - pt2.x()) * (pt1.x() - pt2.x()) +
					(pt1.y() - pt2.y()) * (pt1.y() - pt2.y()) +
					(pt1.z() - pt2.z()) * (pt1.z() - pt2.z()));
			}
			//	计算|a|的值
			template<typename Type>
			inline double Module(const Type& a)
			{
				return sqrt((a.x() * a.x() + a.y() * a.y() + a.z() * a.z()));
			}
			//	计算向量a.b(a和b的点乘)
			template<typename Type>
			double Point_Multi(const Type& a, const Type& b)
			{
				return (a.x() * b.x() + a.y() * b.y() + a.z() * b.z());
			}
			//	计算向量a*b(a和b的叉乘)
			template<typename RetType, typename Type>
			RetType Fork_Multi(const Type& a, const Type& b)
			{
				double x = a.y() * b.z() - a.z() * b.y();
				double y = a.z() * b.x() - a.x() * b.z();
				double z = a.x() * b.y() - a.y() * b.x();
				return RetType(x, y, z);
			}
			//	计算向量的并矢(张量积)
			template<typename RetType, typename Type>
			RetType Tensor_Multi(const Type& a, const Type& b)
			{
				RetType res;
				for (int i = 0; i < 3; ++i)
				{
					for (int j = 0; j < 3; ++j)
					{
						res(i, j) = a(i) * b(j);
					}
				}
				return res;
			}
			//	归一化处理向量a
			template<typename Type>
			void Normalizer(Type& a)
			{
				const double length = Module(a);
				if (length > ERR_VALUE)
				{
					a.x() = a.x() / length;
					a.y() = a.y() / length;
					a.z() = a.z() / length;
				}
			}
			//	计算以pt1, pt2, pt3构成的夹角的cos值, 其中pt2为顶角
			template<typename Type>
			double Calculate_COS(const Type& pt1, const Type& pt2, const Type& pt3)
			{
				Type v1 = pt1 - pt2;
				Type v2 = pt3 - pt2;
				return (Point_Multi(v1, v2) / (Module(v1) * Module(v2)));
			}
			//	计算以向量v1, v2构成的夹角的cos值
			template<typename Type>
			double Calculate_COS(const Type& v1, const Type& v2)
			{
				return (Point_Multi(v1, v2) / (Module(v1) * Module(v2)));
			}
			//	计算三点构成平面的法线
			template<typename Type>
			Type Calculate_Normal(const Type& p1, const Type& p2, const Type& p3)
			{
				Type v1 = p2 - p1;
				Type v2 = p3 - p2;
				return Fork_Multi<Type>(v1, v2);
			}
			//	计算三角形面积
			template<typename Type>
			double Calculate_Area_Triangle(const Type& p1, const Type& p2, const Type& p3)
			{
				double area = 0;

				double a = Distance_2pt(p1, p2);
				double b = Distance_2pt(p2, p3);
				double c = Distance_2pt(p3, p1);
				double p = (a + b + c) / 2.0;

				double s = p * (p - a) * (p - b) * (p - c);
				if (s < 0.0)
				{
					area = 0;
				}
				else
				{
					area = sqrt(s);
				}

				return area;
			}
			//	计算三角形形心
			template<typename Type>
			Type Calculate_ShapeCenter_Triangle(const Type& p1, const Type& p2, const Type& p3)
			{
				Type center;
				center.x() = (p1.x() + p2.x() + p3.x()) / 3.0;
				center.y() = (p1.y() + p2.y() + p3.y()) / 3.0;
				center.z() = (p1.z() + p2.z() + p3.z()) / 3.0;

				return center;
			}
			//	计算扇形面积
			template<typename Type>
			double Calculate_Area_Sector(const Type& c, const Type& p1, const Type& p2)
			{
				double theta = acos(Calculate_COS(p1, c, p2));
				double r = Distance_2pt(c, p1);
				double res = theta * r * r / 2.0;

				return res;
			}
			//	计算扇形形心
			template<typename Type>
			Type Calculate_ShapeCenter_Sector(const Type& c, const Type& p1, const Type& p2)
			{
				Type center;

				Type p3 = (p1 + p2) / 2.0;
				Type v3 = p3 - c;
				Normalizer(v3);
				double r = Distance_2pt(c, p1);
				double theta = acos(Calculate_COS(p1, c, p2));

				double alpha = theta / 2.0;
				double xc = 2 * r * sin(alpha) / (3.0 * alpha);

				center = v3 * xc + c;

				return center;
			}
			//	计算弓形面积
			template<typename Type>
			double Calculate_Area_Arch(const Type& c, const Type& p1, const Type& p2)
			{
				return Calculate_Area_Sector(c, p1, p2) - Calculate_Area_Triangle(c, p1, p2);
			}
			//	计算弓形形心
			template<typename Type>
			Type Calculate_ShapeCenter_Arch(const Type& c, const Type& p1, const Type& p2)
			{
				Type center;

				Type p3 = (p1 + p2) / 2.0;
				Type v3 = p3 - c;
				Normalizer(v3);
				double r = Distance_2pt(c, p1);
				double theta = acos(Calculate_COS(p1, c, p2));

				double alpha = theta / 2.0;
				double xc = 4 * r * pow(sin(alpha), 3) / (3 * (2 * alpha - sin(2 * alpha)));

				center = v3 * xc + c;

				return center;
			}
			//	判断两个点是否重合
			template<typename Type>
			bool IsEqual(const Type& pt1, const Type& pt2)
			{
				return ((pt1.x() - pt2.x()) < ERR_VALUE) &&
					((pt1.y() - pt2.y()) < ERR_VALUE) &&
					((pt1.z() - pt2.z()) < ERR_VALUE);
			}
			//	判断点是否在线段内部
			template<typename Type>
			bool IsPointOnLine(const Type& point, const Type& lpt, const Type& rpt)
			{
				double a = Distance_2pt(point, lpt);
				double b = Distance_2pt(point, rpt);
				double c = Distance_2pt(lpt, rpt);
				if (abs(c - a - b) < ERR_VALUE)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			//	全局坐标点转换到局部坐标点
			template<typename Type, typename TNormal>
			Type GlobalToLocal(const Type& pt, const TNormal& nx, const TNormal& ny, const TNormal& nz)
			{
				Type res;
				res.x() = pt.x() * nx.x() + pt.y() * ny.x() + pt.z() * nz.x();
				res.y() = pt.x() * nx.y() + pt.y() * ny.y() + pt.z() * nz.y();
				res.z() = pt.x() * nx.z() + pt.y() * ny.z() + pt.z() * nz.z();
				return res;
			}
			//	点投影到平面上
			template<typename Type, typename TNormal>
			Type PointProjectToPlane(const Type& pt, const TNormal& planeNormal, const Type& planeOrigial)
			{
				double a = planeNormal.x();
				double b = planeNormal.y();
				double c = planeNormal.z();

				double t = a * planeOrigial.x() + b * planeOrigial.y() + c * planeOrigial.z() - (a * pt.x() + b * pt.y() + c * pt.z());
				t = t / (a * a + b * b + c * c);
				Type res;
				res.x() = pt.x() + a * t;
				res.y() = pt.y() + b * t;
				res.z() = pt.z() + c * t;

				return res;
			}
			//	求点到经过两点直线的距离
			template<typename Type>
			double Distance_pt2line(const Type& Pt, const Type& P1, const Type& P2)
			{
				double m = P2.x() - P1.x();
				double n = P2.y() - P1.y();
				double p = P2.z() - P1.z();
				double t = (m * (Pt.x() - P1.x()) + n * (Pt.y() - P1.y()) + p * (Pt.z() - P1.z())) / (m * m + n * n + p * p);

				Type Pt_c;
				Pt_c.x() = m * t + P1.x();
				Pt_c.y() = n * t + P1.y();
				Pt_c.z() = p * t + P1.z();

				//	如果垂足在线段内部，返回实际距离值；如果垂直不在线段内部，返回一个极大值，表示无意义
				double distance = 1E20;
				if (IsPointOnLine(Pt_c, P1, P2))
				{
					distance = Distance_2pt(Pt, Pt_c);
				}

				return distance;
			}
			//	设置平面多边形的方向，顺时针还是逆时针，bcockwise表示是否顺时针
			template<typename Type>
			void SetPolygonDirection(Type* polygon, int vertex_num, bool bcockwise = true)
			{
				//	错误保护，防止环上面的点构不成平面
				if (vertex_num < 3)
					return;

				int xMinId = 0;
				set<int> xMinId_set;
				double xMin = 1e10;
				for (int i = 0; i < vertex_num; ++i)
				{
					if (xMin > (polygon[i].x()))
					{
						xMin = polygon[i].x();
					}
				}
				for (int j = 0; j < vertex_num; ++j)
				{
					if (polygon[j].x() < (xMin + ERR_VALUE))
					{
						xMinId_set.insert(j);
						xMinId = j;
					}
				}
				double yMax = 1e-10;
				for (set<int>::iterator iter = xMinId_set.begin(); iter != xMinId_set.end(); ++iter)
				{
					int id = *iter;
					if (yMax < polygon[id].y())
					{
						yMax = polygon[id].y();
						xMinId = id;
					}
				}
				Type v1(polygon[xMinId], polygon[(xMinId - 1 + vertex_num) % vertex_num]);
				Type v2(polygon[xMinId], polygon[(xMinId + 1) % vertex_num]);
				double vz = v1.x() * v2.y() - v1.y() * v2.x();

				//	type = 0表示要求顺时针
				//	vz > 0表示该环是顺时针

				//	要求是顺时针，但实际是逆时针，或者是要求是逆时针，但实际是顺时针
				if ((bcockwise && vz < 0) || (!bcockwise && vz > 0))
				{
					reverse(polygon, polygon + vertex_num);
				}
			}
			//	判断点是否在多边形内部，只适合在xy平面上
			//	0表示在内部，1表示在外部，2表示在多边形边界
			enum PolygonPosition
			{
				OutOfPolygon = 1, InPolygon, OnPolygon
			};
			template<typename Type>
			PolygonPosition IsPointInPolygon_plane(const Type& point, const Type* polygon, int vertexNum)
			{
				PolygonPosition res = InPolygon;
				bool oddNode = 0;
				double x = point.x();
				double y = point.y();

				for (int i = 0; i < vertexNum; i++)
				{
					int j = (i + 1) % vertexNum;
					double x1 = polygon[i].x();
					double y1 = polygon[i].y();

					double x2 = polygon[j].x();
					double y2 = polygon[j].y();

					if (((y1 < y) && (y2 >= y)) || ((y1 >= y) && (y2 < y)) && ((x1 < x) || (x2 < x)))
					{
						double x_intersection = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
						if (abs(x_intersection - x) < 1E-6)
						{
							res = OnPolygon;
							return res;
						}
						if (x_intersection < x)
						{
							oddNode = !oddNode;
						}
					}
				}

				if (oddNode)
				{
					res = InPolygon;
				}
				else
				{
					res = OutOfPolygon;
				}
				return res;
			}
			//	直线与多边形的交点
			template<typename Type>
			Type IntersectOfLineAndPolygon(const Type& pt1, const Type& pt2, const Type* polygon, int vertexNum)
			{
				assert(vertexNum >= 3);
				Type pn = Fork_Multi<Type, Type>(polygon[1] - polygon[0], polygon[2] - polygon[1]);
				Type ln = pt2 - pt1;

				double a = pn.x();
				double b = pn.y();
				double c = pn.z();
				double d = -(a * polygon[0].x() + b * polygon[0].y() + c * polygon[0].z());

				double P1_D = -(a * pt1.x() + b * pt1.y() + c * pt1.z() + d);
				double P1_D2 = a * ln.x() + b * ln.y() + c * ln.z();
				Type pt_inter(pt1);
				if (fabs(P1_D2) > ERR_VALUE)
				{
					double n = P1_D / P1_D2;
					pt_inter.x() = pt1.x() + n * (pt2.x() - pt1.x());
					pt_inter.y() = pt1.y() + n * (pt2.y() - pt1.y());
					pt_inter.z() = pt1.z() + n * (pt2.z() - pt1.z());
				}
				return pt_inter;
			}
			// 多边形与过point的X轴平行线的最小距离
			template<typename Type>
			double MinDisPtToPolygonAlongX(const Type& point, const Type* polygon, int vertexNum)
			{
				double x = point.x();
				double y = point.y();

				vector<Type> vecIntersection;
				vecIntersection.clear();

				for (int i = 0; i < vertexNum; i++)
				{
					int j = (i + 1) % vertexNum;
					double x1 = polygon[i].x();
					double y1 = polygon[i].y();

					double x2 = polygon[j].x();
					double y2 = polygon[j].y();

					if (((y1 < y) && (y2 > y)) || ((y1 > y) && (y2 < y)))
					{
						double x_intersection = x1 + (y - y1) * (x2 - x1) / (y2 - y1);
						Type pt_sec(x_intersection, y);
						vecIntersection.push_back(pt_sec);
					}
				}
				double res = 1E10;
				for (int i = 0; i < (int)(vecIntersection.size()); i++)
				{
					double dis = Distance_2pt(point, vecIntersection[i]);
					if (res > dis)
					{
						res = dis;
					}
				}

				return res;
			}
			// 多边形与过point的Y轴平行线的最小距离
			template<typename Type>
			double MinDisPtToPolygonAlongY(const Type& point, const Type* polygon, int vertexNum)
			{
				double x = point.x();
				double y = point.y();

				vector<Type> vecIntersection;
				vecIntersection.clear();

				for (int i = 0; i < vertexNum; i++)
				{
					int j = (i + 1) % vertexNum;
					double x1 = polygon[i].x();
					double y1 = polygon[i].y();

					double x2 = polygon[j].x();
					double y2 = polygon[j].y();

					if (((x1 < x) && (x2 > x)) || ((x1 > x) && (x2 < x)))
					{
						double y_intersection = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
						Type pt_sec(x, y_intersection);
						vecIntersection.push_back(pt_sec);
					}
				}

				double res = 1E10;
				for (int i = 0; i < (int)(vecIntersection.size()); i++)
				{
					double dis = Distance_2pt(point, vecIntersection[i]);
					if (res > dis)
					{
						res = dis;
					}
				}

				return res;
			}
			template<typename Type>
			bool IsTwoRectIntersect_xy_plane(const Type& P1, const Type& P2, const Type& Q1, const Type& Q2)
			{
				//	首先判断X轴方向上面是否能相交
				bool ret = false;
				double xpMin = Min(P1.x(), P2.x());
				double xpMax = Max(P1.x(), P2.x());
				double xqMin = Min(Q1.x(), Q2.x());
				double xqMax = Max(Q1.x(), Q2.x());
				//	如果在X轴方向能够相交，则一定有一线段的最小x值夹在另一线段中间
				bool xret = ((xpMin < (xqMax + ERR_VALUE)) && (xpMin > (xqMin - ERR_VALUE)) ||
					(xqMin < (xpMax + ERR_VALUE)) && (xqMin > (xpMin - ERR_VALUE)));
				//	如果在x轴上面能相交，还需要检测在Y轴上面是否能相交
				if (xret)
				{
					double ypMin = Min(P1.y(), P2.y());
					double ypMax = Max(P1.y(), P2.y());
					double yqMin = Min(Q1.y(), Q2.y());
					double yqMax = Max(Q1.y(), Q2.y());
					bool yret = ((ypMin < (yqMax + ERR_VALUE)) && (ypMin > (yqMin - ERR_VALUE)) ||
						(yqMin < (ypMax + ERR_VALUE)) && (yqMin > (ypMin - ERR_VALUE)));

					ret = xret || yret;
				}

				return ret;
			}
			//	判断点是否在以一条线段为对角线的矩形内部
			template<typename Type>
			bool IsPointInRect_xy_plane(const Type& pt, const Type& Q1, const Type& Q2)
			{
				double xMax = Max(Q1.x(), Q2.x());
				double xMin = Min(Q1.x(), Q2.x());
				double yMax = Max(Q1.y(), Q2.y());
				double yMin = Min(Q1.y(), Q2.y());
				return (pt.x() < (xMax + ERR_VALUE)) && (pt.x() > (xMin - ERR_VALUE)) &&
					(pt.y() < (yMax + ERR_VALUE)) && (pt.y() > (yMin - ERR_VALUE));
			}
			//	根据两点求直线解析方程Ax+By+C=0;
			//	仅适合在x_y平面上
			template<typename Type>
			void MakeLine(const Type& P1, const Type& P2, double& A, double& B, double& C)
			{
				int sign = 1;
				A = P2.y() - P1.y();
				if (A < 0)
				{
					sign = -1;
					A = -A;
				}
				B = (P1.x() - P2.x()) * sign;
				C = (P1.y() * P2.x() - P1.x() * P2.y()) * sign;
			}
			//	判断x_y平面上面两条线段是否相交
			template<typename Type>
			bool IsTowLineIntersect_xy_plane(const Type& P1, const Type& P2, const Type& Q1, const Type& Q2)
			{
				bool res = false;
				double A1, B1, C1;
				double A2, B2, C2;
				Type insec;

				MakeLine(P1, P2, A1, B1, C1);
				MakeLine(Q1, Q2, A2, B2, C2);
				double k = A1 * B2 - A2 * B1;
				if (abs(k) < ERR_VALUE)
				{
					//	两条直线平行，则在直线时，可能节点相交
					if (IsPointOnLine(Q1, P1, P2) || IsPointOnLine(Q2, P1, P2))
					{
						res = true;
					}
				}
				else
				{
					insec.x() = (C2 * B1 - C1 * B2) / k;
					insec.y() = (C1 * A2 - C2 * A1) / k;
					insec.z() = 0;
					if (IsPointOnLine(insec, P1, P2) && IsPointOnLine(insec, Q1, Q2))
					{
						res = true;
					}
				}

				return res;
			}

			//	求点在过两点直线的垂足
			template<typename Type>
			Type Perpendicular_pt(const Type& Pt, const Type& P1, const Type& P2)
			{
				double m = P2.x() - P1.x();
				double n = P2.y() - P1.y();
				double p = P2.z() - P1.z();
				double t = (m * (Pt.x() - P1.x()) + n * (Pt.y() - P1.y()) + p * (Pt.z() - P1.z())) / (m * m + n * n + p * p);

				Type Pt_c;
				Pt_c.x() = m * t + P1.x();
				Pt_c.y() = n * t + P1.y();
				Pt_c.z() = p * t + P1.z();

				return Pt_c;
			}
			//	求圆与线段的交点
			template<typename Type>
			vector<Type> Calculate_Intersect_Point_Cricle_Segment(const Type& center, double r, const Type& p1, const Type& p2)
			{
				vector<Type> res;
				res.clear();

				double x0 = center.x();
				double y0 = center.y();
				double z0 = center.z();

				double x1 = p1.x();
				double y1 = p1.y();
				double z1 = p1.z();

				double x2 = p2.x();
				double y2 = p2.y();
				double z2 = p2.z();

				double a = pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2);
				double b = 2 * ((x2 - x1) * (x1 - x0) + (y2 - y1) * (y1 - y0) + (z2 - z1) * (z1 - z0));
				double c = pow((x1 - x0), 2) + pow((y1 - y0), 2) + pow((z1 - z0), 2) - pow(r, 2);

				double delta = b * b - 4 * a * c;
				if (delta > 0)
				{
					double k1 = (-b - sqrt(delta)) / (2 * a);
					//	k1必须在[0,1]区间，表示交点在线段之间
					if (k1 > -ERR_VALUE && k1 < 1 + ERR_VALUE)
					{
						Type s1;
						s1.x() = x1 + k1 * (x2 - x1);
						s1.y() = y1 + k1 * (y2 - y1);
						s1.z() = z1 + k1 * (z2 - z1);

						res.push_back(s1);
					}
					double k2 = (-b + sqrt(delta)) / (2 * a);
					//	k2必须在[0,1]区间，表示交点在线段之间
					if (k2 > -ERR_VALUE && k2 < 1 + ERR_VALUE)
					{
						Type s2;
						s2.x() = x1 + k2 * (x2 - x1);
						s2.y() = y1 + k2 * (y2 - y1);
						s2.z() = z1 + k2 * (z2 - z1);

						res.push_back(s2);
					}
				}

				return res;
			}
			//	计算两个圆相交的面积
			template<typename Type>
			double Calculate_Intersect_Area_Two_Circle(const Type& c1, double r1, const Type& c2, double r2)
			{
				double r = Min(r1, r2);
				double R = Max(r1, r2);
				double d = Distance_2pt(c1, c2);
				if ((d + ERR_VALUE) > (r + R))
				{
					return 0;
				}
				else if ((d + r) < (R + ERR_VALUE))
				{
					return PI * r * r;
				}
				else
				{
					double ang1 = acos((r * r + d * d - R * R) / (2.0 * r * d));
					double ang2 = acos((R * R + d * d - r * r) / (2.0 * R * d));

					return ang1 * r * r + ang2 * R * R - r * sin(ang1) * d;
				}
			}
			//	计算两个球相交的体积
			template<typename Type>
			double Calculate_Intersect_Volume_Two_Sphere(const Type& c1, double r1, const Type& c2, double r2)
			{
				double r = Min(r1, r2);
				double R = Max(r1, r2);
				double d = Distance_2pt(c1, c2);
				if ((d + ERR_VALUE) > (r + R))
				{
					return 0;
				}
				else if ((d + r) < (R + ERR_VALUE))
				{
					return (4.0 / 3.0) * PI * pow(r, 3);
				}
				else
				{
					double ang1 = acos((r * r + d * d - R * R) / (2.0 * r * d));
					double ang2 = acos((R * R + d * d - r * r) / (2.0 * R * d));

					return (4.0 / 3.0) * ang1 * pow(r, 3) + (4.0 / 3.0) * ang2 * pow(R, 3) - (1.0 / 3.0) * PI * pow(r * sin(ang1), 2) * d;
				}
			}
			//	计算圆与三角形相交的面积
			template<typename Type>
			double Calculate_Intersect_Area_Circle_Triangle(const Type& c, double r, const Type& p1, const Type& p2, const Type& p3, Type& center)
			{
				double area = 0;
				center = (p1 + p2 + p3) / 3;

				vector<Type> pts_in_circle;
				vector<Type> pts_out_circle;
				if (Distance_2pt(c, p1) < r)
				{
					pts_in_circle.push_back(p1);
				}
				else
				{
					pts_out_circle.push_back(p1);
				}
				if (Distance_2pt(c, p2) < r)
				{
					pts_in_circle.push_back(p2);
				}
				else
				{
					pts_out_circle.push_back(p2);
				}
				if (Distance_2pt(c, p3) < r)
				{
					pts_in_circle.push_back(p3);
				}
				else
				{
					pts_out_circle.push_back(p3);
				}

				if ((int)(pts_in_circle.size()) == 0)
				{
					//	对 d1  d2  d3进行排序
					map<double, vector<Type> > d_intsect_pts;
					d_intsect_pts.clear();
					double d1 = Distance_pt2line(c, p2, p3);
					double d2 = Distance_pt2line(c, p3, p1);
					double d3 = Distance_pt2line(c, p1, p2);
					if (d1 < r - ERR_VALUE)
					{
						vector<Type> pts = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, p2, p3);
						if (pts.size() == 2)
						{
							d_intsect_pts.insert(make_pair(d1, pts));
						}
					}
					if (d2 < r - ERR_VALUE)
					{
						vector<Type> pts = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, p3, p1);
						if (pts.size() == 2)
						{
							d_intsect_pts.insert(make_pair(d2, pts));
						}
					}
					if (d3 < r - ERR_VALUE)
					{
						vector<Type> pts = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, p1, p2);
						if (pts.size() == 2)
						{
							d_intsect_pts.insert(make_pair(d3, pts));
						}
					}

					if (d_intsect_pts.size() > 0)
					{
						typename map<double, vector<Type> >::iterator iter = d_intsect_pts.begin();
						Type pt4 = (*iter).second[0];
						Type pt5 = (*iter).second[1];
						area = Calculate_Area_Arch(c, pt4, pt5);
						center = Calculate_ShapeCenter_Arch(c, pt4, pt5);
					}
					if (d_intsect_pts.size() > 1)
					{
						typename map<double, vector<Type> >::iterator iter = d_intsect_pts.begin();
						iter++;
						Type pt6 = (*iter).second[0];
						Type pt7 = (*iter).second[1];

						double old_area = area;
						Type old_sc = center;

						double area_arch_67 = Calculate_Area_Arch(c, pt6, pt7);
						Type center_arch_67 = Calculate_ShapeCenter_Arch(c, pt6, pt7);

						area = old_area - area_arch_67;
						center = (old_area * old_sc - area_arch_67 * center_arch_67) / area;
					}
					if (d_intsect_pts.size() > 2)
					{
						typename map<double, vector<Type> >::iterator iter = d_intsect_pts.begin();
						iter++;
						iter++;

						Type pt8 = (*iter).second[0];
						Type pt9 = (*iter).second[1];

						double old_area = area;
						Type old_sc = center;

						double area_arch_89 = Calculate_Area_Arch(c, pt8, pt9);
						Type center_arch_89 = Calculate_ShapeCenter_Arch(c, pt8, pt9);

						area = old_area - area_arch_89;
						center = (old_area * old_sc - area_arch_89 * center_arch_89) / area;
					}
				}
				else if ((int)(pts_in_circle.size()) == 1)
				{
					//	pt1在圆内部，pt4和pt5为线段p1p2和p1p3与圆的交点				
					const Type& pt1 = pts_in_circle[0];
					Type normal = Calculate_Normal(pt1, pts_out_circle[0], pts_out_circle[1]);
					Type pt2, pt3;
					//	确保P1 -> P2 -> P3的法线为正向
					if (normal.z() > 0)
					{
						pt2 = pts_out_circle[0];
						pt3 = pts_out_circle[1];
					}
					else
					{
						pt2 = pts_out_circle[1];
						pt3 = pts_out_circle[0];
					}

					Type pt4 = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt1, pt2)[0];
					Type pt5 = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt1, pt3)[0];

					double area_arch_45 = Calculate_Area_Arch(c, pt4, pt5);
					double area_triangle_145 = Calculate_Area_Triangle(pt1, pt4, pt5);

					Type center_arch_45 = Calculate_ShapeCenter_Arch(c, pt4, pt5);
					Type center_triangle_145 = Calculate_ShapeCenter_Triangle(pt1, pt4, pt5);

					area = area_arch_45 + area_triangle_145;
					center = (area_arch_45 * center_arch_45 + area_triangle_145 * center_triangle_145) / area;

					//	如果线段23与圆相交，则需要去掉弓形67
					double d = Distance_pt2line(c, pt2, pt3);
					if (d < r - ERR_VALUE)
					{
						vector<Type> pts = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt2, pt3);
						if (pts.size() == 2)
						{
							Type pt6 = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt2, pt3)[0];
							Type pt7 = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt2, pt3)[1];

							double area_arch_67 = Calculate_Area_Arch(c, pt6, pt7);
							Type center_arch_67 = Calculate_ShapeCenter_Arch(c, pt6, pt7);

							double area_section_145 = area;
							Type center_section_145 = center;

							area = area_section_145 - area_arch_67;
							center = (area_section_145 * center_section_145 - area_arch_67 * center_arch_67) / area;
						}
					}
				}
				else if ((int)(pts_in_circle.size()) == 2)
				{
					const Type& pt1 = pts_out_circle[0];
					Type normal = Calculate_Normal(pts_in_circle[0], pt1, pts_in_circle[1]);
					Type pt2;
					Type pt3;
					//	确保P2 -> P1 -> P3的法线为正向
					if (normal.z() > 0)
					{
						pt2 = pts_in_circle[0];
						pt3 = pts_in_circle[1];
					}
					else
					{
						pt2 = pts_in_circle[1];
						pt3 = pts_in_circle[0];
					}

					Type pt4 = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt1, pt2)[0];
					Type pt5 = Calculate_Intersect_Point_Cricle_Segment<Type>(c, r, pt1, pt3)[0];

					double area_triangle_245 = Calculate_Area_Triangle(pt2, pt4, pt5);
					double area_triangle_235 = Calculate_Area_Triangle(pt2, pt3, pt5);
					double area_arch_45 = Calculate_Area_Arch(c, pt4, pt5);

					Type center_triangle_235 = Calculate_ShapeCenter_Triangle(pt2, pt3, pt5);
					Type center_triangle_245 = Calculate_ShapeCenter_Triangle(pt2, pt4, pt5);
					Type center_arch_45 = Calculate_ShapeCenter_Arch(c, pt4, pt5);

					area = area_arch_45 + area_triangle_235 + area_triangle_245;
					center = (area_arch_45 * center_arch_45 + area_triangle_235 * center_triangle_235 + area_triangle_245 * center_triangle_245) / area;
				}
				else if ((int)(pts_in_circle.size()) == 3)
				{
					area = Calculate_Area_Triangle(p1, p2, p3);
					center = Calculate_ShapeCenter_Triangle(p1, p2, p3);
				}

				return area;
			}
			//	计算圆与四边形相交的面积
			template<typename Type>
			double Calculate_Intersect_Area_Circle_Quadrangle(const Type& c, double r, const Type& p1, const Type& p2, const Type& p3, const Type& p4, Type& center)
			{
				Type sc1, sc2;
				double a1 = Calculate_Intersect_Area_Circle_Triangle(c, r, p1, p2, p3, sc1);
				double a2 = Calculate_Intersect_Area_Circle_Triangle(c, r, p1, p3, p4, sc2);
				double area = a1 + a2;
				if (area < ERR_VALUE)
				{
					area = 0.0;
					center = (p1 + p2 + p3 + p4) / 4;
				}
				else
				{
					center = (a1 * sc1 + a2 * sc2) / area;
				}

				return area;
			}
		} //	end of namespace PERIDYNAMIC
	} // end of namespace SAE
} // end of namespace DLUT

#endif