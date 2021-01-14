/** @file gsRMShellBase.hpp

	@brief Provides basic parameters for the RM shell equation.

	This file is part of the G+Smo library.

	This Source Code Form is subject to the terms of the Mozilla Public
	License, v. 2.0. If a copy of the MPL was not distributed with this
	file, You can obtain one at http://mozilla.org/MPL/2.0/.

	Author(s): A. Mantzaflaris, Y. Xia, HS. Wang

	Date:   2020-12-23
*/
#pragma once

#include <gismo.h>
#include <string.h>
#include"gsMyBase/gsMyBase.h"
using namespace Eigen;
using namespace std;

namespace gismo
{
	// 材料属性
	class Material
	{
	public:
		// 默认构造函数
		//Material() {}
		// 构造函数（初始化）
		Material()
		{
			thickness	= 0.01; // 壳的厚度(m)
			E_modulus	= 200E9;// 弹性模量(Pa)
			poisson_ratio = 0.3;  // 泊松比
			rho			= 1.0;  // 密度
			lambda		= 0.0;
			mu			= 0.0;
		}
		// 析构函数
		~Material() {}
	public:
		real_t E_modulus;
		real_t poisson_ratio;
		real_t rho;
		real_t thickness;
		real_t lambda;
		real_t mu;
	};

	// ======================================================



	// =============================================================
	struct gsShellBoundaryCondition
	{
		int patch;
		gsVector<double> point;
		bool is_parametric;
		double value;
		int direction;  // UX=0, UY=1, UZ=2, RX=3, RY=4, RZ=5, ALL=6.

		gsShellBoundaryCondition() {}
		// position is single on double parametric 
		gsShellBoundaryCondition(const int pat,
			const double posU1,
			const double posU2)
		{
			patch = pat;
			is_parametric = true;
			point.resize(2);
			point(0) = posU1;
			point(1) = posU2;
		}
		
		~gsShellBoundaryCondition() {}
	};
	// ===================================================================

	/** @brief
	Struct defining a point together with a scalar or vector load.

	\ingroup Pde
*/
	template<class T>
	struct point_load
	{
		point_load(const gsVector<T>& _point,
			const T             _value,
			int _patch = 0,
			bool _parametric = true)
			:
			patch(_patch), value(_value), point(1), parametric(_parametric)
		{
			point[0] = _value;
		}

		point_load(const gsVector<T>& _point,
			const gsVector<T>& _value,
			int _patch = 0,
			bool _parametric = true)
			:
			patch(_patch), value(_value), point(_point), parametric(_parametric)
		{ }

		int patch;

		gsVector<T> value;

		gsVector<T> point;

		bool parametric;
	};

	/** @brief Class containing a set of points on a multi-patch
		isogeometric domain, together with boundary conditions.

		\ingroup Pde
	*/
	template<class T>
	class gsPointLoads
	{
	public:
		typedef point_load<T> pLoad;
		typedef typename std::vector<pLoad> plContainer;

		typedef typename std::vector<pLoad>::iterator iterator;

		typedef typename std::vector<pLoad>::const_iterator const_iterator;

	public:

		/// Prints the object as a string.
		std::ostream& print(std::ostream& os) const
		{
			os << "gsPointLoads: " << m_pointLoads.size() << "\n";
			return os;
		}

	public:

		/// Default empty constructor
		gsPointLoads()
		{ }

		~gsPointLoads() // Destructor
		{ }

	public:

		void clear()
		{
			m_pointLoads.clear();
		}

		inline pLoad   operator [] (size_t i) const { return m_pointLoads[i]; }
		inline pLoad& operator [] (size_t i) { return m_pointLoads[i]; }

		void addLoad(const gsVector<T>& _point,
			const gsVector<T>& _value,
			int _patch = 0,
			bool _parametric = true)
		{
			m_pointLoads.push_back(pLoad(_point, _value, _patch, _parametric));
		}

		void addLoad(const gsVector<T>& _point,
			const T             _value,
			int _patch = 0,
			bool _parametric = true)
		{
			m_pointLoads.push_back(pLoad(_point, _value, _patch, _parametric));
		}

		size_t numLoads() const { return  m_pointLoads.size(); }

	private:

		plContainer  m_pointLoads; ///< List of Point loads

	}; // class gsPointLoads
	// =============================================================================

/** @brief
	Struct defining a pressure load with a scalar or vector load.

	\ingroup Pde
*/
	template<class T>
	struct distri_load
	{
		distri_load(const gsVector<T>& _value,
			int _patch = 0,
			string _zone = "",
			bool _parametric = true)
			:
			patch(_patch), value(_value), section(_zone),parametric(_parametric)
		{
		}

		int patch;

		gsVector<T> value;
		string section;
		bool parametric;
	};



	/** @brief Class containing a set of distributed loads on a multi-patch
		isogeometric domain, together with boundary conditions.

		\ingroup Pde
	*/
	template<class T>
	class gsDistriLoads
	{
	public:
		typedef distri_load<T> pLoad;
		typedef typename std::vector<pLoad> plContainer;

		typedef typename std::vector<pLoad>::iterator iterator;

		typedef typename std::vector<pLoad>::const_iterator const_iterator;

	public:

		/// Prints the object as a string.
		std::ostream& print(std::ostream& os) const
		{
			os << "gsDistributedLoads: " << m_distriLoads.size() << "\n";
			return os;
		}

	public:

		/// Default empty constructor
		gsDistriLoads()
		{ }

		~gsDistriLoads() // Destructor
		{ }

	public:

		void clear()
		{
			m_distriLoads.clear();
		}

		inline pLoad   operator [] (size_t i) const { return m_distriLoads[i]; }
		inline pLoad& operator [] (size_t i) { return m_distriLoads[i]; }

		void addLoad(const gsVector<T>& _value,
			int _patch = 0,
			string _zone = "",
			bool _parametric = true)
		{
			m_distriLoads.push_back(pLoad(_value, _patch, _zone, _parametric));
		}

		size_t numLoads() const { return  m_distriLoads.size(); }

	public:

		plContainer  m_distriLoads; ///< List of domain loads

	}; // class gsDistriLoads
}
