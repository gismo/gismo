/** @file monitor_poisson_composed_r-adaptivity.cpp

	@brief Tutorial on how to use expression assembler to solve the Poisson equation

	This file is part of the G+Smo library.

	This Source Code Form is subject to the terms of the Mozilla Public
	License, v. 2.0. If a copy of the MPL was not distributed with this
	file, You can obtain one at http://mozilla.org/MPL/2.0/.

	Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsModeling/gsBarrierCore.h>
#include <gsOptimizer/gsGradientDescent.h>
#include <gsOptim/gsOptim.h>

//! [Include namespace]

namespace gismo{
namespace expr{

template<class T>
class detJ_expr : public _expr<detJ_expr<T> >
{
public:

	typedef T Scalar;

private:
	typename gsGeometryMap<T>::Nested_t _G;

	const short_t domainDim;
	const short_t targetDim;

	mutable gsMatrix<Scalar> jac;
	mutable gsMatrix<Scalar> S;
	mutable gsMatrix<Scalar> U_VT;
	mutable Scalar detG;

public:
	enum{ Space = 0, ScalarValued= 1, ColBlocks= 0};

	detJ_expr(const gsGeometryMap<Scalar> & G) 
	: 
	_G(G),
	domainDim(_G.source().domainDim()),
	targetDim(_G.source().targetDim())
	{
	}

#   define Eigen gsEigen
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

	const Scalar & eval(const index_t k) const
	{
		// Jacobian of G (domainDim x targetDim)
        jac = _G.data().values[1].reshapeCol(k, domainDim, targetDim).transpose();
		// SVD gives:
		// - U: targetDim x targetDim
		// - S: targetDim x domainDim. svd.SingularValues() is min(domainDim,targetDim) x 1
		// - V: domainDim x domainDim
		gsEigen::JacobiSVD<gsEigen::MatrixXd> svd(jac, gsEigen::ComputeFullU | gsEigen::ComputeFullV);
		// Compute the Jacobian determinant using the singular value decomposition
		S.setZero(targetDim,domainDim);
		S.diagonal() = svd.singularValues();

		U_VT = svd.matrixU().topRows(domainDim) * S * svd.matrixV().transpose();
		// Compute the determinant
		detG = U_VT.determinant();
		return detG;
	}

	index_t rows() const { return 0; }

	index_t cols() const { return 0; }

	void parse(gsExprHelper<Scalar> & evList) const
	{
		evList.add(_G);
        _G.data().flags |= NEED_DERIV;
		// jac_expr<Scalar>(_G).parse(evList);
		// _G.parse(evList);
	}


	const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
	const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

	void print(std::ostream &os) const { os << "detJ("; _G.print(os); os <<")"; }


};

template<class T> EIGEN_STRONG_INLINE
detJ_expr<T> detJ(const gsGeometryMap<T> & G) { return detJ_expr<T>(G); }

template<class T>
class penJinv_expr : public _expr<penJinv_expr<T> >
{
public:

	typedef T Scalar;

private:
	typename gsGeometryMap<T>::Nested_t _G;

	const Scalar eps;

	const short_t domainDim;
	const short_t targetDim;

	mutable gsMatrix<Scalar> jac;
	mutable gsMatrix<Scalar> S;
	mutable gsMatrix<Scalar> U_VT;
	mutable Scalar chi;
	mutable Scalar detG;
	mutable gsMatrix<Scalar> res;

public:
	enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

	penJinv_expr(const gsGeometryMap<Scalar> & G, const Scalar epsilon) 
	: 
	_G(G),
	eps(epsilon),
	domainDim(_G.source().domainDim()),
	targetDim(_G.source().targetDim())
	{
	}

// #   define Eigen gsEigen
// 	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
// #   undef Eigen

	const gsMatrix<Scalar> & eval(const index_t k) const
	{
		// Jacobian of G (targetDim x domainDim)
        jac = _G.data().values[1].reshapeCol(k, domainDim, targetDim).transpose();

		// SVD gives:
		// - U: targetDim x targetDim
		// - S: targetDim x domainDim. svd.SingularValues() is min(domainDim,targetDim) x 1
		// - V: domainDim x domainDim
		gsEigen::JacobiSVD<gsEigen::MatrixXd> svd(jac, gsEigen::ComputeFullU | gsEigen::ComputeFullV);
		// Compute the Jacobian determinant using the singular value decomposition
		S.setZero(targetDim,domainDim);
		S.diagonal() = svd.singularValues();

		U_VT = svd.matrixU().topRows(domainDim) * S * svd.matrixV().transpose();
		// Compute the determinant
		detG = U_VT.determinant();
		// Compute the penalized determinant
		chi = (detG<0) ? 0.5 * (detG + math::pow(eps + math::pow(detG, 2.0), 0.5)) : 1.;

		// Compute the inverse of the singular values
		S.diagonal().array() = 1.0 / S.diagonal().array(); 
		S.transposeInPlace(); // From now on, S is its inverse

		// Compute the penalized inverse equivalent of the Moore-Penrose inverse		
		res = svd.matrixV() * S * svd.matrixU().transpose() / chi;
		return res;
	}

	index_t rows() const { return domainDim; }

	index_t cols() const { return domainDim; }

	void parse(gsExprHelper<Scalar> & evList) const
	{
		evList.add(_G);
        _G.data().flags |= NEED_DERIV;
		// jac_expr<Scalar>(_G).parse(evList);
		// _G.parse(evList);
	}

	const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
	const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

	void print(std::ostream &os) const { os << "penJinv("; _G.print(os); os <<")"; }


};

template<class T> EIGEN_STRONG_INLINE
penJinv_expr<T> penJinv(const gsGeometryMap<T> & G, const T eps) { return penJinv_expr<T>(G,eps); }

}
}

using namespace gismo;

// template<typename T = real_t>
// class gsExprAsFunction : public gsFunction<T>
// {
// public:
//     template<class E>
//     gsExprAsFunction( const expr::_expr<E> & expr)
//     : m_expr(expr)
//     {}

//     short_t domainDim() const { return m_expr.parDim(); }
//     short_t targetDim() const { return m_expr.rows();   }

//     void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
//     {
//         result.resize(u,this->targetDim());
//         gsExprEvaluator<T> ev;
//         for (index_t k = 0; k!=u.cols(); k++)
//             result.col(k) = ev.eval(m_expr,u.col(k));
//     }

// protected:
//     const expr::_expr<E> & m_expr;
// };




template<typename T = real_t>
class gsOptMesh : public gsOptProblem<T> 
{
	using Base = gsOptProblem<T>;

private:
	typedef typename gsExprAssembler<T>::geometryMap geometryMap;
	typedef typename gsExprAssembler<T>::space space;
	typedef typename gsExprAssembler<T>::solution solution;

public:
	gsOptMesh(  gsFunction<T> & composition,
				const gsGeometry<T> & geometry,
				T eps = 1e-4)
	:
	m_comp(&composition),
	m_geom(geometry),
	m_eps(eps)
	{
		m_numDesignVars = m_comp->nControls();
		m_curDesign.resize(m_numDesignVars,1);
	}

	/// Evaluates the objective function at the given point u.
	T evalObj(const gsAsConstVector<T> &u) const override
	{
		for (index_t k=0; k!=u.rows(); k++)
			m_comp->control(k) = u[k];

		gsComposedGeometry<T> cgeom(*m_comp,m_geom);

		gsMultiPatch<> mp;
		mp.addPatch(cgeom);
		gsMultiBasis<> mb(m_geom.basis());

		m_evaluator.setIntegrationElements(mb);
		geometryMap G = m_evaluator.getMap(mp);

		if (m_geom.domainDim()==m_geom.targetDim())
		{
			auto detG = jac(G).det();
			auto M = 1./detG;
			auto chi = 0.5 * (detG + pow(m_eps + pow(detG, 2.0), 0.5));
			return m_evaluator.integral( (M*jac(G).adj()/chi).sqNorm()*meas(G));
		}
		else
		{
			auto detG = detJ(G);
			auto M = 1./detG;
			return m_evaluator.integral( (M*penJinv(G,m_eps)).sqNorm()*meas(G));
		}
		// else 
		// {
		//     auto jacG = signed svd;
		//     /* 
		//         if (sigma2>0) // smallest one
		//             chi = sigma1*sigma2 // = jac.det
		//         else
		//             jacG = sigma1*sigma2
		//             chi = 0.5 * (jacG + pow(eps.val() + pow(jacG, 2.0), 0.5));

		//      */
		//     auto chi = 0.5 * (jacG + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
		// }
		
	}

	// /// Computes the gradient of the objective function at the given point u
	// // and stores it in result.
	// void gradObj_into(const gsAsConstVector<T> &u,
	//                 gsAsVector<T> &result) const override;

protected:

	// From gsOptProblem
	using Base::m_curDesign;
	using Base::m_numDesignVars;


	gsFunction<T> * m_comp;
	const gsGeometry<T> & m_geom;
	T m_eps;

	mutable gsExprEvaluator<T> m_evaluator;
};


int main(int arg, char *argv[])
{
	//! [Parse command line]
	bool plot = false;
	index_t numRefine  = 2;
	index_t numRefineD = 2;
	index_t numElevate = 0;
	index_t numElevateD= 0;
	index_t maxIt = 100;
	real_t tol_g = 5e-5;
	real_t eps = 1e-4;
	bool slide = true;
	index_t testCase = 0;
	index_t opt = 2;
	std::string input = "domain2d/lake.xml";

	gsCmdLine cmd("Tutorial on solving a Poisson problem.");
	cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
	cmd.addInt( "E", "elevDomain","Number of degree elevation steps to perform for the domain", numElevateD );
	cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
	cmd.addInt( "R", "refDomain", "Number of Uniform h-refinement loops for the domain",  numRefineD );
	cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
	cmd.addReal("g", "tolG", "relative tol", tol_g);
	cmd.addInt( "i", "maxIt", "max num iterations",  maxIt );
	cmd.addReal("", "eps", "eps",  eps );
	cmd.addSwitch("noslide", "Do not slide the boundaries",  slide );
	cmd.addInt( "o", "opt", "Optimizer: 0: gsGradientDescent, 1: gsHLBFGS, 2: gsOptim::LBFGS.",  opt );
	cmd.addString("f", "file", "Input file", input);

	try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
	//! [Parse command line]

	// Load XML file - multi-patch computational domain
	//! [Read geometry]
	// Check if the input file exists
	if (!gsFileManager::fileExists(input)) 
	{
		gsWarn << "The file cannot be found!\n";
		return EXIT_FAILURE;
	}

	// MultiPatch reader
	gsInfo << "Read file \"" << input << "\"\n";
	gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(input);
	gsInfo << " Got" << *mp << " \n";
	//! [Read geometry]

	mp->embed(3);

	gsMultiBasis<> mb(*mp);

	// Basis for the square domain
	gsKnotVector<> kv2({0,0,1,1},1);
	gsTensorBSplineBasis<2> dbasis(kv2,kv2);
	dbasis.degreeElevate(numElevateD);
	for (index_t i = 0; i < numRefineD; i++)
		dbasis.uniformRefine();
	
	gsInfo<<"Mapper basis:\n"<<dbasis<<"\n";

	gsSquareDomain<2,real_t> domain(dbasis);
	domain.options().addSwitch("Slide","",slide);
	domain.applyOptions();

	gsComposedGeometry<real_t> cspline(domain,mp->patch(0));


	gsOptMesh<> optMesh(domain,mp->patch(0),eps);
	gsVector<> controls(domain.nControls());
	 for (size_t k=0; k!=domain.nControls(); k++)
		controls[k] = domain.control(k);


	gsOptimizer<real_t> * optimizer;
	if      (opt==0) // gsGradientDescent
	{
		optimizer = new gsGradientDescent<real_t>(&optMesh);
		optimizer->options().setInt("MaxIterations",maxIt);
		optimizer->options().setInt("Verbose",2);
		optimizer->options().setReal("MinGradientLength",tol_g);

	}
	else if (opt==1) // gsHLBFGS
	{
		optimizer = new gsHLBFGS<real_t>(&optMesh);
		optimizer->options().setInt("MaxIterations",maxIt);
		optimizer->options().setInt("Verbose",2);
		optimizer->options().setReal("tolRelG",tol_g);
	}
	else if (opt==2) //gsOptim::LBFGS
	{
		optimizer = new gsOptim<real_t>::LBFGS(&optMesh);
		optimizer->options().setInt("MaxIterations",maxIt);
		optimizer->options().setInt("Verbose",1);
		optimizer->options().setReal("GradErrTol",tol_g);
	}
	else
	{
		GISMO_ERROR("Unknown optimizer");
	}

	optimizer->solve(controls);
	gsVector<> optSol = optimizer->currentDesign();
	
	for (size_t k=0; k!=optSol.rows(); k++)
		domain.control(k) = optSol[k];

    gsMultiPatch<> cmp;
    cmp.addPatch(cspline);
    // mp.embed(3);
    gsMultiBasis<> cmb(cmp);


	gsExprEvaluator<> ev;
	ev.setIntegrationElements(cmb);
	auto Gold = ev.getMap(*mp);
	auto Gnew = ev.getMap(cmp);

	gsWriteParaview(domain.domain(),"domain",1000,true,true);
	// ev.writeParaview(detJ(G),G,"jacobian_determinant");
	ev.writeParaview(detJ(Gold),Gold,"OldJacobian_determinant");
	ev.writeParaview(detJ(Gnew),Gnew,"NewJacobian_determinant");

	delete optimizer;
	return 0;
}// end main
