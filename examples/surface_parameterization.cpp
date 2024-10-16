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


namespace gismo {
namespace expr {

template<class E0, class E1, class E2>
class ternary_expr; // ternary expression


template<class E0, class E1, class E2>
class ternary_expr : public _expr<ternary_expr<E0, E1, E2> >
{
  typename E0::Nested_t _u;
  typename E1::Nested_t _v;
  typename E2::Nested_t _w;
 public:
  typedef typename E1::Scalar Scalar;

  explicit ternary_expr(_expr<E0> const &u,
                        _expr<E1> const &v,
                        _expr<E2> const &w)
      :
      _u(u),
      _v(v),
      _w(w) {
    GISMO_ASSERT(E0::ScalarValued, "Condition must be scalar valued");
    GISMO_ASSERT((int) E1::ScalarValued == (int) E2::ScalarValued,
                 "Both v and w must be scalar valued (or not).");
    GISMO_ASSERT((int) E1::ColBlocks == (int) E2::ColBlocks,
                 "Both v and w must be colblocks (or not).");
    GISMO_ASSERT((int) E1::Space == (int) E2::Space,
                 "Both v and w must be space (or not), but E1::Space = "
                     << E1::Space << " and E2::Space = " << E2::Space);
    GISMO_ASSERT(_v.rows() == _w.rows(),
                 "Rows of v and w differ. _v.rows() = " << _v.rows()
                                                        << ", _w.rows() = "
                                                        << _w.rows());
    GISMO_ASSERT(_v.cols() == _w.cols(),
                 "Columns of v and w differ. _v.cols() = " << _v.cols()
                                                           << ", _w.cols() = "
                                                           << _w.cols());
    GISMO_ASSERT(_v.rowVar() == _w.rowVar(), "rowVar of v and w differ.");
    GISMO_ASSERT(_v.colVar() == _w.colVar(), "colVar of v and w differ.");
  }
 public:
  enum {
    ScalarValued = E1::ScalarValued,
    ColBlocks = E1::ColBlocks,
    Space = E1::Space
  }; // == E2::Space

//  const Scalar eval(const index_t k) const { return (_u.eval(k) > 0 ? _v.eval
//  (k) : _w.eval(k)); }

  const Temporary_t eval(const index_t k) const
  {
    return (_u.eval(k) > 0 ? _v.eval(k) : _w.eval(k));
  }

  // { res = eval_impl(_u,_v,_w,k); return  res;}

  index_t rows() const { return _v.rows(); }
  index_t cols() const { return _v.cols(); }
  void parse(gsExprHelper<Scalar> &evList) const {
    _u.parse(evList);
    _v.parse(evList);
    _w.parse(evList);
  }

  const gsFeSpace<Scalar> &rowVar() const { return _v.rowVar(); }
  const gsFeSpace<Scalar> &colVar() const { return _v.colVar(); }

};


/// Ternary ternary_expr
template<class E0, class E1, class E2>
EIGEN_STRONG_INLINE
ternary_expr<E0, E1, E2> ternary(const E0 &u,
                                 const E1 &v,
                                 const E2 &w)
{
  return ternary_expr<E0, E1, E2>(u, v, w);
}

} // namespace expr
}// namespace gismo

using namespace gismo;


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
			// auto chi = 0.5 * (detG + pow(m_eps + pow(detG, 2.0), 0.5));

            // Compute the chi part
            gsConstantFunction<T> epsilon(m_eps, m_comp->domainDim());
            auto eps = m_evaluator.getVariable(epsilon);
            auto chiPPart = eps * ((detG.val() - eps.val()).exp());

            // Ternary operation to compute chi and chip
            auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());

			return m_evaluator.integral( (M*jac(G).adj()/chi).sqNorm()*meas(G));
		}
		else
		{
            auto fform = jac(G).tr()*jac(G);
            auto detG = pow(fform.det().val(),0.5); //jacobian determinant for a surface, i.e. the measure
			auto M = 1./detG;

            // Compute the chi part
            gsConstantFunction<T> epsilon(m_eps, m_comp->domainDim());
            auto eps = m_evaluator.getVariable(epsilon);
            auto chiPPart = eps * ((detG.val() - eps.val()).exp());

            // Ternary operation to compute chi and chip
            auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());

            auto invJacMat = fform.adj()/chi; // inverse of jacobian matrix with 'determinant' replaced
			return m_evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
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

	// mp->embed(3);

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
	ev.writeParaview((jac(Gold).tr()*jac(Gold)).det(),Gold,"OldJacobian_determinant");
	ev.writeParaview((jac(Gnew).tr()*jac(Gnew)).det(),Gnew,"NewJacobian_determinant");

	delete optimizer;
	return 0;
}// end main
