#include <gismo.h>
#include <Eigen/Dense>
#include <gsModeling/gsFitting.h>
#include <gsModeling/gsFitting.hpp>


namespace gismo {

template <class T>
class gsFastFitting : public gsFitting<T>
{
public:
    gsFastFitting(gsMatrix<T> const & param_values,
                    gsMatrix<T> const & points,
                    gsBasis<T> & basis)
        : gsFitting<T>(param_values, points, basis)
    {}
    void changeParam(const gsMatrix<T> &uv, const gsMatrix<T> &xyz);
    void FindGridPoint1D(const gsMatrix<T> &ugrid, const T& curr_param, index_t& k2) const;
    void FindGridPoint(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, const gsMatrix<T>& curr_param, index_t &k1, index_t &k2);
    void GridProjection(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, gsMatrix<index_t>& weights);
    void BuildLookupTable(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, gsMatrix<T>& Table, const gsMatrix<index_t>& weights);
    void assembleSystem(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, gsMatrix<T>& Table, gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B);
    void compute(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, const bool condcheck);
    void computeAllProjectedErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz, const gsMatrix<T>& ugrid, const gsMatrix<T>& vgrid);
    void computeProjectedAverageErrors(const gsMatrix<index_t>& weights,const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid);
    void plotErrors(const std::string & fname) const;
};


template<class T>
void gsFastFitting<T>::changeParam(const gsMatrix<T> &uv, const gsMatrix<T> &xyz)
{
    GISMO_ASSERT( uv.cols() == xyz.cols() && uv.rows() == 2 && xyz.rows() == 3,
                  "Wrong input");
    this -> m_param_values = uv;
    this -> m_points = xyz.transpose();
}


template<class T>
void gsFastFitting<T>::FindGridPoint1D(const gsMatrix<T> &ugrid, const T& curr_param, index_t& k2) const
{
    // Could be improved with a faster search method
    if ((ugrid(0)+ugrid(1))/2 >= curr_param)
    {
        k2 = 0;
        return;
    }

    for (index_t i= 1; i< ugrid.cols()-1 ; i++)
    {
        if ((ugrid(i-1)+ugrid(i))/2< curr_param && (ugrid(i)+ugrid(i+1))/2 >= curr_param)
        {
            k2 = i;
            return;
        }
    }

    k2 = ugrid.cols()-1;
}

template<class T>
void gsFastFitting<T>::FindGridPoint(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, const gsMatrix<T>& curr_param, index_t& k1, index_t& k2)
{
    // Search function, k is the global index of the nearest neighbour of curr_param, with k=(k1,k2)
    FindGridPoint1D(ugrid,curr_param(0),k1);
    FindGridPoint1D(vgrid,curr_param(1),k2);
}

template<class T>
void gsFastFitting<T>::GridProjection(const gsMatrix<T> &ugrid, const gsMatrix<T> &vgrid, gsMatrix<index_t>& weights)
{
    // For now the identity:
    index_t n1,n2;
    n1 = weights.rows();
    n2 = weights.cols();
    weights.setZero(n1,n2);

    const int dimension = this->m_points.cols();
    const int numpoints = this->m_points.rows();

    // Initialize new m_params, m_points
    gsMatrix<T> params_temp, points_temp;
    params_temp.setZero(2,n1*n2);
    points_temp.setZero(n1*n2,dimension);
    gsMatrix<T> curr_param,curr_point;

    //gsTensorBSplineBasis<2,T> *tensorbasis = (static_cast<gsTensorBSplineBasis<2,T>*>(this->m_basis));
    for (index_t i=0; i<numpoints; i++)
    {
        // we have to check which points belong to which grid point. New variables:
        // m_params: ugrid x vgrid as tensor product structure -- Done in the next step
        // m_points: init with 0, add p_i to corresponding col i, where gridpoint (u_i1,v_i2) (of ugrid, vgrid) are the nearest neighbour of m_params
        // weights: increase weights(i) by one for each m_points added to i

        curr_param = this->m_param_values.col(i);
        curr_point = this->m_points.row(i);
        index_t k,k1=0,k2=0;
        FindGridPoint(ugrid,vgrid,curr_param,k1,k2);
        k = k2 * n1 + k1;
        //k = tensorbasis->index(k1,k2);

        points_temp.row(k) += curr_point;
        weights(k1,k2) += 1;
    }
    // For now, adjust m_param_values:+= curr_point;
    for (index_t i2= 0; i2<n2; i2++)
    {
        for (index_t i1= 0; i1<n1; i1++)
        {
            //index_t iglobal = tensorbasis->index(i1,i2);
            index_t iglobal = i2 * n1 + i1;
            params_temp(0,iglobal)=ugrid(i1);
            params_temp(1,iglobal)=vgrid(i2);
        }
    }

    this -> m_param_values = params_temp;
    this -> m_points = points_temp;
}

template<class T>
void gsFastFitting<T>::BuildLookupTable(const gsMatrix<T>& ugrid, const gsMatrix<T>& vgrid, gsMatrix<T>& Table, const gsMatrix<index_t>& weights)
{
    //Initialize look-up Table
    //const int num_basis=this->m_basis->size();
    //const int dimension=this->m_points.cols();

    // Get 1D basis:
    gsBasis<T>* ubasis_tmp;
    ubasis_tmp = &(this->m_basis->component(0));
    //vbasis = this->m_basis->component(1);
    gsBSplineBasis<T> ubasis = *(static_cast<gsBSplineBasis<T>*>(ubasis_tmp));

    index_t N1;
    //index_t N2;
    N1 = ubasis.size();
    //N2 = num_basis/N1;

    int p;
    p = this->m_basis->degree(0);
    int n1,n2;
    n1 = vgrid.cols();  // gridpoints in u-direction
    n2 = ugrid.cols();  // gridpoints in v-direction
    Table.setZero(N1*N1,n2);

    // pre-evaluate ubasis in grid-points
    gsMatrix<unsigned> actives;
    gsMatrix<real_t> uvalues;
    ubasis.active_into(ugrid,actives);
    ubasis.eval_into(ugrid,uvalues);

    // Build look-up table, i.e. elements h_{k2,i2,j2} in Table
    index_t flops = 0;
    for (index_t k1=0; k1<n1; k1++)
    {
        for (index_t i=0; i<p+1; i++)
        {
            for (index_t j=0; j<p+1; j++)
            {
                real_t bi,bj;
                index_t i1, j1;
                bi = uvalues(i,k1);
                i1 = actives(i,k1);
                bj = uvalues(j,k1);
                j1 = actives(j,k1);
                index_t ind_temp = j1*N1+i1;   // lexicographical running index
                for (index_t k2=0; k2<n2; k2++)
                {
                    if(weights(k1,k2)!=0)
                    {
                        // Judge i < j and use symmetry of bi*bj -> more efficiently?
                        Table(ind_temp,k2) += bi*bj*weights(k1,k2);
                        flops += 3;
                    }
                }
            }
        }
    }
    // Output number flops:
    gsInfo << "Counted flops fast fitting--PART 1: " << flops << std::endl;
}

template<class T>
void gsFastFitting<T>::assembleSystem(const gsMatrix<T>& ugrid, const gsMatrix<T>& vgrid, gsMatrix<T>& Table,gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B)
{
    //Initialize look-up Table
    const int num_basis=this->m_basis->size();
    //const int dimension=this->m_points.cols();

    // Get 1D basis:
    gsBasis<T>* vbasis_tmp;
    vbasis_tmp = &(this->m_basis->component(1));
    //ubasis = this->m_basis->component(0);
    gsBSplineBasis<T> vbasis = *(static_cast<gsBSplineBasis<T>*>(vbasis_tmp));

    index_t N1;
    index_t N2;
    N2 = vbasis.size();
    N1 = num_basis/N2;

    int p;
    p = this->m_basis->degree(1);
    int n1;
    n1 = vgrid.cols();  // gridpoints in u-direction
    int n2;
    n2 = ugrid.cols();  // gridpoints in v-direction
    A_mat.setZero();

    // pre-evaluate ubasis in grid-points
    gsMatrix<unsigned> actives;
    gsMatrix<real_t> vvalues;
    vbasis.active_into(vgrid,actives);
    vbasis.eval_into(vgrid,vvalues);

    gsTensorBSplineBasis<2,T> *tensorbasis = (static_cast<gsTensorBSplineBasis<2,T>*>(this->m_basis));
    // Assemble matix A, right hand side b
    index_t flops = 0;
    for (index_t k2=0; k2<n1; k2++)
    {
        for (index_t i2=0; i2<p+1; i2++)
        {
            for (index_t j2=0; j2<p+1; j2++)
            {
                real_t bi,bj;
                index_t ui, uj;
                bi = vvalues(i2,k2);
                ui = actives(i2,k2);
                bj = vvalues(j2,k2);
                uj = actives(j2,k2);
                for (index_t i1=0; i1<N1; i1++)
                {
                    for (index_t j1=std::max(0,i1-p); j1<std::min(i1+p+1,N1); j1++)
                    {
                        index_t i_lexiglobal = tensorbasis->index(i1,ui);
                        index_t j_lexiglobal = tensorbasis->index(j1,uj);
                        index_t ind_temp = j1*N1+i1;
                        A_mat(i_lexiglobal,j_lexiglobal) += bi*bj*Table(ind_temp,k2);   // Symmetry of bi*bj -> could be done more efficiently?
                        flops += 3;
                    }
                }
            }
        }
    }
    // Info for the flop count:
    gsInfo << "Counted flops fast fitting--PART 2: " << flops << std::endl;

    gsMatrix<T> value, curr_point;
    gsMatrix<unsigned> active;
    for(index_t k = 0; k < n1*n2; ++k)
    {
        curr_point = this->m_param_values.col(k);

        //computing the values of the basis functions at the current point
        this->m_basis->eval_into(curr_point, value);

        // which functions have been computed i.e. which are active
        this->m_basis->active_into(curr_point, active);

        const index_t numActive = active.rows();

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii = active.at(i);
            m_B.row(ii) += value.at(i) * this->m_points.row(k);
        }
    }

    // For testing purpose to compare A_mat:
    gsFileData<> fd;
    gsMatrix<T> C=A_mat.toDense();
    fd << C;
    fd.dump("FastFittingMatrix");
    gsFileData<> fb;
    fb << m_B;
    fb.dump("FastFittingb");
}


template<class T>
void gsFastFitting<T>::compute(const gsMatrix<T>& ugrid, const gsMatrix<T>& vgrid, const bool condcheck)
{
    // Wipe out previous result
    if ( this->m_result )
        delete this->m_result;

    // Initialization of variables num_basis, dimension
    const int num_basis=this->m_basis->size();
    const int dimension=this->m_points.cols();

    //left side matrix
    gsSparseMatrix<T> A_mat(num_basis, num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here, to do: improve
    int nonZerosPerCol = 1;
    for (int i = 0; i < this->m_basis->dim(); ++i)
        nonZerosPerCol *= ( 2 * this->m_basis->degree(i) + 1 ) * 4;
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    m_B.setZero(); // ensure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    // Initialize weights and look-up table
    gsMatrix<index_t> weights (ugrid.cols(),vgrid.cols());
    //weights.setOnes(ugrid.rows(),vgrid.rows());
    gsMatrix<T> Table;
    GridProjection(ugrid,vgrid,weights);
    //weights.setOnes();

    gsStopwatch time;
    time.restart();
    BuildLookupTable(ugrid,vgrid,Table,weights);

    assembleSystem( ugrid, vgrid, Table, A_mat, m_B);
    time.stop();
    gsInfo<<"Assembly time                     : "<< time <<"\n";
    // To do: include regularization later after projection step
    //applySmoothing(lambda, A_mat, dreg);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>:: BiCGSTABILUT solver( A_mat );
    gsMatrix<T> x;
    if (condcheck)
    {
        //real_t cnum;
        //cnum = gsSolverUtils<T>::conditionNumber(A_mat);
        //gsInfo << "Condition number: " << cnum << std::endl;

        //solver.compute(A_mat);
        Eigen::SparseMatrix<double> I(A_mat.rows(),A_mat.cols());
        I.setIdentity();
        auto A_inv = solver.solve(I);
        gsInfo << "Condition number: " << A_mat.norm() * A_inv.norm() << std::endl;
        x = A_inv * m_B;
    }

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        std::cerr<<  "The preconditioner failed. Aborting.";// << std::endl;
        this->m_result = NULL;
        return;
    }

    // Solves for many right hand side  columns
    //else
    //    x = solver.solve(m_B); //toDense()
    x = solver.solve(m_B); //toDense()
    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve
    this->m_result = this->m_basis->makeGeometry( give(x) ).release();
    computeProjectedAverageErrors(weights,ugrid,vgrid);
}


template<class T>
void gsFastFitting<T>::computeAllProjectedErrors(const gsMatrix<T>& uv, const gsMatrix<T>& xyz, const gsMatrix<T>& ugrid, const gsMatrix<T>& vgrid)
{
    this->m_pointErrors.clear();

    gsMatrix<T> uv_temp (2,uv.cols());
    index_t k1,k2;
    for (index_t k=0; k<uv.cols(); k++)
    {
        FindGridPoint(ugrid,vgrid,uv.col(k),k1,k2);
        uv_temp(0,k) = ugrid(k1);
        uv_temp(1,k) = vgrid(k2);

        if(xyz(0,k) == 0 && xyz(1,k)==0 && xyz(2,k)==0)
            gsInfo << "zero at " << k << std::endl;
    }

    changeParam(uv_temp,xyz);
    this -> computeErrors();

    // change back param and points?
}


template<class T>
void gsFastFitting<T>::computeProjectedAverageErrors(const gsMatrix<index_t>& weights, const gsMatrix<T>& ugrid, const gsMatrix<T>& vgrid)
{
    this->m_pointErrors.clear();

    gsMatrix<T> val_i;
    bool c = false;
    //index_t count = 0;
    this->m_result->eval_into(this->m_param_values, val_i);

    index_t k1,k2;
    FindGridPoint(ugrid,vgrid,this->m_param_values.col(0),k1,k2);
    if (weights(k1,k2) != 0)
    {
        this->m_pointErrors.push_back( (this->m_points.row(0)/weights(k1,k2) - val_i.col(0).transpose()).norm() );
        this->m_max_error = this->m_min_error = this->m_pointErrors.back();
        c = true;
        //count += 1;
    }

    for (index_t i = 1; i < this->m_points.rows(); i++)
    {
        FindGridPoint(ugrid,vgrid,this->m_param_values.col(i),k1,k2);
        if (weights(k1,k2) != 0)
        {
            const T err = (this->m_points.row(i)/weights(k1,k2) - val_i.col(i).transpose()).norm() ;
            this->m_pointErrors.push_back(err);

            if (c==false)
            {
                this->m_max_error = this->m_min_error = this->m_pointErrors.back();
                c = true;
            }

            if ( err > this->m_max_error ) this->m_max_error = err;
            if ( err < this->m_min_error ) this->m_min_error = err;
            //count += 1;
        }
    }
}

template<class T>
void gsFastFitting<T>::plotErrors(const std::string & fname) const
{
    // Plot Errors
    const std::vector<real_t>& eval_field = this->pointWiseErrors();
    gsMatrix<real_t> bigMatrix(4, eval_field.size());
    for(size_t i=0; i<eval_field.size(); i++)
    {
        bigMatrix(0,i) = this->m_param_values(0,i);
        bigMatrix(1,i) = this->m_param_values(1,i);
        bigMatrix(2,i) = 0;
        bigMatrix(3,i) = eval_field[i];
    }
    gsWriteParaviewPoints(bigMatrix, fname);
}


} //namespace gismo
