

#pragma once

namespace gismo
{



template <typename T>
class gsSpringPatch
{
    /**
       
       boundary: 0, 1, .., parDim-1
     */
    gsSpringPatch( const gsMultiPatch<T> & boundary)
    : m_boundary(&boundary)
    {
        GISMO_ASSERT( m_boundary->numCurves() == 2*m_boundary->parDim(), 
                      "Expecting "<<2*m_boundary->parDim()<<" boundaries");
    }

    ~gsSpringPatch()
    {
        delete m_result;
    }

    void compute();

    const gsGeometry<T> & result() { return m_result; }

private:

    const gsMultiPatch<T> * m_boundary;
        
    gsGeometry<T> * m_result;

}; // gsSpringPatch


///////////////////////////////////////
///////////////////////////////////////

template <typename T>
class gsSpringPatch<T>::compute()
{
    // check codim == 1
    
    // Construct tensor-basis 
    //switch (dim)
    // switch derived type..?

    // init coefs with boundary
    
    const int dd = 2 * ( m_boundary[0]->parDim()+1 );
    // get strides
    
    // Init/reserve sparse mat and rhs
    
    // * use the cube of interior points
    //nextCubePoint
    
    // For all interior nodes // rows
    //     A(i,i) = dd;
    //  For all neighbors (strides)
    //      if interior
    //       A(i,j) = -1;
    //      else
    //       b(i) += coefs;
    
    // sol = 
    
    // For all interior nodes
    //  basis.index( )
    //    write out sol to coefs
    
    // Construct geometry
    
}





}// namespace gismo
