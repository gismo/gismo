/** @file gsHBSplineBasis.h

    @brief Provides implementation of HBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once 

#include <gsHSplines/gsHTensorBasis.h>

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{

/*
template<short_t d, class T>
gsHBSplineBasis<d,T>::gsHBSplineBasis(gsBSplineBasis<T> &  bsbasis)
    : gsHTensorBasis<d,T>( gsTensorBSplineBasis<d,T>(&bsbasis) )
{
    //  Note: The compiler adds automatically a return statement
    //  at the end of each constructor.  Throwing an exception
    //  causes this return statement to be unreachable, and
    //  warning 4702 is emitted.  To stop this warning we add
    //  "bsbasis.dim()==1", which is not known at compile time
    GISMO_ASSERT(d==1, "Wrong dimension");
}
*/

template<short_t d, class T>
typename gsHBSplineBasis<d,T>::BoundaryBasisType * gsHBSplineBasis<d,T>::basisSlice(index_t dir_fixed,T par ) const
{
    GISMO_ASSERT(d-1>=0,"d must be greater or equal than 1");
    GISMO_ASSERT(dir_fixed>=0 && static_cast<index_t>(dir_fixed)<d,"cannot fix a dir greater than dim or smaller than 0");
    const boxSide side(dir_fixed,0);
    const typename gsTensorBSplineBasis<d,T>::BoundaryBasisType::uPtr bBSplineBasis =
            this->m_bases[0]->boundaryBasis(side);
    typename gsHBSplineBasis<d,T>::BoundaryBasisType* bBasis =
            new typename gsHBSplineBasis<d,T>::BoundaryBasisType(*bBSplineBasis);

    if(d!=1)
    {
        std::vector<index_t> boxes;
        this->getBoxesAlongSlice(dir_fixed,par,boxes);
        bBasis->refineElements(boxes);
    }

    return bBasis;
}

template<short_t d, class T>
std::ostream & gsHBSplineBasis<d,T>::print(std::ostream &os) const
{
    os << "Hierarchical B-spline ";
    gsHTensorBasis<d,T>::printBasic(os);
    //this->printCharMatrix(os);
    return os;
}


template<short_t d, class T>
void gsHBSplineBasis<d,T>::initialize()
{
    // Sets everything related to gsHTensorBasis
    // this->update_structure(); // base class update
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    int lvl = this->levelOf(i);
    this->m_bases[lvl]->evalSingle_into( 
        this->m_xmatrix[lvl][ i - this->m_xmatrix_offset[lvl] ], u, result);
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    int lvl = this->levelOf(i);
    this->m_bases[lvl]->derivSingle_into( 
        this->m_xmatrix[lvl][ i - this->m_xmatrix_offset[lvl] ], u, result);
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    int lvl = this->levelOf(i);
    this->m_bases[lvl]->deriv2Single_into( 
        this->m_xmatrix[lvl][ i - this->m_xmatrix_offset[lvl] ], u, result);
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    gsVector<index_t> nact;
    gsMatrix<index_t> act;
    // Get the number of active functions
    this->numActive_into(u, nact);
    // Set everything to zero
    result.setZero(nact.maxCoeff(), u.cols());
    
    gsMatrix<T> curr_res;
    for(index_t j = 0; j < u.cols(); j++) // for all points
    {
        const gsMatrix<T> & curr_u =  u.col(j);
        
        // Compute the indices of active functions on curr_u
        gsHBSplineBasis<d,T>::active_into(curr_u, act);
        
        // for all actives on the point
        for(index_t i = 0; i != act.rows(); i++) 
        {
            gsHBSplineBasis<d,T>::evalSingle_into( act(i,0), curr_u, curr_res );
            result(i,j) = curr_res(0,0);
        }
    }
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    gsVector<index_t> nact;
    gsMatrix<index_t> act;
    // Get the number of active functions
    this->numActive_into(u, nact);
    // Set everything to zero
    result.setZero( d*nact.maxCoeff(), u.cols() );
    
    gsMatrix<T> curr_res;
    for(index_t j = 0; j < u.cols(); j++) // for all points
    {
        const gsMatrix<T> & curr_u =  u.col(j);
        
        // Compute the indices of active functions on curr_u
        gsHBSplineBasis<d,T>::active_into(curr_u, act);
        
        // for all actives on the point
            for(index_t i = 0; i < act.rows(); i++) 
            {
                gsHBSplineBasis<d,T>::derivSingle_into( act(i,0), curr_u, curr_res );
                result.block(i*d,j,d,1).noalias() = curr_res;
            }
    }
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
{
    gsVector<index_t> nact;
    gsMatrix<index_t> act;

    const int blocksize = d*(d + 1)/2;

    // Get the number of active functions
    this->numActive_into(u, nact);
    result.resize( blocksize * nact.maxCoeff(), u.cols() );
    // Set everything to zero
    result.setZero();

    // Current result for u.col(j)
    gsMatrix<T> curr_res;

    for(index_t j = 0; j < u.cols(); j++) // for all points
    {
        const gsMatrix<T> & curr_u =  u.col(j);
        
        // Compute the indices of active functions on curr_u
        gsHBSplineBasis<d,T>::active_into(curr_u, act);
        
        // for all actives at the point
        for(index_t i = 0; i < act.rows(); i++) 
        {
            gsHBSplineBasis<d,T>::deriv2Single_into( act(i,0), curr_u, curr_res );
            result.block(i*blocksize,j,blocksize,1) = curr_res;
        }
    }
}

template<short_t d, class T>
void gsHBSplineBasis<d,T>::transferbyLvl (std::vector<gsSparseMatrix<T> >& result)
{
    result.clear();
    //gsVector<index_t> level;
    //gsMatrix<index_t> b1, b2; //boxes in highest level numbering
    //this->m_tree.getBoxesInLevelIndex(b1, b2, level);//return boxes in level indices
    tensorBasis curTensorlvl = this->tensorLevel(0);
    //std::vector< gsSparseMatrix<T,RowMajor> > transfer(this->maxLevel());
    gsSparseMatrix<T,RowMajor> transfer;
    std::vector<std::vector<T> > knots(d);

    std::vector<CMatrix> xmatLvl_i, xmatLvl_i1;

    this->setActiveToLvl(0, xmatLvl_i );

    for(unsigned i = 1; i <= this->maxLevel(); ++i)
    {
        for(unsigned dim = 0; dim != d; ++dim)
        {
            const gsKnotVector<T> & ckv = m_bases[i-1]->knots(dim);
            const gsKnotVector<T> & fkv = m_bases[i  ]->knots(dim);
            ckv.symDifference(fkv, knots[dim]);
            
            // equivalent (dyadic ref.):
            // ckv.getUniformRefinementKnots(1, knots[dim]);
        }

        curTensorlvl.refine_withTransfer(transfer, knots);

        this->setActiveToLvl(i, xmatLvl_i1);

        const gsSparseMatrix<T> crs = this->coarsening(xmatLvl_i, xmatLvl_i1, transfer);
        result.push_back(crs);

        xmatLvl_i.swap(xmatLvl_i1);
    }
}


template<short_t d, class T>
gsSparseMatrix<T> gsHBSplineBasis<d,T>::coarsening( const std::vector<gsSortedVector<index_t> >& old, const std::vector<gsSortedVector<index_t> >& n, const gsSparseMatrix<T,RowMajor> & transfer) const
{
    int size1= 0;int size2 = 0;
    int glob_numb = 0;//continous numbering of hierarchical basis
    for(unsigned int i =0; i< old.size();i++){//count the number of basis functions in old basis
        size1 += old[i].size();
    }
    for(unsigned int i =0; i< n.size();i++){//count the number of basis functions in new basis
        size2 += n[i].size();
    }
    gsSparseMatrix<T> result(size2,size1);

    //gsMatrix<T> transferDense = transfer;
    gsSparseMatrix<T,ColMajor> temptransfer =  transfer;
//    gsDebug<<"size1 "<<size1<<" size2 "<<size2<<std::endl;
//    gsDebug<<"temptransfer.rows"<< temptransfer.rows()<<" cols "<<temptransfer.cols();
    for (unsigned int i = 0; i < old.size(); i++)//iteration through the levels of the old basis
    {
        // find starting index of level i in new basis
        int start_lv_i = 0;
        for(unsigned int l =0; l < i; l++)
        {
            start_lv_i += n[l].size();
        }

        for (unsigned int j = 0; j < old[i].size();j++)//iteration through the basis functions in the given level
        {
//            gsDebug<<"function "<<j<<std::endl;
            const unsigned old_ij = old[i][j];  // tensor product index

            if( n[i].bContains(old_ij) )//it he basis function was not refined
            {
                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
            }
            else
            {
//          TODO delete if the spafce matrix version is tested
//                for(int k = 0; k < transferDense.rows(); k++)//basis function was refined->looking for the coeficinets from the global transfer matrix
//                {
//                    if(transferDense(k, old_ij) != 0)//if the coefficient is non zero we find the coresponding function in n
//                    {
//                        const int pos = start_lv_i + n[i].size() + std::distance(n[i+1].begin(), n[i+1].find_it_or_fail(k));
//                        result(pos,glob_numb) = transferDense(k, old_ij);
//                    }
//                }
                for(typename gsSparseMatrix<T,ColMajor>::InnerIterator k(temptransfer,old_ij); k; ++k)//basis function was refined->looking for the coeficinets from the global transfer matrix
                {

                    //if(transferDense(k, old_ij) != 0)//if the coefficient is non zero we find the coresponding function in n
                    //{
                    const int pos = start_lv_i + n[i].size() + std::distance(n[i+1].begin(), n[i+1].find_it_or_fail(k.row()));
                    result(pos,glob_numb) = k.value();//transferDense(k, old_ij);
                    //}
                }
            }
            glob_numb++;
        }
    }
    return result;
}

template<short_t d, class T>
gsSparseMatrix<T> gsHBSplineBasis<d,T>::coarsening_direct( const std::vector<gsSortedVector<index_t> >& old, const std::vector<gsSortedVector<index_t> >& n, const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const
{
    int size1= 0;int size2 = 0;
    int glob_numb = 0;//continous numbering of hierarchical basis
    for(unsigned int i =0; i< old.size();i++){//count the number of basis functions in old basis
        size1 += old[i].size();
    }
    for(unsigned int i =0; i< n.size();i++){//count the number of basis functions in new basis
        size2 += n[i].size();
    }
    gsSparseMatrix<T> result(size2,size1);


//    std::vector<gsMatrix<T> > transferDense;// = transfer;
//    transferDense.resize(transfer.size());
//    for (unsigned int i = 0; i < transfer.size();i++){
//        transferDense[i] = transfer[i];
//    }
    std::vector<gsSparseMatrix<T,ColMajor> > temptransfer;// = transfer;
    temptransfer.resize(transfer.size());
    for (unsigned int i = 0; i < transfer.size();i++){
        temptransfer[i] = transfer[i];
    }

    for (unsigned int i = 0; i < old.size(); i++)//iteration through the levels of the old basis
    {
        // find starting index of level i in new basis
        int start_lv_i = 0;
        for(unsigned int l =0; l < i; l++)
        {
            start_lv_i += n[l].size();
        }

        
        for (unsigned int j = 0; j < old[i].size();j++)//iteration through the basis functions in the given level
        {
            
            start_lv_i = 0;
            for(unsigned int l =0; l < i; l++)
            {
                start_lv_i += n[l].size();
            }
            const unsigned old_ij = old[i][j];  // tensor product index

            if( n[i].bContains(old_ij) )//it he basis function was not refined
            {
                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
            }
            else
            {
                std::vector<lvl_coef> coeffs;
                lvl_coef temp;
                temp.pos = old_ij;
                temp.coef = 1;
                temp.lvl = i;
                coeffs.push_back(temp);
                for (size_t coeff_index = 0; coeff_index < coeffs.size(); ++coeff_index)
                {
                    const lvl_coef coeff = coeffs[coeff_index];

                    start_lv_i = 0;
                    for(unsigned int l =0; l < coeff.lvl; l++)
                    {
                        start_lv_i += n[l].size();
                    }

                    for(typename gsSparseMatrix<T,ColMajor>::InnerIterator k(temptransfer[coeff.lvl],coeff.pos); k; ++k)//basis function was refined->looking for the coeficinets from the global transfer matrix
                    {
                        if(n[coeff.lvl+1].bContains(k.row()))
                        {
                            //gsDebug<<"j: "<<j<<" "<<"oldij"<<old_ij<<"k.row() "<<k.row()<<"   ";
                            //gsDebug<<"j: "<<j<<" "<<"globnumb "<<glob_numb<<" "<<" coef "<<  coeff.coef * temptransfer[coeff.lvl](k.row(), coeff.pos)<<" ";
                            const int pos = start_lv_i + n[coeff.lvl].size() + std::distance(n[coeff.lvl+1].begin(), n[coeff.lvl+1].find_it_or_fail(k.row()));
                            //gsDebug<<"pos:"<<pos<<std::endl;
                            //double ppp =  transferDense[coeff.lvl](k, coeff.pos);
                            result(pos,glob_numb) += coeff.coef * temptransfer[coeff.lvl](k.row(), coeff.pos);//transferDense[coeff.lvl](k, coeff.pos);
                        }else
                        {
                            temp.pos = k.row();
                            temp.coef = temptransfer[coeff.lvl](k.row(), coeff.pos) * coeff.coef;
                            temp.lvl = coeff.lvl+1;
                            coeffs.push_back(temp);
                        }
                    }

//                    for(unsigned int ii = 0; ii < coeffs.size();ii++){
//                        gsDebug<<"( "<< coeffs[ii].pos<<" , "<<coeffs[ii].coef<<" , "<<coeffs[ii].lvl<<" )"<<std::endl;
//                    }
                }
            }
            glob_numb++;
        }
    }
    return result;
}

template<short_t d, class T>
gsSparseMatrix<T> gsHBSplineBasis<d,T>::coarsening_direct2( const std::vector<gsSortedVector<index_t> >& old,
                                                     const std::vector<gsSortedVector<index_t> >& n,
                                                     const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const
{
    int size1 = 0, size2 = 0;
    int glob_numb = 0;//continous numbering of hierarchical basis
    for(unsigned int i =0; i< old.size();i++){//count the number of basis functions in old basis
        size1 += old[i].size();
    }
    for(unsigned int i =0; i< n.size();i++){//count the number of basis functions in new basis
        size2 += n[i].size();
    }
    gsSparseMatrix<T> result(size2,size1);


//    std::vector<gsMatrix<T> > transferDense;// = transfer;
//    transferDense.resize(transfer.size());
//    for (unsigned int i = 0; i < transfer.size();i++){
//        transferDense[i] = transfer[i];
//    }
    std::vector<gsSparseMatrix<T,ColMajor> > temptransfer;// = transfer;
    temptransfer.resize(transfer.size());
    for (unsigned int i = 0; i < transfer.size();i++){
        temptransfer[i] = transfer[i];
        //gsDebug<<"transfer"<<i<<"\n"<<transfer[i]<<std::endl;
    }

    for (unsigned int i = 0; i < old.size(); i++)//iteration through the levels of the old basis
    {
        // find starting index of level i in new basis
        int start_lv_i = 0;
        for(unsigned int l =0; l < i; l++)
        {
            start_lv_i += n[l].size();
        }


        for (unsigned int j = 0; j < old[i].size();j++)//iteration through the basis functions in the given level
        {
            //gsDebug<<"j = "<<j<<std::endl;
            start_lv_i = 0;
            for(unsigned int l =0; l < i; l++)
            {
                start_lv_i += n[l].size();
            }
            const unsigned old_ij = old[i][j];  // tensor product index

            if( n[i].bContains(old_ij) )//it he basis function was not refined
            {
                result(start_lv_i + std::distance(n[i].begin(), n[i].find_it_or_fail(old_ij) ),glob_numb ) = 1;//settign the coefficient of the not refined basis function to 1
            }
            else
            {

                gsMatrix<index_t, d, 2> supp(d, 2);
                this->m_bases[i]->elementSupport_into(old_ij, supp);//this->support(start_lv_i+old_ij);
                //gsDebug<<"supp "<< supp<<std::endl;
                //unsigned max_lvl = math::min<index_t>( this->m_tree.query4(supp.col(0),supp.col(1), i), transfer.size() ) ;
                gsSparseVector<T,RowMajor> t(this->m_bases[i]->size());
                t.setZero();
                t[old_ij] = 1;
                //gsDebug<<"j = "<<j<<std::endl;
                //gsDebug<<"nsize"<<n.size()<<std::endl;
                for(unsigned int k = i+1; k < n.size();k++){
                    start_lv_i = 0;
                    for(unsigned int l =0; l < k-1; l++)
                    {
                        start_lv_i += n[l].size();
                    }
//                    gsDebug<<"k "<<k<<std::endl;
//                    gsDebug<<"nk"<<std::endl;
//                    for(int a = 0; a < n[k].size();a++){
//                        gsDebug<<n[k][a]<<" ";
//                    }
                    //gsDebug<<std::endl;
                    gsSparseVector<T,RowMajor> M;
                    M.setZero();
                    M = temptransfer[k-1] * t.transpose();

                    //gsDebug<<"M\n"<<M<<std::endl;
                    for(int l = 0 ; l < M.size();l++)
                    //for(typename gsSparseVector<T,RowMajor>::InnerIterator l(M); l; ++l)//basis function was refined->looking for the coeficinets from the global transfer matrix
                    {
                        if(M[l]!=0)
                            if(n[k].bContains(l))
                            {
                                //gsDebug<<"j: "<<j<<" "<<"oldij"<<old_ij<<" "<<"l:"<<l<<"    ";
                                //gsDebug<<"k "<<k<<" j: "<<j<<" "<<"globnumb "<<glob_numb<<" "<<"l:"<<l<<" "<<M[l]<<"  ";
                                const int pos = start_lv_i + n[k-1].size() + std::distance(n[k].begin(), n[k].find_it_or_fail(l));
                                //gsDebug<<"pos:"<<pos<<std::endl;
                                //double ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                                result(pos,glob_numb) = M[l];//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                                M[l] = 0;
                                //const int pos = start_lv_i + n[coeffs[0].lvl].size() + std::distance(n[coeffs[0].lvl+1].begin(), n[coeffs[0].lvl+1].find_it_or_fail(k.row()));
                                //double ppp =  transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                                //result(pos,glob_numb) += coeffs[0].coef * temptransfer[coeffs[0].lvl](k.row(), coeffs[0].pos);//transferDense[coeffs[0].lvl](k, coeffs[0].pos);
                            }
                    }
                    t = M;
                    //gsDebug<<"M after \n"<<M<<std::endl;
                }
            }
            glob_numb++;
        }

    }

    return result;
}
namespace internal
{

/// Get a Hierarchical B-spline basis from XML data
template<short_t d, class T>
class gsXml< gsHBSplineBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsHBSplineBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "HBSplineBasis"+ (d>1 ? to_string(d):""); }
    
    static gsHBSplineBasis<d,T> * get (gsXmlNode * node)
    {
        return getHTensorBasisFromXml< gsHBSplineBasis<d,T> > (node);
    }
  
    static gsXmlNode * put (const gsHBSplineBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putHTensorBasisToXml< gsHBSplineBasis<d,T> > (obj, data);
    }
};

}


}// namespace gismo
