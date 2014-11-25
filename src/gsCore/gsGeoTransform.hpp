/*
* gsGeoTransform.hpp created on 05.08.2014
*
* Author: Andrea Bressan
*
* This file is part of the G+SMO library
*/

#pragma once

#include <gsCore/gsBasisEvaluator.h>

// TODO flag system should be improved in order to avoid computing useless stuff

namespace gismo {

template <typename T, int ParDim, int TarDim, int AmbDim>
struct gsGeoNoTransform
{
    static unsigned addAuxiliaryFlags (unsigned flags)
    {
        if (TarDim!=1 && flags&NEED_HESSIAN)
        {
            GISMO_ASSERT(false, "can compute Hessian of scalar functions only.\n");
        }
        if (ParDim!=TarDim)
        {
            if (flags & (NEED_CURL | NEED_DIV) )
            {
                GISMO_ASSERT(false, "For curl and div the dimension of parametric domain and target domain must agree.\n");
            }
        }
        if (flags & (NEED_CURL | NEED_DIV | NEED_JACOBIAN) )
            flags |= NEED_GRAD;
        if (flags & (NEED_HESSIAN | NEED_LAPLACIAN) )
            flags |= NEED_2ND_DER;
        return flags;
    }
    static unsigned getGeometryFlags (unsigned flags)
    {
        return 0 & flags;
    }


    static void computeValues  (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>           (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        result.setZero(numA*TarDim,numP);

        int start=0;
        for(int i=0; i<TarDim;++i)
        {
            for( int j=0; j< activeNum[i] ; ++j)
            {
                result.row(start+TarDim*j+i)=b_values[i].row(j);
            }
            start+=TarDim*activeNum[i];
        }
    }
    static void computeJacobians   (        const gsBasisEvaluator<T>    *b_eval,
                                            const gsGeometryEvaluator<T> *g_Eval,
                                            const gsMatrix<T>             (&b_values)[TarDim],
                                            const int                     activeNum[TarDim],
                                                  gsMatrix<T>            &result
                                            )
    {
        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        result.setZero(numA*TarDim,ParDim*numP);

        for (int a=0; a<numA; ++a)
        {
            for (int c=0; c<TarDim; ++c)
            {
                for (int p=0; p<numP; ++p)
                {
                result.template block<1,ParDim>(a*TarDim+c,p*ParDim) = b_eval->derivs().template block<ParDim,1>(a*ParDim*TarDim+c*ParDim,p).transpose() ;
                }
            }
        }
    }

    static void computeGrads   (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        result.resize(numA*TarDim*ParDim,numP);

        int ao =0; // active origin
        int aoc=0; // active origin control
        for (int a=0; a<numA; ++a, ++aoc )
        {
            if (aoc==activeNum[ao])
            {
                ++ao;
                aoc=0;
            }
            for (int c=0; c<TarDim; ++c)
            {
                if (ao!=c)
                {
                    result.middleRows(a*ParDim*TarDim+c*ParDim,ParDim).setZero();
                }
                else
                {
                    result.middleRows(a*ParDim*TarDim+c*ParDim,ParDim)=
                        b_values[ao].middleRows(activeNum[ao]+aoc*ParDim,ParDim);
                }
            }
        }
    }

    static void computeDivs    (            const gsBasisEvaluator<T>    *b_eval,
                                            const gsGeometryEvaluator<T> *g_Eval,
                                            const gsMatrix<T>             (&b_values)[TarDim],
                                            const int                     activeNum[TarDim],
                                                  gsMatrix<T>            &result
                                            )
    {
        GISMO_ASSERT(ParDim==TarDim,"DIV exists only for tangent vectors");

        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();

        result.setZero(numA,numP);

        for( int j=0; j< numA ; ++j)
        {
            for(int i=0; i<TarDim;++i)
            {
            // TODO specialize for 2D and 3D so that we loop only once
            // and we can remove the setZero
            // 2D result.row(j)=b_eval->derivs.row(ParDim*(TarDim*j+i))+b_eval->derivs.row(ParDim*(TarDim*j+i)+1);
                result.row(j)+=b_eval->derivs().row(ParDim*(TarDim*j+i)+i);
            }
        }
    }
    static void computeCurls   (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        GISMO_ASSERT(ParDim==TarDim && TarDim==3,"CURL exists only for 3D tangent vectors");
        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        if (ParDim==3)
        {
            result.resize(TarDim*numA,numP);
            for (int a=0; a<numA; ++a)
                {
                    result.row( TarDim*a) =
                        b_eval->derivs().row(TarDim*ParDim*a+5)- b_eval->derivs().row(TarDim*ParDim*a+7);
                    result.row( TarDim*a+1) =
                        b_eval->derivs().row(TarDim*ParDim*a+6)- b_eval->derivs().row(TarDim*ParDim*a+2);
                    result.row( TarDim*a+2) =
                        b_eval->derivs().row(TarDim*ParDim*a+1)- b_eval->derivs().row(TarDim*ParDim*a+3);
                }
        }

    }




    /**
     * @brief computeSecDers
     * @param b_eval
     * @param g_Eval
     * @param activeNum
     * @param result resul( i*TarDim*TarDim + j*TarDim + k, pt ) =
     * second derivative of (possibly vector-valued) function with index \em i,
     * with respect to the <em>j</em>-th and <em>k</em>-th variable,
     * evaluated at point with index \em pt ...most likely... have to check again...
     */
    // Comments by S. Kleiss, to be verified by Andrea whether I
    // understood correctly
    static void computeSecDers (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        // number of evaluation points:
        const int numP = b_values[0].cols();
        // total number of active (possibly vector-valued) functions
        const int numA = b_eval->actives().rows();

        // for each component, there are ParDim*ParDim mixed second derivatives:
        const int secDerSize1D = ParDim*ParDim;
        // offset for all the function values and gradients that are also stored
        // in b_values
        const int secDerOff1D  = 1+ParDim;
        // for each active function, the number of second derivatives is secDerSize:
        // number of mixed derivatives per component times number of components.
        const int secDerSize   = secDerSize1D*TarDim;

        result.resize(numA*secDerSize,numP);

        // active origin, loops over the components of the (vector-valued) function
        int ao =0;
        // component active index
        int ac_id=0;
        int at_id=0; // total active index
        for (; at_id<numA; ++at_id, ++ac_id )
        {
            if (ac_id==activeNum[ao])
            {
                // when all active functions of the component have been touched...
                ++ao;     // ...proceed to the next component, and...
                ac_id=0;  // ...reset the active index of that component to zero
            }
            for (int c=0; c<TarDim; ++c)
            {
                if (ao!=c)
                {
                    // if the current targent dimension does not equal the
                    // active component, fill up with zeros.
                    result.middleRows(at_id*secDerSize+c*secDerSize1D,secDerSize1D).setZero();
                }
                else
                {
                    // otherwise, compute the second derivatives
                    result.middleRows(at_id*secDerSize+c*secDerSize1D,secDerSize1D)=
                        b_values[ao].middleRows(activeNum[ao]*secDerOff1D+ac_id*secDerSize1D,secDerSize1D);
                }
            }
        }
    }

    static void computeLaplacians    (  const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();

        const int secDerSize1D = ParDim*ParDim;
        const int secDerSize   = secDerSize1D*TarDim;

        result.resize(numA*TarDim,numP);

        for (int a=0; a<numA; ++a)
        {
            for (int c=0; c<TarDim; ++c)
            {
                result.row(a*TarDim+c) = b_eval->derivs2().middleRows(a*secDerSize,TarDim).colwise().sum();
            }
        }

    }
};


template <typename T, int ParDim, int TarDim, int AmbDim>
struct gsGeoGradPreservingTransform
        : public gsGeoNoTransform<T,ParDim,TarDim, AmbDim>
{
    static unsigned addAuxiliaryFlags (unsigned flags)
    {
        if (TarDim!=1 && flags&NEED_HESSIAN)
        {
            GISMO_ASSERT(false, "can compute Hessian of scalar functions only.\n");
        }
        if (ParDim!=TarDim)
        {
            if (flags & (NEED_CURL | NEED_DIV) )
            {
                GISMO_ASSERT(false, "For curl and div the dimension of parametric domain and target domain must agree.\n");
            }
        }
        if (flags & (NEED_CURL | NEED_DIV | NEED_JACOBIAN) )
            flags |= NEED_GRAD;
        if (flags & (NEED_HESSIAN | NEED_LAPLACIAN) )
            flags |= NEED_2ND_DER;
        return flags;
    }
    static unsigned getGeometryFlags (unsigned flags)
    {
        unsigned geoFlags=0;
        if (flags & NEED_GRAD)
            geoFlags |= NEED_GRAD_TRANSFORM;
        if (flags & NEED_2ND_DER )
            geoFlags |= NEED_2ND_DER;
        return geoFlags;
    }

    /**
     * @brief computeSecDers
     *
     *
     *
     * @param b_eval
     * @param g_Eval
     * @param activeNum
     * @param result
     */
    static void computeSecDers (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        // TODO implement transformation for 2nd derivatives we have the formula
        GISMO_ASSERT(false, "waiting for Stefan to do it :)");

        //qwe

/*
 *
*/
    }

    static void computeGrads   (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        result.resize(numA*TarDim*AmbDim,numP);

        for (int p=0; p<numP; ++p)
        {
            int ag_id = 0;
            for (int b_id=0; b_id<TarDim; ++b_id)
            {
                for ( int a_id=0; a_id<activeNum[b_id]; ++a_id)
                {
                    for (int c=0; c<TarDim; ++c)
                    {
                        if (b_id !=c)
                        {
                            result.template block<AmbDim,1>((a_id + ag_id)*AmbDim*TarDim+c*AmbDim,p).setZero();
                        }
                        else
                        {
                            result.template block<AmbDim,1>((a_id + ag_id)*AmbDim*TarDim+c*AmbDim,p) =
                                g_Eval->gradTransform(p)*b_values[b_id].template block<ParDim,1>(activeNum[b_id]+a_id*ParDim,p);
                        }
                    }
                }
                ag_id += activeNum[b_id];
            }
        }
    }

};


template <typename T, int ParDim, int TarDim, int AmbDim>
struct gsGeoDivPreservingTransform
        : public gsGeoNoTransform<T,ParDim,TarDim, AmbDim>
{

    static unsigned addAuxiliaryFlags (unsigned flags)
    {
        if (ParDim != TarDim)
        {
            GISMO_ASSERT(false, "Can not preserve divergence if the target and parameter dimension disagree.\n");
        }
        if ( flags & NEED_HESSIAN)
        {
            GISMO_ASSERT(false, "can compute Hessian of scalar functions only.\n");
        }
        if (flags & (NEED_CURL | NEED_DIV | NEED_JACOBIAN) )
            flags |= NEED_GRAD;
        if (flags & (NEED_HESSIAN | NEED_LAPLACIAN) )
            flags |= NEED_2ND_DER;
        return flags;
    }
    static unsigned getGeometryFlags (unsigned flags)
    {
        unsigned geoFlags=0;
        if (flags & NEED_GRAD)
            geoFlags |= NEED_GRAD_TRANSFORM | NEED_2ND_DER;
        if (flags & NEED_2ND_DER )
            geoFlags |=NEED_2ND_DER;
        return geoFlags;
    }

    static void computeValues  (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {

        //TODO make assert if g_Eval is null

        // assemble value matrix one point per column
        // each group of TarDim rows contains the transformed values
        // of the basis vector


        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        result.setZero(numA*TarDim,numP);


        const gsVector<T> & det = g_Eval->measures();
        const gsMatrix<T> & jacs = g_Eval->jacobians(); //TarDim X TarDim*numP
        //NB The jacobian might need to be transpoes

        int start=0;
        for(int comp = 0; comp < TarDim; ++comp)
        {
            for( int j=0; j< activeNum[comp] ; ++j)
            {
                for( int p=0; p< numP ; ++p)
                {
                    result.template block<TarDim,1>(start+TarDim*j, p) +=
                            jacs.template block<TarDim,1>(0,TarDim*p+comp)*b_values[comp](j,p)/det(p);
                }
            }
            start+=TarDim*activeNum[comp];
        }
    }

    static void computeGrads   (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        //Assumtions:
        //GeoDim = ParDim = TargetDim

        const int numP=b_values[0].cols();
        const int numA=b_eval->actives().rows();
        const gsMatrix<T> & jacs    = g_Eval->jacobians(); //TarDim X TarDim*numP
        const gsMatrix<T> & invJacs = g_Eval->gradTransforms().transpose();
        const gsMatrix<T> & secDers = g_Eval->derivs2();
        const gsVector<T> & det     = g_Eval->measures();
        T detSigned;
        //gsMatrix<T> DJacx(ParDim,ParDim);
        result.setZero(numA*TarDim*ParDim,numP);

        // Number of second derivatives
        const index_t k2 = (ParDim + (ParDim*(ParDim-1))/2);
        index_t k1 = 0;                         //Index for Jacobian
        int ao =0; // active origin
        int aoc=0; // active origin control



        for (index_t k=0; k < numP; ++k)
        {
            k1 = k*ParDim;

            detSigned = g_Eval->orientation() * det[k];
            //equivalent     :g_Eval->jacobian(k).determinant();
            //equivalent     : g_Eval->jacobians().middleCols(k1,ParDim);
            //equivalent (2d): jacs(0,0+k1)*jacs(1,1+k1)-jacs(0,1+k1)*jacs(1,0+k1);

            //JS2: Verified for 2D case (with NumP = 1)!
            if(TarDim==2)
            {
                //JS2: CODE GOTEN FROM GEOPDES SOFTWARE

                //xu = jacs(0,0+k1)  //xv = jacs(0,1+k1)
                //yu = jacs(1,0+k1)  //yv = jacs(1,1+k1)

                //xuu = secDers(0+0*k2,k) //yuu = secDers(0+1*k2,k)
                //xvv = secDers(1+0*k2,k) //yvv = secDers(1+1*k2,k)
                //xuv = secDers(2+0*k2,k) //yuv = secDers(2+1*k2,k)

                //det = det(k)

                //wh = u1h = b_values[0](aoc,k)
                //zh = u2h = b_values[1](aoc,k)
                //whu = b_values[0](activeNum[0]+aoc*ParDim+0,k)
                //whv = b_values[0](activeNum[0]+aoc*ParDim+1,k)
                //zhu = b_values[1](activeNum[1]+aoc*ParDim+0,k)
                //zhv = b_values[1](activeNum[1]+aoc*ParDim+1,k)

                int tmpInd = 0;

                // under construction
                //DJacx(0,0) = secDers(0+0*k2,k);
                //DJacx(1,0) = secDers(0+1*k2,k);
                //DJacx(0,1) = secDers(2+0*k2,k);
                //DJacx(1,1) = secDers(2+1*k2,k);
                //T detdx = detSigned * ( invJacs.middleCols(k1,TarDim) * DJacx * invJacs.middleCols(k1,TarDim)).trace() ;

// /*
                T detdx = jacs(1,0+k1)*(jacs(1,0+k1)*secDers(1+0*k2,k) + jacs(0,1+k1)*secDers(2+1*k2,k) - jacs(0,0+k1)*secDers(1+1*k2,k));
                detdx  += jacs(1,1+k1)*(jacs(1,1+k1)*secDers(0+0*k2,k) + jacs(0,0+k1)*secDers(2+1*k2,k) - jacs(0,1+k1)*secDers(0+1*k2,k));
                detdx  +=-2*jacs(1,0+k1)*jacs(1,1+k1)*secDers(2+0*k2,k);
                detdx  /= detSigned;
                //detdx = (yu.*(yu.*xvv + xv.*yuv - xu.*yvv) + ...
                //         yv.*(yv.*xuu + xu.*yuv - xv.*yuu) - 2*yu.*yv.*xuv)./det;
//*/

                T detdy = jacs(0,0+k1)*(jacs(1,1+k1)*secDers(2+0*k2,k) + jacs(0,0+k1)*secDers(1+1*k2,k) - jacs(1,0+k1)*secDers(1+0*k2,k));
                detdy  += jacs(0,1+k1)*(secDers(2+0*k2,k)*jacs(1,0+k1) + secDers(0+1*k2,k)*jacs(0,1+k1) - jacs(1,1+k1)*secDers(0+0*k2,k));
                detdy  +=-2*jacs(0,0+k1)*jacs(0,1+k1)*secDers(2+1*k2,k);
                detdy  /= detSigned;
                //detdy = (xu.*(yv.*xuv + xu.*yvv - yu.*xvv) + ...
                //         xv.*(xuv.*yu + yuu.*xv - yv.*xuu) - 2*xu.*xv.*yuv)./det;

                for (int a=0; a<numA; ++a, ++aoc )
                {
                    if (aoc==activeNum[ao])
                    {
                        ++ao;
                        aoc=0;
                        tmpInd = 2;
                    }
                    
                    //result(a*ParDim*TarDim+0*ParDim+0,k) = ...

                    //v1x
                    result(a*ParDim*TarDim+0*ParDim+0,k) += (-detdx*jacs(0,ao+k1)*b_values[ao](aoc,k) \
                            + jacs(1,1+k1)*(secDers(tmpInd+0*k2,k)*b_values[ao](aoc,k) + jacs(0,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+0,k) ) \
                            - jacs(1,0+k1)*(secDers(2-ao  +0*k2,k)*b_values[ao](aoc,k) + jacs(0,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+1,k) ))\
                            / (detSigned*det(k));

                    //v1y
                    result(a*ParDim*TarDim+0*ParDim+1,k) += (-detdy*jacs(0,ao+k1)*b_values[ao](aoc,k) \
                            + jacs(0,0+k1)*(secDers(2-ao  +0*k2,k)*b_values[ao](aoc,k) + jacs(0,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+1,k) ) \
                            - jacs(0,1+k1)*(secDers(tmpInd+0*k2,k)*b_values[ao](aoc,k) + jacs(0,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+0,k) ))\
                            / (detSigned*det(k));

                    //v2x
                    result(a*ParDim*TarDim+1*ParDim+0,k) += (-detdx*jacs(1,ao+k1)*b_values[ao](aoc,k) \
                            + jacs(1,1+k1)*(secDers(tmpInd+1*k2,k)*b_values[ao](aoc,k) + jacs(1,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+0,k) ) \
                            - jacs(1,0+k1)*(secDers(2-ao  +1*k2,k)*b_values[ao](aoc,k) + jacs(1,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+1,k) ))\
                            / (detSigned*det(k));

                    //v2y
                    result(a*ParDim*TarDim+1*ParDim+1,k) += (-detdy*jacs(1,ao+k1)*b_values[ao](aoc,k) \
                            + jacs(0,0+k1)*(secDers(2-ao  +1*k2,k)*b_values[ao](aoc,k) + jacs(1,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+1,k) ) \
                            - jacs(0,1+k1)*(secDers(tmpInd+1*k2,k)*b_values[ao](aoc,k) + jacs(1,ao+k1)*b_values[ao](activeNum[ao]+aoc*ParDim+0,k) ))\
                            / (detSigned*det(k));

                    /*
                    v1x = (-(xu.*wh + xv.*zh).*detdx +...
                                yv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu)...
                               -yu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv))./det2;

                    v1y = (-(xu.*wh + xv.*zh).*detdy +...
                                xu.*(xuv.*wh + xvv.*zh + xu.*whv + xv.*zhv)...
                               -xv.*(xuu.*wh + xuv.*zh + xu.*whu + xv.*zhu))./det2;

                    v2x = (-(yu.*wh + yv.*zh).*detdx +...
                                yv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu)...
                               -yu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv))./det2;

                    v2y = (-(yu.*wh + yv.*zh).*detdy +...
                                xu.*(yuv.*wh + yvv.*zh + yu.*whv + yv.*zhv)...
                               -xv.*(yuu.*wh + yuv.*zh + yu.*whu + yv.*zhu))./det2;
                    */
                }
                ao = 0;
                aoc=0;
            }

            else if(TarDim==3) //3D case
            {
                gsMatrix<T> xd(TarDim,TarDim), yd(TarDim,TarDim), zd(TarDim,TarDim);
                gsVector<T> determinantDerivative;
                determinantDerivative.setZero(TarDim);

                //might need to devide by the signed determinant
                /*detSigned = jacs(0,0+k1)*(jacs(1,1+k1)*jacs(2,2+k1)-jacs(1,2+k1)*jacs(2,1+k1))
                        -   jacs(0,1+k1)*(jacs(1,0+k1)*jacs(2,2+k1)-jacs(1,2+k1)*jacs(2,0+k1))
                        +   jacs(0,2+k1)*(jacs(1,0+k1)*jacs(2,1+k1)-jacs(1,1+k1)*jacs(2,0+k1));
                gsDebug <<"detSigned: "  << detSigned <<"\n";*/

                for (index_t i = 0; i < TarDim; ++i)
                {
                    xd(0,i) = secDers(0+0*k2,k)*invJacs(0,i+k1) + secDers(3+0*k2,k)*invJacs(1,i+k1) + secDers(4+0*k2,k)*invJacs(2,i+k1);
                    xd(1,i) = secDers(3+0*k2,k)*invJacs(0,i+k1) + secDers(1+0*k2,k)*invJacs(1,i+k1) + secDers(5+0*k2,k)*invJacs(2,i+k1);
                    xd(2,i) = secDers(4+0*k2,k)*invJacs(0,i+k1) + secDers(5+0*k2,k)*invJacs(1,i+k1) + secDers(2+0*k2,k)*invJacs(2,i+k1);
                    yd(0,i) = secDers(0+1*k2,k)*invJacs(0,i+k1) + secDers(3+1*k2,k)*invJacs(1,i+k1) + secDers(4+1*k2,k)*invJacs(2,i+k1);
                    yd(1,i) = secDers(3+1*k2,k)*invJacs(0,i+k1) + secDers(1+1*k2,k)*invJacs(1,i+k1) + secDers(5+1*k2,k)*invJacs(2,i+k1);
                    yd(2,i) = secDers(4+1*k2,k)*invJacs(0,i+k1) + secDers(5+1*k2,k)*invJacs(1,i+k1) + secDers(2+1*k2,k)*invJacs(2,i+k1);
                    zd(0,i) = secDers(0+2*k2,k)*invJacs(0,i+k1) + secDers(3+2*k2,k)*invJacs(1,i+k1) + secDers(4+2*k2,k)*invJacs(2,i+k1);
                    zd(1,i) = secDers(3+2*k2,k)*invJacs(0,i+k1) + secDers(1+2*k2,k)*invJacs(1,i+k1) + secDers(5+2*k2,k)*invJacs(2,i+k1);
                    zd(2,i) = secDers(4+2*k2,k)*invJacs(0,i+k1) + secDers(5+2*k2,k)*invJacs(1,i+k1) + secDers(2+2*k2,k)*invJacs(2,i+k1);

                    determinantDerivative[i] += xd(0,i)*invJacs(0,0+k1)+xd(1,i)*invJacs(1,0+k1)+xd(2,i)*invJacs(2,0+k1);
                    determinantDerivative[i] += yd(0,i)*invJacs(0,1+k1)+yd(1,i)*invJacs(1,1+k1)+yd(2,i)*invJacs(2,1+k1);
                    determinantDerivative[i] += zd(0,i)*invJacs(0,2+k1)+zd(1,i)*invJacs(1,2+k1)+zd(2,i)*invJacs(2,2+k1);

                }

                // under construction
                //DJacx(0,0) = secDers(0+0*k2,k);
                //DJacx(1,0) = secDers(0+1*k2,k);
                //DJacx(2,0) = secDers(0+2*k2,k);
                //DJacx(0,1) = secDers(3+0*k2,k);
                //DJacx(1,1) = secDers(3+1*k2,k);
                //DJacx(2,1) = secDers(3+2*k2,k);
                //DJacx(0,2) = secDers(5+0*k2,k);
                //DJacx(1,2) = secDers(5+1*k2,k);
                //DJacx(2,2) = secDers(5+2*k2,k);
                //determinantDerivative[0] = 
                    //detSigned * 
                //( invJacs.middleCols(k1,TarDim) * DJacx * invJacs.middleCols(k1,TarDim)).trace() ;

                for (int a=0; a<numA; ++a, ++aoc )
                {
                    if (aoc==activeNum[ao])
                    {
                        ++ao;
                        aoc=0;
                    }

                    T tmpGradJ1,tmpGradJ2,tmpGradJ3;

                    tmpGradJ1 = b_values[ao](activeNum[ao]+aoc*ParDim+0,k)*invJacs(0,0) \
                            + b_values[ao](activeNum[ao]+aoc*ParDim+1,k)*invJacs(1,0) \
                            + b_values[ao](activeNum[ao]+aoc*ParDim+2,k)*invJacs(2,0);

                    tmpGradJ2 = b_values[ao](activeNum[ao]+aoc*ParDim+0,k)*invJacs(0,1) \
                            + b_values[ao](activeNum[ao]+aoc*ParDim+1,k)*invJacs(1,1) \
                            + b_values[ao](activeNum[ao]+aoc*ParDim+2,k)*invJacs(2,1);

                    tmpGradJ3 = b_values[ao](activeNum[ao]+aoc*ParDim+0,k)*invJacs(0,2) \
                            + b_values[ao](activeNum[ao]+aoc*ParDim+1,k)*invJacs(1,2) \
                            + b_values[ao](activeNum[ao]+aoc*ParDim+2,k)*invJacs(2,2);

                    result(a*ParDim*TarDim+0*ParDim+0,k) += \
                            (-determinantDerivative[0]*jacs(0,ao+k1)*b_values[ao](aoc,k) \
                            + xd(ao,0)*b_values[ao](aoc,k) + jacs(0,ao+k1)*tmpGradJ1)/det(k);

                    result(a*ParDim*TarDim+0*ParDim+1,k) += \
                            (-determinantDerivative[1]*jacs(0,ao+k1)*b_values[ao](aoc,k) \
                            + xd(ao,1)*b_values[ao](aoc,k) + jacs(0,ao+k1)*tmpGradJ2)/det(k);

                    result(a*ParDim*TarDim+0*ParDim+2,k) += \
                            (-determinantDerivative[2]*jacs(0,ao+k1)*b_values[ao](aoc,k) \
                            + xd(ao,2)*b_values[ao](aoc,k) + jacs(0,ao+k1)*tmpGradJ3)/det(k);


                    result(a*ParDim*TarDim+1*ParDim+0,k) +=  \
                            (-determinantDerivative[0]*jacs(1,ao+k1)*b_values[ao](aoc,k) \
                            + yd(ao,0)*b_values[ao](aoc,k) + jacs(1,ao+k1)*tmpGradJ1)/det(k);

                    result(a*ParDim*TarDim+1*ParDim+1,k) += \
                            (-determinantDerivative[1]*jacs(1,ao+k1)*b_values[ao](aoc,k) \
                            + yd(ao,1)*b_values[ao](aoc,k) + jacs(1,ao+k1)*tmpGradJ2)/det(k);

                    result(a*ParDim*TarDim+1*ParDim+2,k) += \
                            (-determinantDerivative[2]*jacs(1,ao+k1)*b_values[ao](aoc,k) \
                            + yd(ao,2)*b_values[ao](aoc,k) + jacs(1,ao+k1)*tmpGradJ3)/det(k);


                    result(a*ParDim*TarDim+2*ParDim+0,k) += \
                            (-determinantDerivative[0]*jacs(2,ao+k1)*b_values[ao](aoc,k) \
                            + zd(ao,0)*b_values[ao](aoc,k) + jacs(2,ao+k1)*tmpGradJ1)/det(k);

                    result(a*ParDim*TarDim+2*ParDim+1,k) += \
                            (-determinantDerivative[1]*jacs(2,ao+k1)*b_values[ao](aoc,k) \
                            + zd(ao,1)*b_values[ao](aoc,k) + jacs(2,ao+k1)*tmpGradJ2)/det(k);

                    result(a*ParDim*TarDim+2*ParDim+2,k) += \
                            (-determinantDerivative[2]*jacs(2,ao+k1)*b_values[ao](aoc,k) \
                            + zd(ao,2)*b_values[ao](aoc,k) + jacs(2,ao+k1)*tmpGradJ3)/det(k);
                }
            }
            else
            {
                GISMO_ERROR("Geometric dimention not valid");
            }
        }
        ao  = 0;
        aoc = 0;
    }

    /*
    static void computeDivs   (         const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        // divergence is computed on the parametric domain and then mapped should be
        // as easy as multiplying by the determinant of the Jacobiam i.e. g_Eval->measures()
        // maybe sign issues with for orientation inverting maps
    }*/


    static void computeSecDers (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        GISMO_ASSERT(false, "NOT IMPLEMENTED. DO IT YOURSELF. PLEASE.");
    }

};



template <typename T, int ParDim, int TarDim, int AmbDim>
struct gsGeoCurlPreservingTransform
        : public gsGeoNoTransform<T,ParDim,TarDim, AmbDim>
{

    static void computeSecDers (        const gsBasisEvaluator<T>    *b_eval,
                                        const gsGeometryEvaluator<T> *g_Eval,
                                        const gsMatrix<T>             (&b_values)[TarDim],
                                        const int                     activeNum[TarDim],
                                              gsMatrix<T>            &result
                                        )
    {
        GISMO_ASSERT(false, "NOT IMPLEMENTED. DO IT YOURSELF. PLEASE.");
    }
};


} // namespace gismo
