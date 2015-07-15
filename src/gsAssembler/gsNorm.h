/** @file gsNorm.h

    @brief Provides generic routines for computing function and error norms

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once


namespace gismo
{

/** @brief The gsNorm class provides to generic routines for
 * computing function norms and distances as well as error estimates
 * element-wise or globally
 *
 * \ingroup Assembler
*/
template <class T>
class gsNorm
{
public:
//typedef typename NormVisitor::Scalar_t T;

    /// Constructor using a multipatch domain
    gsNorm(const gsField<T> & _field1,
           const gsFunction<T> & _func2) 
    : m_zeroFunction(gsConstantFunction<T>(T(0.0),_field1.parDim())), patchesPtr( &_field1.patches() ),
      field1(&_field1), func2(&_func2)
    { }

    /// Constructor using a multipatch domain
    gsNorm(const gsField<T> & _field1)
    : m_zeroFunction(gsConstantFunction<T>(T(0.0),_field1.parDim())), patchesPtr( &_field1.patches() ),
      field1(&_field1), func2(&m_zeroFunction)
    { }


    void setField(const gsField<T> & _field1)
    {
        field1 = &_field1;
    }

    /** \brief Main function for norm-computation.
     *
     * The computed value can be accessed by value().
     *
     * \param[in] visitor The Norm-visitor to be used.
     * \param[in] storeElWise Flag indicating whether the
     * element-wise norms should be stored. See also
     * elementNorms().
     * \param[in] side To be used, if the visitor will
     * iterate over a side of the domain.
     */
    template <class NormVisitor>
    void apply(NormVisitor & visitor, bool storeElWise = false, boxSide side = boundary::none)
    {
        if ( storeElWise )
            m_elWise.clear();

        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights
        gsQuadRule<T> QuRule; // Reference Quadrature rule

        // Evaluation flags for the Geometry map
        unsigned evFlags(0);

        m_value = T(0.0);
        for (unsigned pn=0; pn < patchesPtr->nPatches(); ++pn )// for all patches
        {
            const gsGeometry<T> & func1 = field1->igaFunction(pn);

            // Initialize visitor
            visitor.initialize(func1.basis(), QuRule, evFlags);

            // Initialize geometry evaluator
            typename gsGeometry<T>::Evaluator geoEval(
                patchesPtr->patch(pn).evaluator(evFlags));
            
            typename gsBasis<T>::domainIter domIt = func1.basis().makeDomainIterator(side);
            for (; domIt->good(); domIt->next())
            {
                // Map the Quadrature rule to the element
                QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Evaluate on quadrature points
                visitor.evaluate(*geoEval, func1, *func2, quNodes);
                
                // Accumulate value from the current element (squared)
                const T result = visitor.compute(*domIt, *geoEval, quWeights);
                m_value += result;
                if ( storeElWise )
                    m_elWise.push_back( math::sqrt(result) );
            }
        }

        m_value = math::sqrt(m_value);
    }
    
public:

    /// Return the multipatch
    const gsMultiPatch<T> & patches() const { return *patchesPtr; }

    /** @brief Returns the computed norm values element-wise
     *
     *
     * \returns The SQUARED element-norms in a gsVector.\n
     * It is assumed that they were actually computed by providing
     * the proper flag \em storeEltWise in the call of apply().
     *
     * \remarks The order of the element-norms is "defined"
     * firstly by the numbering of the patches, and secondly
     * by the order in which the gsDomainIterator iterates
     * over the elements of the mesh! As of now (02.Dec.2014),
     * there is no way to access a particularly numbered element
     * on a mesh directly; you have to run the corresponding
     * gsDomainIterator again to reach the respective element.
     *
     */
    const std::vector<T> & elementNorms() const { return m_elWise; }

    /// @brief Returns the computed norm value
    T value() const { return m_value; }

private:
    // constant function zero that is used  when no function is specified
    // to calculate the norm and not the distance
    // maybe its possible to optimize this
    const gsConstantFunction<T> m_zeroFunction;

protected:

    const gsMultiPatch<T> * patchesPtr;

    const gsField<T>    * field1;

    const gsFunction<T> * func2;
  
protected:

    std::vector<T> m_elWise;
    T              m_value;
};

} // namespace gismo

