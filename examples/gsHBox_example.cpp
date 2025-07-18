/** @file gsHBox_example

    @brief Demonstrates functionality of the gsHBox

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (TU Delft 2019-...)
*/

#include <iostream>

#include <gismo.h>
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementHelper.h>

// #include <gsUtils/gsCombinatorics.h>

using namespace gismo;

// template <short_t d, class T>
// class gsHElementMarker
// {

// public:
//     struct CompareElementErrorPair;

// public:
//     typedef          gsHElement<d,T>                                    element_t;
//     typedef typename gsHElement<d,T>::level_t                           level_t;
//     typedef typename gsHElement<d,T>::Compare                           CompareElement;
//     typedef typename std::set<element_t,typename element_t::Compare>    HElementContainer;
//     typedef          T                                                  error_t;
//     typedef typename std::vector<std::pair<element_t, error_t>>         HElementErrorContainer;
//     // typedef typename std::map<element_t, index_t, CompareElement>       HElementMap;
//     // typedef typename std::vector<

// public:
//     struct CompareElementErrorPair
//     {
//         bool operator()(const std::pair<element_t, error_t> & a, const std::pair<element_t, error_t> & b) const
//         {
//             CompareElement comp;
//             return a.second < b.second || (a.second == b.second && comp(a.first,b.first));
//         };
//     };

// public:

//     gsHElementMarker(const gsHTensorBasis<d,T> & basis, gsOptionList options = defaultOptions())
//     :
//     m_basis(basis),
//     m_helper(basis)
//     {
//     }

//     gsOptionList & options() { return m_options; }

//     static gsOptionList defaultOptions()
//     {
//         gsOptionList options;
//         options.addInt("CoarsenRule","Rule used for coarsening: 1=GARU, 2=PUCA, 3=BULK.",1);
//         options.addInt("RefineRule","Rule used for refinement: 1=GARU, 2=PUCA, 3=BULK.",1);
//         options.addReal("CoarsenParam","Parameter used for coarsening",0.1);
//         options.addReal("RefineParam","Parameter used for refinement",0.1);
//         options.addInt("MaxLevel","Maximum refinement level",3);
//         options.addInt("Admissibility","Admissibility region, 0=T-admissibility (default), 1=H-admissibility",0);
//         options.addSwitch("Admissible","Mark the admissible region",true);
//         options.addInt("Jump","Jump parameter m",2);
//         // options.addInt("Verbose","Verbosity level",0);
//         return options;
//     }

//     void setErrors(const std::vector<error_t> & errors)
//     {
//         GISMO_ASSERT(errors.size() == m_basis.numElements(),
//                      "The size of the errors vector must match the number of elements in the basis domain.");

//         m_elementErrors.resize(m_basis.numElements());
//         element_t element;
//         level_t elementLevel;
//         for (auto it = m_basis.domain()->beginAll(); it!= m_basis.domain()->endAll(); ++it)
//         {
//             // Create a pair of element and its associated error
//             elementLevel = static_cast<const gsHDomainIterator<T,d> *>(&domIt)->getLevel();
//             element = m_helper.toElement(it->lowerCorner(), it->upperCorner(), elementLevel, it->patch());
//             // m_elementErrors[elem] = it.id();
//             m_elementErrors[it.id()] = std::make_pair(element, errors[it.id()]);
//         }
//         // Sort the element errors based on the error value in non-decreasing order
//         std::stable_sort(m_elementErrors.begin(), m_elementErrors.end(),CompareElementErrorPair());
//     }

//     HElementContainer markRef() const
//     {

//     }

// private:

//     HElementContainer _markRef_threshold() const
//     {
//         HElementContainer result;
//         T threshold = m_options.askReal("RefineParam",0.1) * m_elementErrors.back().second;
//         for (std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it)
//         {
//             // If the error is below the threshold, stop the iteration
//             if (it->second < threshold)
//                 break;

//             // If the level of the element is larger than the maximum level, skip it
//             if (it->first.level() >= m_options.askInt("MaxLevel",10))
//                 continue;

//             // Add the element to the container
//             result.insert(it->first);
//         }

//         if (m_options.askSwitch("Admissible",true))
//         {
//             // Mark the admissible region
//             result = m_helper.markAdmissible(result, m_options.askInt("Jump",2));
//         }
//     }

//     HElementContainer _markCrs_threshold(const HElementContainer & refined = {}) const
//     {
//         HElementContainer result;
//         T threshold = m_options.askReal("CoarsenParam",0.1) * m_elementErrors.back().second;
//         for (std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it)
//         {
//             // If the error is above the threshold, stop the iteration
//             if (it->second > threshold)
//                 break;

//             // If the element or its siblings are in the refined set, skip it
//             if (refined.empty())
//             {
//                 HElementContainer siblings = m_helper.getSiblings(it->first,false);
//                 if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
//                     continue;
//             }

//             // Add the element to the container
//             result.insert(it->first);
//         }

//         if (m_options.askSwitch("Admissible",true))
//         {
//             GISMO_ERROR("TODO");
//             // // Mark the admissible region
//             // result = m_helper.markAdmissible(result, m_options.askInt("Jump",2));
//         }
//         return result;
//     }

//     HElementContainer _markRef_percentage() const
//     {
//         HElementContainer result;
//         T percentage = m_options.askReal("RefineParam",0.1);
//         index_t numElements = m_elementErrors.size();
//         index_t numToMark = static_cast<index_t>(math::floor(percentage * numElements));
//         index_t numMarked = 0;
//         for (std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it, ++numMarked)
//         {
//             // If we have marked enough elements, stop the iteration
//             if (numMarked >= numToMark)
//                 break;

//             // If the level of the element is larger than the maximum level, skip it
//             if (it->first.level() >= m_options.askInt("MaxLevel",10))
//                 continue;

//             // Add the element to the container
//             result.insert(it->first);
//         }

//         if (m_options.askSwitch("Admissible",true))
//         {
//             // Mark the admissible region
//             result = m_helper.markAdmissible(result, m_options.askInt("Jump",2));
//         }
//         return result;
//     }

//     HElementContainer _markCrs_percentage(const HElementContainer & refined = {}) const
//     {
//         HElementContainer result;
//         T percentage = m_options.askReal("CoarsenParam",0.1);
//         index_t numElements = m_elementErrors.size();
//         index_t numToMark = static_cast<index_t>(math::floor(percentage * numElements));
//         index_t numMarked = 0;
//         for (std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it, ++numMarked)
//         {
//             // If we have marked enough elements, stop the iteration
//             if (numMarked >= numToMark)
//                 break;

//             // If the element or its siblings are in the refined set, skip it
//             if (!refined.empty())
//             {
//                 HElementContainer siblings = m_helper.getSiblings(it->first,false);
//                 if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
//                     continue;
//             }

//             // Add the element to the container
//             result.insert(it->first);
//         }

//         if (m_options.askSwitch("Admissible",true))
//         {
//             GISMO_ERROR("TODO");
//             // // Mark the admissible region
//             // result = m_helper.markAdmissible(result, m_options.askInt("Jump",2));
//         }
//     }

//     // NOTE: This does not take the extra error contributions due to admissibility into account.
//     HElementContainer markRef_fraction() const
//     {
//         // Compute the total error
//         T cummulErrMarked = T(0);
//         T totalError = std::accumulate(m_elementErrors.begin(), m_elementErrors.end(), T(0),[](T sum, const std::pair<element_t, error_t> & elem) { return sum + elem.second; });
//         T errorMarkSum = m_options.askReal("RefineParam",0.1) * totalError;
//         for (std::vector<std::pair<element_t, error_t>>::const_reverse_iterator it = m_elementErrors.rbegin(); it != m_elementErrors.rend(); ++it, ++numMarked)
//         {
//             // If the cumulative error exceeds the threshold, stop the iteration
//             if (cummulErrMarked >= errorMarkSum)
//                 break;

//             // If the level of the element is larger than the maximum level, skip it
//             if (it->first.level() >= m_options.askInt("MaxLevel",10))
//                 continue;

//             // Add the element to the container
//             result.insert(it->first);
//             cummulErrMarked += it->second;
//         }

//         if (m_options.askSwitch("Admissible",true))
//         {
//             // Mark the admissible region
//             result = m_helper.markAdmissible(result, m_options.askInt("Jump",2));
//         }
//         return result;
//     }

//     HElementContainer markCrs_fraction(const HElementContainer & refined = {}) const
//     {
//         // Compute the total error
//         T cummulErrMarked = T(0);
//         T totalError = std::accumulate(m_elementErrors.begin(), m_elementErrors.end(), T(0),[](T sum, const std::pair<element_t, error_t> & elem) { return sum + elem.second; });
//         T errorMarkSum = m_options.askReal("CoarsenParam",0.1) * totalError;
//         for (std::vector<std::pair<element_t, error_t>>::const_iterator it = m_elementErrors.begin(); it != m_elementErrors.end(); ++it, ++numMarked)
//         {
//             // If the cumulative error exceeds the threshold, stop the iteration
//             if (cummulErrMarked >= errorMarkSum)
//                 break;

//             // If there are no refinement elements provided, we don't check siblings
//             if (!refined.empty())
//             {
//                 // Get the siblings of the element
//                 HElementContainer siblings = m_helper.getSiblings(it->first,false);
//                 // If any of the siblings is not active, the element cannot be coarsened (since the elements in that sibling are lilely finer)
//                 if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return !m_helper.isActive(elem); }))
//                     continue;
//                 // If the element or its siblings are in the refined set, skip it
//                 if (std::any_of(siblings.begin(), siblings.end(),[&refined](const element_t & elem) { return refined.find(elem) != refined.end(); }))
//                     continue;
//             }

//             // Add the element to the container
//             result.insert(it->first);
//             cummulErrMarked += it->second;
//         }

//         if (m_options.askSwitch("Admissible",true))
//         {
//             // Loop over the marked elements
//             for (const auto & elem : result)
//             {
//                 bool erase = false;
//                 // Compute the coarsening extension
//                 HElementContainer coarseningExtension = m_helper.getCextension(elem, m_options.askInt("Jump",2));
//                 // For all elements in the coarsening extension, check if the level is the same as in the basis or finer
//                 for (const auto & coarseningElem : coarseningExtension)
//                 {
//                     if (m_helper.levelInBasis(coarseningElem) >= coarseningElem.level())
//                     {
//                         // If the level is the same or finer, the coarsening element is not admissible
//                         erase = true;
//                         break;
//                     }
//                 }

//                 if (erase)
//                 {
//                     // If the element is not admissible, erase it from the result
//                     result.erase(elem);
//                     continue;
//                 }

//                 // If the parents of any cells in the extension overlap with marked refinement cells, coarsening causes a problem
//                 for (auto & coarseningElem : coarseningExtension)
//                 {
//                     for ( const auto & refElem : refined )
//                     {
//                         if (m_helper.contains(coarseningElem, refElem))
//                         {
//                             // If the coarsening element contains a refinement element, the coarsening is not admissible
//                             erase = true;
//                             break;
//                         }
//                     }
//                     if (erase)
//                         break;
//                 }
//                 if (erase)
//                 {
//                     // If the element is not admissible, erase it from the result
//                     result.erase(elem);
//                     continue;
//                 }
//             }
//         }
//     }

// protected:
//     const gsHTensorBasis<d,T> & m_basis; ///< The basis of the elements.
//     gsHElementHelper<d,T> m_helper; ///< Helper for element operations.
//     HElementErrorContainer m_elementErrors; ///< Container for elements and their associated errors.

//     gsOptionList m_options; ///< Options for the marker.
// };


int main(int argc, char *argv[])
{
    index_t degree    = 1;
    index_t m         = 2;
    gsCmdLine cmd("Example of gsHBox.");
    cmd.addInt("m","jump",
               "parameter m", m);
    cmd.addInt("p","degree",
               "Spline degree", degree);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mpBspline, mp;

    gsTensorBSpline<2,real_t> bspline = *gsNurbsCreator<>::BSplineSquare(1,0,0);
    if (degree>1) bspline.degreeElevate(degree-1);

    // for (index_t k = 0; k < 6; ++k)
    //     bspline.uniformRefine();

    mpBspline.addPatch(bspline);

    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (size_t k=0; k!=mpBspline.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mpBspline.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.addPatch(thb);
    }

    std::vector<index_t> boxes(5);

    // Initial refinement
    boxes[0] = 1;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = 2;
    boxes[4] = 2;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 2;
    boxes[1] = 0;
    boxes[2] = 0;
    boxes[3] = 4;
    boxes[4] = 2;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 2;
    boxes[1] = 0;
    boxes[2] = 2;
    boxes[3] = 2;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 2;
    boxes[1] = 2;
    boxes[2] = 2;
    boxes[3] = 4;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 3;
    boxes[1] = 2;
    boxes[2] = 0;
    boxes[3] = 6;
    boxes[4] = 4;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 3;
    boxes[1] = 4;
    boxes[2] = 4;
    boxes[3] = 8;
    boxes[4] = 8;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 4;
    boxes[1] = 6;
    boxes[2] = 4;
    boxes[3] = 8;
    boxes[4] = 6;
    mp.patch(0).refineElements(boxes);

    boxes[0] = 4;
    boxes[1] = 10;
    boxes[2] = 12;
    boxes[3] = 12;
    boxes[4] = 14;
    mp.patch(0).refineElements(boxes);

    gsWriteParaview(mp,"init",1,true);
    gsInfo<<"Initial basis constructed:\n"<<mp.basis(0)<<"\n";

    gsHTensorBasis<2,real_t> * basis = dynamic_cast<gsHTensorBasis<2,real_t> *>(&mp.basis(0));
    gsHBoxContainer<2> markedRef, markedRef2, markedCrs;
    gsHBox<2> cell;

    gsVector<index_t,2> low,upp;
    index_t lvl;
    low <<7,5;
    upp <<8,6;
    lvl = 4;
    cell = gsHBox<2>(low,upp,lvl,basis);
    markedRef.add(cell);
    gsWriteParaview(markedRef,"markedRef");
    gsInfo<<"Added one refinement box (markedRef.pvd):\n"<<markedRef<<"\n";

    markedRef.markAdmissible(m);
    gsWriteParaview(markedRef,"refCell_Admissible");
    gsInfo<<"Admissibility region of the refined cell plotted in refCell_Admissible.pvd\n";

    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<"------------------------------------------------------------------------\n";

    typedef typename gsHElementHelper<2,real_t>::HElementContainer HElementContainer;
    typedef typename gsHElementHelper<2,real_t>::level_t           level_t;

    gsHElementHelper<2,real_t> helper(*basis);
    gsHElement<2,real_t> hcell(low,upp, 4);
    gsMatrix<> hcellMat = helper.toBox(hcell);
    gsWriteParaview(hcellMat,"hcell",(real_t)hcell.level());
    gsInfo<<"HElement created: "<<hcell<<"\n";
    cell.computeCoordinates();
    gsInfo<<"Original HBox:\n"<<cell<<"\n";

    gsHElement<2,real_t> parent = helper.getParent(hcell);
    gsMatrix<> parentMat = helper.toBox(parent);
    gsWriteParaview(parentMat,"parent",(real_t)parent.level());
    gsInfo<<"Parent of the HElement: "<<parent<<"\n";
    gsHBox<2> cellParent = cell.getParent();
    cellParent.computeCoordinates();
    gsInfo<<"Parent of the HBox:\n"<<cellParent<<"\n";

    gsHElement<2,real_t> ancestor = helper.getAncestor(hcell, 2);
    gsMatrix<> ancestorMat = helper.toBox(ancestor);
    gsWriteParaview(ancestorMat,"ancestor",(real_t)ancestor.level());
    gsInfo<<"Ancestor of the HElement (jump 2): "<<ancestor<<"\n";
    gsHBox<2> cellAncestor = cell.getAncestor(2);
    cellAncestor.computeCoordinates();
    gsInfo<<"Ancestor of the HBox (jump 2):\n"<<cellAncestor<<"\n";

    HElementContainer suppExt = helper.explode(helper.getSupportExtension(hcell),hcell.level());
    gsWriteParaview(helper.toBoxes(suppExt),"supportExtension");
    auto suppExtMat = helper.toBoxesAndLevels(suppExt);
    gsWriteParaview(suppExtMat.first,"suppExt",gsVector<real_t>(suppExtMat.second.cast<real_t>()));
    gsInfo<<"Exploded support extension of the HElement at level "<<hcell.level()<<":\n";
    for (const auto & elem : suppExt)
    {
        gsInfo<<elem<<"\n";
    }
    typename gsHBox<2>::Container suppExtCell = cell.getSupportExtension();
    gsInfo<<"Support extension of the HBox:\n";
    for (const auto & elem : suppExtCell)
    {
        elem.computeCoordinates();
        gsInfo<<elem<<"\n";
    }

    HElementContainer marked = {hcell};
    HElementContainer admissibleRegion = helper.markAdmissible(marked, m);
    auto admissibleRegionMat = helper.toBoxesAndLevels(admissibleRegion);
    gsWriteParaview(admissibleRegionMat.first,"admissibleRegion",gsVector<real_t>(admissibleRegionMat.second.cast<real_t>()));

    gsInfo<<"Marked admissible HElement:\n";
    for (const auto & elem : admissibleRegion)
    {
        gsInfo<<elem<<"\n";
    }

    std::vector<index_t> refBoxes = helper.toRefBox(hcell);
    gsInfo<<"Refinement boxes of the HElement: "<<gsAsVector<index_t>(refBoxes).transpose()<<"\n";
    std::vector<index_t> refBoxesCell = cell.toRefBox();
    gsInfo<<"Refinement boxes of the HBox: "<<gsAsVector<index_t>(refBoxesCell).transpose()<<"\n";
    gsInfo<<"Difference: "<<gsAsVector<index_t>(refBoxes).transpose() - gsAsVector<index_t>(refBoxesCell).transpose()<<"\n";

    std::vector<index_t> admissibilityRegionBoxes = helper.toRefBoxes(admissibleRegion);
    gsInfo<<"Admissibility boxes of the HElement: "<<gsAsVector<index_t>(admissibilityRegionBoxes).transpose()<<"\n";
    std::vector<index_t> admissibilityRegionBoxesCell = markedRef.toRefBoxes();
    gsInfo<<"Admissibility boxes of the HBox: "<<gsAsVector<index_t>(admissibilityRegionBoxesCell).transpose()<<"\n";
    gsInfo<<"Difference: "<<gsAsVector<index_t>(admissibilityRegionBoxes).transpose() - gsAsVector<index_t>(admissibilityRegionBoxesCell).transpose()<<"\n";

    // return 0;

    ////
    // // Get all elements in the mesh
    // HElementContainer allElements = helper.toElements();
    // gsInfo<<"All elements in the mesh:\n";
    // for (const auto & elem : allElements)
    // {
    //     gsInfo<<elem<<"\n";
    // }
    // // Make a random error vector, for testing purposes
    // gsVector<real_t> errorVec(allElements.size());
    // errorVec.setRandom();
    // std::vector<real_t> errs(errorVec.data(), errorVec.data() + errorVec.size());
    // // The gsExprEvaluator will provide such a vector in the ordering of the domain iterator. However, the elements in allElements are sorted in a different way (to make the class faster).
    // std::vector<gsHElement<2,real_t>> domainElements(basis->numElements());
    // std::vector<index_t> index(domainElements.size(),0);
    // for (const auto & domIt : basis->domain()->allElements())
    // {
    //     index_t lvl = static_cast<const gsHDomainIterator<real_t,2> *>(&domIt)->getLevel();
    //     domainElements[domIt.id()] = helper.toElement(domIt.lowerCorner(), domIt.upperCorner(), lvl, domIt.patch());
    //     index[domIt.id()] = domIt.id();
    //     gsDebug<<"Element "<<domainElements[domIt.id()]<<" with id "<<domIt.id()<<" has error "<<errs[domIt.id()]<<"\n";
    // }

    // typename gsHElement<2,real_t>::Compare comp;
    // std::sort(index.begin(), index.end(),[&domainElements,&comp](index_t a, index_t b) { return comp(domainElements[a],domainElements[b]); });

    // index_t i = 0;
    // for (auto it = allElements.begin(); it != allElements.end(); ++it, ++i)
    // {
    //     gsInfo<<"Element "<<*it<<" has error "<<errs[index[i]]<<"\n";
    // }




    return 0;

    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<"----------------------------STRESS TEST---------------------------------\n";
    gsInfo<<"------------------------------------------------------------------------\n";
    gsStopwatch timer;
    const index_t N = 100; // Number of runs
    gsInfo<<"constructors():\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHBox<2> cellTest(low,upp,lvl,basis);
    gsInfo<<"[gsHBox] Time for "<<N<<" constructors: "<<timer.stop()<<" seconds.\n";

    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHElement<2,real_t> hcellTest(low,upp, 4);
    gsInfo<<"[gsHElement] Time for "<<N<<" constructors: "<<timer.stop()<<" seconds.\n";
    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<".getParent() test:\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHBox<2> cellParentTest = cell.getParent();
    gsInfo<<"[gsHBox] Time for "<<N<<" .getParent(): "<<timer.stop()<<" seconds.\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHElement<2,real_t> hcellParentTest = helper.getParent(hcell);
    gsInfo<<"[gsHElement] Time for "<<N<<" .getParent(): "<<timer.stop()<<" seconds.\n";
    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<".getAncestor() test:\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHBox<2> cellAncestorTest = cell.getAncestor(2);
    gsInfo<<"[gsHBox] Time for "<<N<<" .getAncestor(2): "<<timer.stop()<<" seconds.\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHElement<2,real_t> hcellAncestorTest = helper.getAncestor(hcell, 2);
    gsInfo<<"[gsHElement] Time for "<<N<<" .getAncestor(2): "<<timer.stop()<<" seconds.\n";
    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<".getSupportExtension() test:\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        gsHBoxContainer<2,real_t> suppExtCellTest = cell.getSupportExtension();
    gsInfo<<"[gsHBox] Time for "<<N<<" .getSupportExtension(): "<<timer.stop()<<" seconds.\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
        HElementContainer suppExtHCellTest = helper.explode(helper.getSupportExtension(hcell),hcell.level());
    gsInfo<<"[gsHElement] Time for "<<N<<" .getSupportExtension(): "<<timer.stop()<<" seconds.\n";
    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<".markAdmissible() test:\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
    {
        gsHBoxContainer<2> markedTest(cell);
        markedTest.markTadmissible(m);
    }
    gsInfo<<"[gsHBox] Time for "<<N<<" .markAdmissible(): "<<timer.stop()<<" seconds.\n";

    timer.restart();
    for (index_t i = 0; i < N; ++i)
        HElementContainer markedAdmissibleTest = helper.markAdmissible(marked, m);
    gsInfo<<"[gsHElement] Time for "<<N<<" .markAdmissible(): "<<timer.stop()<<" seconds.\n";
    gsInfo<<"------------------------------------------------------------------------\n";
    gsInfo<<".toRefBox(.markAdmissible()) test:\n";
    timer.restart();
    for (index_t i = 0; i < N; ++i)
    {
        gsHBoxContainer<2> markedTest(cell);
        markedTest.markTadmissible(m);
        markedTest.toRefBoxes();
    }
    gsInfo<<"[gsHBox] Time for "<<N<<" .markAdmissible(): "<<timer.stop()<<" seconds.\n";

    timer.restart();
    for (index_t i = 0; i < N; ++i)
        helper.toRefBoxes(helper.markAdmissible(marked, m));
    gsInfo<<"[gsHElement] Time for "<<N<<" .toRefBox(.markAdmissible()): "<<timer.stop()<<" seconds.\n";

    gsInfo<<"-------------------------------END--------------------------------------\n";


    // low <<4,2;
    // upp <<5,3;
    // lvl = 3;
    // cell = gsHBox<2>(low,upp,lvl,basis);
    // markedCrs.add(cell);
    // gsWriteParaview(cell,"crsCell");
    // gsInfo<<"Added one coarsening box (markedCrs.pvd):\n"<<markedCrs<<"\n";

    // gsHBoxContainer<2> Cextension(cell.getParent().getCextension(m));
    // gsWriteParaview(Cextension,"crsCell_Cextension");
    // gsInfo<<"Coarsening extension plotted in crsCell_Cextension.pvd\n";

    // gsHBoxContainer<2> Cneighborhood(cell.getParent().getCneighborhood(m));
    // gsWriteParaview(Cneighborhood,"crsCell_Cneighborhood");
    // gsInfo<<"Coarsening neighborhood plotted in crsCell_Cneighborhood.pvd\n";
    return 0;
}

