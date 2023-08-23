/** @file gsXmlGenericUtils.hpp

    @brief Provides implementation of generic XML functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConstantFunction.h>

namespace gismo {

namespace internal {

template<class Object>
Object * getTensorBasisFromXml ( gsXmlNode * node)
{
    // Note: we do not check for the Object tag here, to allow some
    // freedom (eg. BSplineBasis instead of TensorBSplineBasis1
    GISMO_ASSERT( !strcmp( node->name(),"Basis"), "Invalid tag when \"Basis\" was expected.");
    
	// Component container
	std::vector<typename Object::CoordinateBasis* > bb;

    // Special case of reading a 1D tensor basis
    if ( !strcmp(node->first_attribute("type")->value(),
                 gsXml<typename Object::CoordinateBasis>::type().c_str() ) ) // to do in derived
    {
        bb.push_back( gsXml<typename Object::CoordinateBasis>::get(node) );
        return Object::New( bb );
    }

    const unsigned d = Object::Dim;
    gsXmlNode * tmp = node->first_node("Basis");
    if ( tmp )
    {
        for ( unsigned i = 0; i!=d; ++i)
        {
            bb.push_back( gsXml<typename Object::CoordinateBasis>::get(tmp) );
            tmp =  tmp->next_sibling("Basis");        
        }
    }
    else
    {
        GISMO_ASSERT( d == 1, "Wrong data in the xml file.");
        bb.push_back( gsXml<typename Object::CoordinateBasis>::get(node) );
    }

    //gsDebugVar( bb.size() );
    return Object::New( bb );
}

template<class Object>
gsXmlNode * putTensorBasisToXml ( Object const & obj, gsXmlTree & data)
{
    // Write the component bases
    static const unsigned d = Object::Dim;
    if (d==1)
        return gsXml<typename Object::CoordinateBasis>::put(obj.component(0), data);

    // Add a new node (without data)
    gsXmlNode* tp_node = internal::makeNode("Basis" , data);        
    tp_node->append_attribute( makeAttribute("type", 
                                             internal::gsXml<Object>::type().c_str(), data) );
    
	for ( unsigned i = 0; i!=d; ++i )
    {
		gsXmlNode* tmp = 
		    internal::gsXml< typename Object::CoordinateBasis >::put(obj.component(i), data );
		tmp->append_attribute( makeAttribute("index", i, data) );
		tp_node->append_node(tmp);
    }
    
    // All set, return the basis
    return tp_node;
}

template<class Object>
Object * getRationalBasisFromXml ( gsXmlNode * node)
{
    GISMO_ASSERT( ( !strcmp( node->name(),"Basis") )
                  &&  ( !strcmp(node->first_attribute("type")->value(),
                                internal::gsXml<Object>::type().c_str() ) ),
                  "Something is wrong with the XML data: There should be a node with a "
                  <<internal::gsXml<Object>::type().c_str()<<" Basis.");

    // Read source basis
    gsXmlNode * tmp = node->first_node("Basis");
    typename Object::SourceBasis * src = gsXml<typename Object::SourceBasis>::get(tmp) ;

    // Read weights
    tmp = node->first_node("weights");
    gsMatrix<typename Object::Scalar_t> weights;
    gsXmlAttribute* format = tmp->first_attribute("format");
    std::string format_flag = format ? format->value() : "ASCII";
    getMatrixFromXml<typename Object::Scalar_t>(tmp, src->size(), 1, weights,
                                                format_flag);
    return new Object(src, give(weights));
}

template<class Object>
Object * getHTensorBasisFromXml ( gsXmlNode * node)
{
    GISMO_ASSERT( ( !strcmp( node->name(),"Basis") )
                  &&  ( !strcmp(node->first_attribute("type")->value(),
                                internal::gsXml<Object>::type().c_str() ) ),
                  "Something is wrong with the XML data: There should be a node with a "<<
                  internal::gsXml<Object>::type().c_str()<<" Basis.");

    typedef typename Object::Scalar_t T;
    static const int d = Object::Dim;
    
    // Read max level
    //unsigned lvl = atoi( node->first_attribute("levels")->value() );
    gsXmlNode * tmp = node->first_node("Basis");
    GISMO_ASSERT( tmp , "Expected to find a basis node.");
    
    // Read the Tensor-product basis
    gsTensorBSplineBasis<d,T> * tp = 
        gsXml<gsTensorBSplineBasis<d,T> >::get(tmp);
    
    // Initialize the HBSplineBasis
    std::istringstream str;
    
    // Insert all boxes
    unsigned c;
    std::vector<index_t> all_boxes;
    for (tmp = node->first_node("box"); 
         tmp; tmp = tmp->next_sibling("box"))
    {
        all_boxes.push_back(atoi( tmp->first_attribute("level")->value() ));
        str.clear();
        str.str( tmp->value() );
        for( unsigned i = 0; i < 2*d; i++)
        {
            str>> c;
            all_boxes.push_back(c);
        }
    }

    gsXmlAttribute * manualLevels = node->first_attribute("manualLevels");
    bool ml = manualLevels && !strcmp(manualLevels->value(),"true");
    Object * hbs = new Object(*tp, ml);
    delete tp;

    if (ml)
    {
        index_t lvl = 1;
        const gsXmlAttribute * id_at;
        for (gsXmlNode * child = node->first_node("Basis");
             child; child = child->next_sibling("Basis"))
        {
            id_at = child->first_attribute("level");
            if (id_at && atoi(id_at->value()) == lvl )
            {
                ++lvl;
                auto tb = memory::make_unique(
                    internal::gsXml<gsTensorBSplineBasis<d,T> >::get(child) );
                hbs->addLevel( give(*tb) );
            }
        }
    }

    hbs->refineElements(all_boxes);
    return hbs;
}

template<class Object>
gsXmlNode * putHTensorBasisToXml ( Object const & obj, gsXmlTree & data)
{
    //typedef typename Object::Scalar_t T;
    const int d = obj.dim();

    // Add a new node (without data)
    gsXmlNode* tp_node = internal::makeNode("Basis" , data);
    
    tp_node->append_attribute( makeAttribute("type",
                                             internal::gsXml<Object>::type().c_str(), data) );

    //tp_node->append_attribute( makeAttribute( "levels",2 ,data )); // deprecated
  
    gsXmlNode * tmp;
    if (obj.manualLevels())
    {
        tp_node->append_attribute( makeAttribute("manualLevels","true", data) );
        for (index_t l = 0; l != obj.numLevels(); l++)
        {
            tmp = putTensorBasisToXml(obj.tensorLevel(l), data);
            tmp->append_attribute( makeAttribute("level", to_string(l), data ) );
            tp_node->append_node(tmp);
        }
    }
    else
    {
        tp_node->append_attribute( makeAttribute("manualLevels","false", data) );
        // Write the component bases
        tmp = putTensorBasisToXml(obj.tensorLevel(0), data);
        tp_node->append_node(tmp);
    }

    
    //Output boxes
    gsMatrix<index_t> box(1,2*d);

    for( typename Object::hdomain_type::const_literator lIter = 
             obj.tree().beginLeafIterator(); lIter.good() ; lIter.next() )
    {
        if ( lIter->level > 0 )
        {
            box.leftCols(d)  = lIter.lowerCorner().transpose();
            box.rightCols(d) = lIter.upperCorner().transpose();
       
            tmp = putMatrixToXml( box, data, "box" );
           
            tmp->append_attribute( makeAttribute("level", to_string(lIter->level), data ) );
            tp_node->append_node(tmp);
        }
    }

/*
// Write box history (deprecated)
typename Object::boxHistory const & boxes = obj.get_inserted_boxes();

for(unsigned int i = 0; i < boxes.size(); i++)
{
box.leftCols(d)  = boxes[i].lower.transpose();
box.rightCols(d) = boxes[i].upper.transpose();

tmp = putMatrixToXml( box, data, "box" );
tmp->append_attribute( makeAttribute("level", to_string(boxes[i].level), data ) );
tp_node->append_node(tmp);
}
*/

    // All set, return the basis
    return tp_node;
}


template<class Object>
gsXmlNode * putRationalBasisToXml ( Object const & obj, gsXmlTree & data)
{
    // Add a new node
    gsXmlNode* rat_node = internal::makeNode("Basis" , data);        
    rat_node->append_attribute( makeAttribute("type",
                                              internal::gsXml< Object >::type().c_str(), data) );
	
    // Write the source basis
	gsXmlNode* tmp = 
        internal::gsXml< typename Object::SourceBasis >::put(obj.source(), data );
    rat_node->append_node(tmp);
    
    // Write the weights
    tmp = putMatrixToXml( obj.weights(), data, "weights" );
    rat_node->append_node(tmp);
	
	// All done, return the node
	return rat_node;
}




/*
template<class Object>
Object * getById(gsXmlNode * node, const int & id)
{
    std::string tag = internal::gsXml<Object>::tag();
    for (gsXmlNode * child = node->first_node(tag.c_str());
         child; child = child->next_sibling(tag.c_str()))
    {
        if (  atoi(child->first_attribute("id")->value() ) == id )
            return internal::gsXml<Object>::get(child);
    }
    std::cerr<<"gsXmlUtils Warning: "<< internal::gsXml<Object>::tag() 
             <<" with id="<<id<<" not found.\n";
    return NULL;
}
*/


// Helper to get the tag of an id
// std::string getTag(gsXmlNode * node, const int & id)
// {
//     for (gsXmlNode * child = node->first_node(); 
// 	 child; child = child->next_sibling() )
// 	if (  atoi(child->first_attribute("id")->value() ) == id )
// 	    return child->name();
    
//     std::cerr<<"gsXmlUtils Warning: Tag with id="<<id<<" not found.\n";
//     return "";
// }




/// Helper to fetch geometries
//template<class Object> 
//Object * getGeometryFromXml ( gsXmlNode * node);
template<class Object>
Object * getGeometryFromXml ( gsXmlNode * node)
{
    //gsWarn<<"Reading "<< gsXml<Object>::type() <<" Geometry..\n";
    assert ( ( !strcmp( node->name(),"Geometry") ) &&
             ( !strcmp(node->first_attribute("type")->value(), gsXml<Object>::type().c_str() ) ) );

    gsXmlNode * tmp = node->first_node("Basis");

    // to do: avoid copy object here (remove Ptr)
    typename Object::Basis::Ptr b( gsXml<typename Object::Basis>::get(tmp) );

    //gsWarn<<"Read basis from node "<< tmp <<", got "<< *b <<"\n";

    tmp = node->first_node("coefs");
    GISMO_ASSERT( tmp, "Did not find any coefficients for "<< gsXml<Object>::type().c_str() );
    gsXmlAttribute * at_geodim = tmp->first_attribute("geoDim");
    GISMO_ASSERT( at_geodim , "geoDim attribute not found in Geometry XML tag");
    unsigned geoDim = atoi(at_geodim->value() ) ;

    //gsWarn<<"Read mat "<< b->size()<<"x"<< geoDim <<"\n";

    // Read the Coefficients and store them in a matrix
    gsMatrix<typename Object::Scalar_t> coefficient_matrix;
    gsXmlAttribute* format = tmp->first_attribute("format");
    std::string format_flag = format ? format->value() : "ASCII";
    getMatrixFromXml<typename Object::Scalar_t>(
        tmp, b->size(), geoDim, coefficient_matrix, format_flag);


    gsXmlAttribute* coef_order = tmp->first_attribute("order");
    if (nullptr != coef_order)
        if (!strcmp(coef_order->value(), "coordinates")) {
            coefficient_matrix.transposeInPlace();
            coefficient_matrix.resize(b->size(), geoDim);
        }

    // Looking for transformations
    tmp = node->first_node("transform");
    gsMatrix<typename Object::Scalar_t> a;
    if ( tmp )
    {
        for (gsXmlNode * tr = tmp->first_node();  
             tr; tr= tr->next_sibling() )
        {
            std::string val( tr->name() );

            if (val == "translation")
            {
                getMatrixFromXml<typename Object::Scalar_t>(tmp, 3, 1 ,a);
                // coefficient_matrix->rowwise() += a->transpose(); // TO DO
            }
            if (val ==  "rotation" ) // 3d
            {
                getMatrixFromXml<typename Object::Scalar_t>(tmp, 4, 1, a);
                gsEigen::Transform<typename Object::Scalar_t,3,gsEigen::Affine> 
                    rot( gsEigen::AngleAxis<typename Object::Scalar_t> 
                         ( a(3,0), a.template block<3,1>(0,0).normalized() ) );
                coefficient_matrix = (coefficient_matrix. rowwise().homogeneous() * 
                     rot.matrix().transpose() ).leftCols(3) ;

            }
            if (val == "scale")
            {
            
            }
            else
            {
                gsWarn<< "Unidentified transform tag in XML.\n";
            }
        }
    }

    Object * result = new Object(*b,coefficient_matrix);
    return result;
}

/// Helper to put geometries to XML
//template<class Object>
//gsXmlNode * putGeometryToXml ( Object const & obj, gsXmlTree & data);
template<class Object>
gsXmlNode * putGeometryToXml ( Object const & obj, gsXmlTree & data)
{
    // Make a new XML Geometry node 
    gsXmlNode * bs = internal::makeNode("Geometry", data);
    bs->append_attribute( makeAttribute("type", 
                                        internal::gsXml<Object>::type().c_str(), data) );

    // Add the basis
    gsXmlNode* tmp = 
	    internal::gsXml< typename Object::Basis >::put(obj.basis(), data);        
	if ( ! tmp )
    {
	    gsWarn<<"XML Warning: Writing basis failed.\n";
	    return NULL;
    }

    bs->append_node(tmp);

    // Write the coefficient matrix
    tmp = putMatrixToXml( obj.coefs(), data, "coefs" );
    tmp->append_attribute( makeAttribute("geoDim", obj.geoDim(), data) );
    bs->append_node(tmp);
    return bs;
}

template < class T >
gsXmlNode * putFunctionExprToXml(const gsFunctionExpr<T> & obj, gsXmlNode * result, gsXmlTree & data)
{
    std::string typeStr = gsXml<gsFunctionExpr<T> >::type();
    gsXmlAttribute * type = internal::makeAttribute("type", typeStr, data);
    result->append_attribute(type);
    gsXmlAttribute * dim = internal::makeAttribute("dim", obj.domainDim(),
            data);
    result->append_attribute(dim);
    // set value
    const short_t tdim = obj.targetDim();
    if ( tdim == 1)
    {
        result->value( makeValue(obj.expression(), data) );
    }
    else
    {
        gsXmlNode * cnode;
        for (short_t c = 0; c!=tdim; ++c)
        {
            cnode = makeNode("c", obj.expression(c), data);
            result->append_node(cnode);
        }
    }
    return result;
}

template < class T >
gsXmlNode * putConstantFunctionToXml(const gsConstantFunction<T> & obj,gsXmlNode * result, gsXmlTree & data)
{
    std::string typeStr = "FunctionExpr";
    gsXmlAttribute * type = internal::makeAttribute("type", typeStr, data);
    result->append_attribute(type);
    gsXmlAttribute * dim = internal::makeAttribute("dim", obj.domainDim(),
            data);
    result->append_attribute(dim);

    // set value
    gsMatrix<T> value = obj.value();
    result->value( makeValue( value, data, true) );
    return result;
}

template < class T >
gsXmlNode * putFunctionToXml ( const typename gsFunctionSet<T>::Ptr & obj, gsXmlTree & data, int index)
{
    gsXmlNode * result = internal::makeNode("Function", data);
    if (typeid(*obj) == typeid(gsFunctionExpr<T> ))
    {
        gsFunctionExpr<T> * ptr2 =
                dynamic_cast<gsFunctionExpr<T> *>(obj.get());
        result = putFunctionExprToXml(*ptr2, result, data);
    }
    else if (typeid(*obj) == typeid(gsConstantFunction<T> ))
        {
        gsConstantFunction<T> * ptr2 =
                dynamic_cast<gsConstantFunction<T> *>(obj.get());
        result = putConstantFunctionToXml(*ptr2, result, data);
    }
    gsXmlAttribute * indexNode = internal::makeAttribute("index", index,
            data);
    result->append_attribute(indexNode);
    return result;
}





}// end namespace internal

}// end namespace gismo
