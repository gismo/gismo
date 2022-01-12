/** @file gsTensor.h

    @brief Provides declarations of the gsTensor class which is a
    wrapper of the LibTorch torch::Tensor class (https://pytorch.org)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

#include <gsTorch/gsTensorOptions.h>

namespace gismo {

  template <typename T>
  class gsTensor : public torch::Tensor
  {
  public:
    // Import constructors from base class
    using torch::Tensor::Tensor;
    
    // Default constructor
    gsTensor()
      : torch::Tensor(torch::empty({}, caffe2::TypeMeta::Make<T>())) {}

    // Copy constructor (from torch::Tensor)
    gsTensor(const torch::Tensor& other)
      : torch::Tensor(other.to(caffe2::TypeMeta::Make<T>()))
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
    }
    
    // Move constructor (from torch::Tensor)
    gsTensor(torch::Tensor&& other)
      : torch::Tensor(other.to(caffe2::TypeMeta::Make<T>()))
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
    }

    // Copy constructor (from gsTensor)
    gsTensor(const gsTensor&) = default;

    // Move constructor (from gsTensor)
    gsTensor(gsTensor&&) = default;

    // Constructor from dimension array
    gsTensor(at::IntArrayRef size)
      : torch::Tensor(torch::empty(size, caffe2::TypeMeta::Make<T>())) {}
    
  public:
    // Copy assignment operator (from torch::Tensor)
    gsTensor& operator=(const torch::Tensor& other)
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
      torch::Tensor::operator=(other.to(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Move assignment operator (from torch::Tensor)
    gsTensor& operator=(torch::Tensor&& other)
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
      torch::Tensor::operator=(other.to(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Copy assignment operator (from gsTensor)
    gsTensor& operator=(const gsTensor& other)
    {
      torch::Tensor::operator=(other);
      return *this;
    }

    // Move assignment operator (from gsTensor)
    gsTensor& operator=(gsTensor&& other)
    {
      torch::Tensor::operator=(other);
      return *this;
    }
    
  public:
    // Query functions
    bool isActive() const;
    bool isPassive() const;
    bool isStrided() const;
    bool isSparse() const;
    bool isSparseCsr() const;
    bool isCPU() const;
    bool isCUDA(int index=0) const;
    bool isPinnedMemory() const;
    bool isMemoryFormatPreserve() const;
    bool isMemoryFormatContiguous() const;
    bool isMemoryFormatChannelsLast() const;
    bool isMemoryFormatChannelsLast3d() const;

    // Setter functions
    gsTensor& setRequiresGrad();
    gsTensor& unsetRequiresGrad();
  };

  namespace internal
  {
    /// @brief Get a gsTensor from XML data
    template<typename T>  
    class gsXml< gsTensor<T> >
    {
    private:
      gsXml() { }
      typedef gsTensor<T> Object;
    public:
      GSXML_COMMON_FUNCTIONS(Object);
      static std::string tag ()  { return "Tensor"; }
      static std::string type () { return "Tensor"; }
      
      GSXML_GET_POINTER(Object);
      
      static void get_into (gsXmlNode * node, Object & obj)
      {
        gsXmlAttribute * attrib = node->first_attribute("binary");

        if (attrib != NULL && attrib->value()) {

          gsXmlNode * data = node->first_node("base64");

          if (data != NULL) {
            std::stringstream sstream;
            sstream << util::from_base64(data->value());
            torch::load(obj, sstream);
          } else {
            GISMO_ERROR("XML object does not provide element \"base64\".");
          }
        } else {
          GISMO_ERROR("XML object does not provide attribute \"binary\".");
        }                        
      }
      
      static gsXmlNode * put (const Object & obj, gsXmlTree & data )
      {
        gsXmlNode * node = makeNode("Tensor", data);

        // Binary encoding
        node->append_attribute( makeAttribute("binary", true, data) );
        
        std::stringstream sstream;
        torch::save(obj, sstream);
        std::string str = util::to_base64(sstream.str());
        
        node->append_node( data.allocate_node( rapidxml::node_element,
                                               data.allocate_string("base64"),
                                               data.allocate_string(str.c_str(),
                                                                    str.size()) ) );        
        return node;
      }
    };
  } // namespace internal
  
} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensor.hpp)
#endif
