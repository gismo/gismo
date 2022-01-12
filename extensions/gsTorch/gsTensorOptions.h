/** @file gsTensorOptions.h

    @brief Provides declarations of the gsTensorOptions class which is a
    wrapper of the LibTorch torch::TensorOptions class (https://pytorch.org)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#pragma once

// LibTorch has \a index_t as type. We must therefore disable it
// temporarily to prevent compiler errors.
#pragma push_macro("index_t")
#undef index_t
#include <torch/torch.h>
#pragma pop_macro("index_t")

#include <gsIO/gsXml.h>

namespace gismo {

  /// @brief Tensor options
  ///
  /// gsTensorOptions is a G+Smo wrapper for \a torch::TensorOptions
  /// and provides all of its functionality. In contrast to \a
  /// torch::TensorOptions, its data type is specified by the template
  /// parameter \a T and cannot be changed at runtime.
  template <typename T>
  class gsTensorOptions : public torch::TensorOptions
  {
  public:
    // Import constructors from base class
    using torch::TensorOptions::TensorOptions;

    // Delete constructor that modifies the data type (dtype) since
    // this is set by the compile-time parameter T and cannot be
    // modified to be consistent with G+Smo design paradigms
    gsTensorOptions(at::ScalarType) = delete;
    
    // Default constructor
    gsTensorOptions()
      : torch::TensorOptions(caffe2::TypeMeta::Make<T>()) {}

    // Copy constructor (from torch::TensorOptions)
    gsTensorOptions(const torch::TensorOptions& other)
      : torch::TensorOptions(other.dtype(caffe2::TypeMeta::Make<T>()))
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
    }

    // Move constructor (from torch::TensorOptions)
    gsTensorOptions(torch::TensorOptions&& other)
      : torch::TensorOptions(other.dtype<T>())
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
    }
    
    // Copy constructor (from gsTensorOptions)
    gsTensorOptions(const gsTensorOptions&) = default;

    // Move constructor (from gsTensorOptions)
    gsTensorOptions(gsTensorOptions&&) = default;

  public:
    // Delete dtype<T>() function from base class since the data type
    // (dtype) is set by the compile-time parameter T and cannot be
    // modified to be consistent with G+Smo design paradigms
    template<typename>
    auto dtype() = delete;

    template<typename Arg>
    auto dtype(Arg&& arg)
    { return torch::TensorOptions::dtype(std::forward<Arg>(arg)); }

    caffe2::TypeMeta dtype() const noexcept
    { return torch::TensorOptions::dtype(); }
    
  public:
    // Copy assignment operator (from torch::TensorOptions)
    gsTensorOptions& operator=(const torch::TensorOptions& other)
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
      torch::TensorOptions::operator=(other.dtype(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Move assignment operator (from torch::TensorOptions)
    gsTensorOptions& operator=(torch::TensorOptions&& other)
    {
      if (other.dtype() != caffe2::TypeMeta::Make<T>())
        gsWarn << "Implicit conversion of dtype from "
               << other.dtype() << " to "
               << caffe2::TypeMeta::Make<T>() << ".\n";
      torch::TensorOptions::operator=(other.dtype(caffe2::TypeMeta::Make<T>()));
      return *this;
    }

    // Copy assignment operator (from gsTensorOptions)
    gsTensorOptions& operator=(const gsTensorOptions& other)
    {
      torch::TensorOptions::operator=(other);
      return *this;
    }

    // Move assignment operator (from gsTensorOptions)
    gsTensorOptions& operator=(gsTensorOptions&& other)
    {
      torch::TensorOptions::operator=(other);
      return *this;
    }
    
    
  public:
    // Setter functions
    gsTensorOptions& setRequiresGrad();
    gsTensorOptions& unsetRequiresGrad();
    gsTensorOptions& setStrided();
    gsTensorOptions& setSparse();
    gsTensorOptions& setSparseCsr();
    gsTensorOptions& setCPU();
    gsTensorOptions& setCUDA(int index=0);
    gsTensorOptions& setPinnedMemory();
    gsTensorOptions& unsetPinnedMemory();
    gsTensorOptions& setMemoryFormatPreserve();
    gsTensorOptions& setMemoryFormatContiguous();
    gsTensorOptions& setMemoryFormatChannelsLast();
    gsTensorOptions& setMemoryFormatChannelsLast3d();
    gsTensorOptions& unsetMemoryFormat();

  public:
    // Comparison operators
    bool operator==(const gsTensorOptions& other) const;
    bool operator!=(const gsTensorOptions& other) const;
  };

  namespace internal
  {
    /// @brief Get a gsTensorOptions from XML data
    template<typename T>  
    class gsXml< gsTensorOptions<T> >
    {
    private:
      gsXml() { }
      typedef gsTensorOptions<T> Object;
    public:
      GSXML_COMMON_FUNCTIONS(Object);
      static std::string tag ()  { return "TensorOptions"; }
      static std::string type () { return "TensorOptions"; }
      
      GSXML_GET_POINTER(Object);
      
      static void get_into (gsXmlNode * node, Object & obj)
      {
        obj = obj.device((c10::DeviceType)atoi(node->first_attribute("device_type")->value()),
                         atoi(node->first_attribute("device_index")->value()));
        obj = obj.layout((c10::Layout)atoi(node->first_attribute("layout")->value()));
        obj = obj.pinned_memory(atoi(node->first_attribute("pinned_memory")->value()));
        obj = obj.requires_grad(atoi(node->first_attribute("requires_grad")->value()));
      }
      
      static gsXmlNode * put (const Object & obj, gsXmlTree & data )
      {
        gsXmlNode * node = makeNode("TensorOptions", data);

        node->append_attribute( makeAttribute("device_type", util::to_string((int)obj.device().type()), data) );
        node->append_attribute( makeAttribute("device_index", util::to_string((int)obj.device().index()), data) );
        node->append_attribute( makeAttribute("layout", util::to_string((int)obj.layout()), data) );
        node->append_attribute( makeAttribute("pinned_memory", util::to_string(obj.pinned_memory()), data) );
        node->append_attribute( makeAttribute("requires_grad", util::to_string(obj.requires_grad()), data) );
        
        return node;
      }
    };
  } // namespace internal
  
} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorOptions.hpp)
#endif

