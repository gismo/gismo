/** @file gsIOUtils.h

    @brief Input and output Utilities
    
    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Speh, J. Zwar
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsHSplines/gsHDomainIterator.h>
#include <gsHSplines/gsHTensorBasis.h>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace gismo {

/// \brief Returns the computational mesh of \a basis.
///
/// \param[in] basis
/// \param mesh
/// \param[in] n number of samples per element side
template <class T>
void makeMesh(const gsBasis<T>& basis, gsMesh<T>& mesh, int n = 0) {
  const unsigned d = basis.dim();

  typedef typename gsMesh<T>::VertexHandle Vertex;
  typename gsBasis<T>::domainIter domIter = basis.makeDomainIterator();

  // variables for iterating over a cube (element is a cube)
  const gsVector<unsigned> zeros = gsVector<unsigned>::Zero(d);
  const gsVector<unsigned> ones = gsVector<unsigned>::Ones(d);
  gsVector<unsigned> cur;

  // maps integer representation of a vertex into pointer to the
  // vertex coordinates
  std::vector<Vertex> map(1ULL << d);

  // neighbour[i] are integer representations of certain neighbours of
  // vertex i (i counts in lexicographics order over all vertices)
  std::vector<std::vector<unsigned> > neighbour(1ULL << d,
                                                std::vector<unsigned>());

  cur.setZero(d);
  int counter = 0;
  do {
    // set neighbour
    for (unsigned dim = 0; dim < d; dim++) {
      if (cur(dim) == 0) {
        const unsigned tmp = counter | (1 << dim);
        neighbour[counter].push_back(tmp);
      }
    }
    counter++;

  } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));

  gsVector<T> vertex(d);

  for (; domIter->good(); domIter->next()) {
    const gsVector<T>& low = domIter->lowerCorner();
    const gsVector<T>& upp = domIter->upperCorner();
    const T vol = domIter->volume();

    vertex.setZero();
    cur.setZero();
    counter = 0;

    // add points to the mesh
    do {
      // get appropriate coordinate of a point
      for (unsigned dim = 0; dim < d; dim++) {
        vertex(dim) = (cur(dim) ? upp(dim) : low(dim));
      }

      Vertex v = mesh.addVertex(vertex);
      v->data = vol;
      map[counter++] = v;

    } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));

    // add edges to the mesh (connect points)
    for (size_t index = 0; index != neighbour.size(); index++) {
      const std::vector<unsigned>& v = neighbour[index];

      for (size_t ngh = 0; ngh != v.size(); ngh++) {
        // add more vertices (n) for better physical resolution
        mesh.addLine(map[index], map[v[ngh]], n);
        // mesh.addEdge( map[index], map[v[ngh]] );
      }
    }

    // idea: instead of edges add the faces to the mesh
    // mesh->addFace( mesh.vertices().back(),
    //                *(mesh.vertices().end()-3),
    //                *(mesh.vertices.end()-4),
    //                *(mesh.vertices.end()-2)
    //     );
  }
}

namespace internal {

/// Look at function gismo::makeHierarchicalMesh
template <short_t d, typename T>
void makeHierarchicalMesh(const gsHTensorBasis<d, T>& basis,
                          std::vector<gsMesh<T> >& meshes, int n = 0) {
  // prepare meshes
  meshes.clear();
  for (unsigned i = 0; i < basis.maxLevel() + 1; i++) {
    meshes.push_back(gsMesh<T>());
  }

  // variables for iterating over a cube (element is a cube)
  const gsVector<unsigned> zeros = gsVector<unsigned>::Zero(d);
  const gsVector<unsigned> ones = gsVector<unsigned>::Ones(d);
  gsVector<unsigned> cur;

  // neighbour[i] are integer representations of certain neighbours of
  // vertex i (i counts in lexicographics order over all vertices)
  std::vector<std::vector<unsigned> > neighbour(1 << d,
                                                std::vector<unsigned>());

  cur.setZero(d);
  int counter = 0;
  do {
    // set neighbour
    for (unsigned dim = 0; dim < d; dim++) {
      if (cur(dim) == 0) {
        const unsigned tmp = counter | (1 << dim);
        neighbour[counter].push_back(tmp);
      }
    }
    counter++;

  } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));

  // maps integer representation of a vertex into pointer to the
  // vertex coordinates
  typedef typename gsMesh<T>::VertexHandle Vertex;
  std::vector<Vertex> map(1 << d);

  gsVector<T> vertex(d);

  gsHDomainIterator<T, d> domIter(basis);

  for (; domIter.good(); domIter.next()) {
    int level = domIter.getLevel();

    const gsVector<T>& low = domIter.lowerCorner();
    const gsVector<T>& upp = domIter.upperCorner();

    vertex.setZero();
    cur.setZero();
    counter = 0;

    // add points to the mesh
    do {
      // get appropriate coordinate of a point
      for (unsigned dim = 0; dim < d; dim++) {
        vertex(dim) = (cur(dim) ? upp(dim) : low(dim));
      }

      meshes[level].addVertex(vertex);
      map[counter++] = meshes[level].vertices().back();

    } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));

    // add edges to the mesh (connect points)
    for (size_t index = 0; index != neighbour.size(); index++) {
      const std::vector<unsigned>& v = neighbour[index];

      for (size_t ngh = 0; ngh != v.size(); ngh++) {
        // add more vertices (n) for better physical resolution
        meshes[level].addLine(map[index], map[v[ngh]], n);
      }
    }
  }
}

}  // end namespace internal

/// Constructs a series of meshes, each mesh presents on level in hierarchical
/// basis.
///
/// \param basis hierarchocal tensor basis
/// \param meshes we return meshes via reference to std::vector
/// \param n how many vertices is used in presentation of one element
///
/// \result success of the function
template <typename T>
bool makeHierarchicalMesh(const gsBasis<T>& basis,
                          std::vector<gsMesh<T> >& meshes, int n = 0) {
  const gsHTensorBasis<1, T>* hBasis1 =
      dynamic_cast<const gsHTensorBasis<1, T>*>(&basis);

  if (hBasis1 != NULL) {
    internal::makeHierarchicalMesh<1, T>(*hBasis1, meshes, n);
    return true;
  }

  const gsHTensorBasis<2, T>* hBasis2 =
      dynamic_cast<const gsHTensorBasis<2, T>*>(&basis);

  if (hBasis2 != NULL) {
    internal::makeHierarchicalMesh<2, T>(*hBasis2, meshes, n);
    return true;
  }

  const gsHTensorBasis<3, T>* hBasis3 =
      dynamic_cast<const gsHTensorBasis<3, T>*>(&basis);

  if (hBasis3 != NULL) {
    internal::makeHierarchicalMesh<3, T>(*hBasis3, meshes, n);
    return true;
  }

  return false;
}

/**
 * @brief Encode for base64 export
 *
 * Static class to provide functionality to encode any type of vector into an
 * ascii string (uncompressed) and its inverse operation.
 *
 */
class Base64 {
 private:
  /// Alias for one byte type
  using ByteRepresentation = unsigned char;

  /// Look up table
  static const char char_encode_table(const unsigned& index) {
    // Create static array in function to avoid use of c++17 static member
    // declaration
    static const std::array<char, 64> encode_table{
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'};
    GISMO_ASSERT((index < 64), "Requested index out of range. Input invalid");
    return encode_table[index];
  }

  /**
   * @brief Reverse the encoding table in an array
   *
   * Chars are casted to unsigned types for indexing, function will only be
   * called once
   *
   * @return std::array<unsigned, 256>
   */
  static const std::array<unsigned, 256> ReverseCharEncodeTable_() {
    std::array<unsigned, 256> et_reversed{};
    // Fill remainding array entries with values outside range to throw
    // assertions easily
    et_reversed.fill(256);
    const unsigned n_chars_{64};
    for (unsigned i{}; i < n_chars_; i++) {
      et_reversed[static_cast<unsigned>(char_encode_table(i))] = i;
    }
    return et_reversed;
  }

  /// Lookup Table for Decoding B64 string
  static const unsigned char_decode_table(const unsigned& index) {
    static const std::array<unsigned, 256> decode_table{
        ReverseCharEncodeTable_()};
    GISMO_ASSERT((index < 256), "Requested index out of range. Input invalid");
    GISMO_ASSERT((encode_table[index] == 256),
                 "Invalid decode type, this should never occur!");
    return decode_table[index];
  }

  // Check if string fulfills minimum requirements for B64 encoded strings
  static bool isValidBase64String(const std::string& s) {
    // Check if size is feasible
    bool is_valid{s.size() % 4 == 0};
    // Check if string only contains alpha-numeric chars or `/`, `+` or `=`
    return is_valid && std::all_of(s.begin(), s.end(), [](const char& c) {
             return isalnum(c) || c == '+' || c == '/' || c == '=';
           });
  }

 public:
  /**
   * @brief Actual encoding routine
   *
   * @tparam BaseType type of individual data entries
   * @param data_vector data to be encoded
   * @return std::string encoded data
   */
  template <typename BaseType>
  static std::string Encode(const std::vector<BaseType>& data_vector) {
    const ByteRepresentation* vector_as_bytes =
        reinterpret_cast<const ByteRepresentation*>(&data_vector[0]);

    // Number of bytes for an entry
    constexpr const std::size_t length_of_entry{sizeof(BaseType{})};
    // Minimum number of bytes required
    const std::size_t minimum_n_bytes_required =
        length_of_entry * data_vector.size();
    // Number of bytes must be divisible by three
    const std::size_t additional_padding_bytes =
        (3 - minimum_n_bytes_required % 3) % 3;
    // Required groups of three
    const std::size_t number_of_groups =
        (minimum_n_bytes_required + additional_padding_bytes) / 3;

    // Initialize return value
    std::string encoded_string;
    encoded_string.reserve(number_of_groups * 4);

    // Loop over bytes and decode them
    for (std::size_t i_group{}; i_group < number_of_groups; i_group++) {
      const std::size_t buffer_index = i_group * 3;
      std::array<ByteRepresentation, 3> buffer{};
      buffer[0] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 0]
                      : 0;
      buffer[1] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 1]
                      : 0;
      buffer[2] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 2]
                      : 0;

      // Transform bytes into chars using above private encoder table
      encoded_string.push_back(char_encode_table(((buffer[0] & 0xfc) >> 2)));
      encoded_string.push_back(char_encode_table(((buffer[0] & 0x03) << 4) +
                                                 ((buffer[1] & 0xf0) >> 4)));
      encoded_string.push_back(char_encode_table(((buffer[1] & 0x0f) << 2) +
                                                 ((buffer[2] & 0xc0) >> 6)));
      encoded_string.push_back(char_encode_table(((buffer[2] & 0x3f) << 0)));
    }

    // Replace trailing invalid data with `=`
    for (size_t i = 0; i < additional_padding_bytes; ++i) {
      encoded_string.push_back('=');
    }

    return encoded_string;
  }

  /**
   * @brief Reading a B64 string, transforming it into a vector of a specific
   * type
   *
   * @tparam OutputType target type
   */
  template <typename OutputType>
  static std::vector<OutputType> Decode(const std::string& base64string) {
    // Check validity of string
    GISMO_ASSERT(isValidBase64String(base64string.size),
                 "Validity check failed");

    // Init return value
    const std::size_t number_of_groups{base64string.size() / 4};
    constexpr const std::size_t length_of_entry{sizeof(OutputType{})};
    const std::size_t number_of_output_values{(number_of_groups * 3) /
                                              length_of_entry};
    std::vector<OutputType> return_value;
    return_value.resize(number_of_output_values);

    // Access as byte stream
    ByteRepresentation* vector_as_bytes =
        reinterpret_cast<ByteRepresentation*>(&return_value[0]);

    // Start the reverse process
    for (std::size_t i_group{}; i_group < number_of_groups; i_group++) {
      const std::size_t buffer_index = i_group * 4;
      std::array<unsigned, 4> buffer{};
      for (unsigned i{}; i < 4; i++) {
        buffer[i] = base64string[buffer_index + i] != '='
                        ? char_decode_table(static_cast<unsigned>(
                              base64string[buffer_index + i]))
                        : 255;
      }

      // Write bytes into vector
      if (buffer[1] != 255) {
        vector_as_bytes[i_group * 3] =
            ((buffer[0] & 0x3f) << 2) + ((buffer[1] & 0x30) >> 4);
      }
      if (buffer[2] != 255) {
        vector_as_bytes[i_group * 3 + 1] =
            ((buffer[1] & 0x0f) << 4) + ((buffer[2] & 0x3c) >> 2);
      }
      if (buffer[3] != 255) {
        vector_as_bytes[i_group * 3 + 2] =
            ((buffer[2] & 0x03) << 6) + ((buffer[3] & 0x3f) >> 0);
      }
    }
    return return_value;
  }
};

}  // namespace gismo
