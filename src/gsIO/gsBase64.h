#pragma once

#include <gsCore/gsDebug.h>
#include <gsMatrix/gsMatrix.h>

#include <algorithm>
#include <array>
#include <numeric>
#include <string>
#include <vector>

namespace gismo {

/**
 * @brief Encode for base64 export
 *
 * Static class to provide functionality to encode any type of vector into an
 * ascii string (uncompressed) and its inverse operation.
 *
 */
class GISMO_EXPORT Base64 {
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
        GISMO_ASSERT((index < 64),
                     "Requested index out of range. Input invalid");
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
        GISMO_ASSERT((index < 256),
                     "Requested index out of range. Input invalid");
        GISMO_ASSERT((decode_table[index] != 256),
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

    /**
     * @brief Trim trailing and preceding whitespaces
     *
     * @param s string to be processed
     * @return std::string
     */
    static std::string trimWhitespaces(const std::string& s) {
        const std::string delimiters(" \n\t");
        size_t first = s.find_first_not_of(delimiters);
        if (std::string::npos == first) {
            GISMO_ERROR("Empty string cannot be converted into data-vector");
        }
        size_t last = s.find_last_not_of(delimiters);
        return s.substr(first, (last - first + 1));
    }

    /**
     * @brief Actual encoding routine where byte-stream is transformed and
     * encoded
     *
     * @param byte_vector_ptr pointer to start of an array (contiguous)
     * @param minimum_n_bytes_required number_of_bytes in byte-stream
     * @return std::string encoded data
     */
    static std::string Encode_(const ByteRepresentation* byte_vector_ptr,
                               const std::size_t& minimum_n_bytes_required) {
        // Padding blocks
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
                            ? byte_vector_ptr[buffer_index + 0]
                            : 0;
            buffer[1] = buffer_index < minimum_n_bytes_required
                            ? byte_vector_ptr[buffer_index + 1]
                            : 0;
            buffer[2] = buffer_index < minimum_n_bytes_required
                            ? byte_vector_ptr[buffer_index + 2]
                            : 0;

            // Transform bytes into chars using above private encoder table
            encoded_string.push_back(
                char_encode_table(((buffer[0] & 0xfc) >> 2)));
            encoded_string.push_back(char_encode_table(
                ((buffer[0] & 0x03) << 4) + ((buffer[1] & 0xf0) >> 4)));
            encoded_string.push_back(char_encode_table(
                ((buffer[1] & 0x0f) << 2) + ((buffer[2] & 0xc0) >> 6)));
            encoded_string.push_back(
                char_encode_table(((buffer[2] & 0x3f) << 0)));
        }

        // Replace trailing invalid data with `=`
        for (size_t i = 0; i < additional_padding_bytes; ++i) {
            encoded_string[number_of_groups * 4 - i - 1] = '=';
        }

        // Safety check
        GISMO_ASSERT(
            isValidBase64String(encoded_string),
            "Something went wrong in B64 encoding, please write an issue");
        return encoded_string;
    }

    /**
     * @brief Copy a read vector of type BaseType into a gsMatrix with
     * ScalarType TargetType.
     *
     * This is required as they might differ from each other, e.g., double into
     * float, uint into int, etc.. Further, most matrices in Gismo (like
     * coefficients) are stored in a different order (colwise) than the input
     * stream, which prohibits the use of writing directly on the pointer
     *
     * @tparam BaseType Type as encoded into the input file
     * @tparam TargetType Type as requested from target
     * @param base_vector decoded std::vector with base type
     * @param result gsMatrix passed as reference
     */
    template <typename BaseType, typename TargetType>
    static void CopyIntoGsMatrix(const std::vector<BaseType>& base_vector,
                                 gsMatrix<TargetType>& result) {
        // Check for size
        const unsigned rows = result.rows();
        const unsigned cols = result.cols();
        // Size check
        if (base_vector.size() != (rows * cols)) {
            GISMO_ERROR(
                "Input array has the wrong size or could not be converted");
        }
        // Converting into gsMatrix (manipulating directly on gsMatrix is
        // more efficient if T=InputType)
        for (unsigned i = 0; i < rows; ++i) {
            for (unsigned j = 0; j < cols; ++j) {
                result(i, j) =
                    static_cast<TargetType>(base_vector[i * cols + j]);
            }
        }
    }

    /**
     * @brief Cast a vector of a base type into a vector of TargetType
     *
     * @tparam BaseType Type as encoded into the input file
     * @tparam TargetType Type as requested from target
     * @param base_vector decoded std::vector with base type
     * @param result decoded std::vector with target type
     */
    template <typename BaseType, typename TargetType>
    static void CopyIntoVector(const std::vector<BaseType>& base_vector,
                               std::vector<TargetType>& result) {
        std::transform(base_vector.cbegin(), base_vector.cend(),
                       std::back_inserter(result), [](const BaseType& c) {
                           return static_cast<TargetType>(c);
                       });
    }

   public:
    /**
     * @brief Helper routine for std::vector data
     *
     * @tparam BaseType type of individual data entries
     * @param data_vector data to be encoded
     * @return std::string encoded data
     */
    template <typename BaseType>
    static std::string Encode(const std::vector<BaseType>& data_vector) {
        const ByteRepresentation* vector_as_bytes_ptr =
            reinterpret_cast<const ByteRepresentation*>(&data_vector[0]);

        // Number of bytes for an entry
        constexpr const std::size_t length_of_entry{sizeof(BaseType{})};
        // Minimum number of bytes required
        const std::size_t minimum_n_bytes_required =
            length_of_entry * data_vector.size();

        return Encode_(vector_as_bytes_ptr, minimum_n_bytes_required);
    }

    /**
     * @brief Helper routine for gsMatrix Types (non-sparse)
     *
     * @tparam BaseType type of individual data entries
     * @param data_vector data to be encoded
     * @param row_wise Encode in row_wise style
     * @return std::string encoded data
     */
    template <typename BaseType>
    static std::string Encode(const gsMatrix<BaseType>& data_vector,
                              const bool& row_wise = true) {
        GISMO_ASSERT(std::is_arithmetic<BaseType>::value, // can be static
                      "Encoding is unsafe for non-arithmetic types.");
        // We need to ensure that the export is in the demanded export order, if
        // the data is only a vector (i.e. either col or row is 1) or if the
        // storage scheme is coherent with the demanded export order
        if ((data_vector.cols() == 1) || (data_vector.rows() == 1) ||
            (gsMatrix<BaseType>::IsRowMajor == row_wise)) {
            const ByteRepresentation* vector_as_bytes_ptr =
                reinterpret_cast<const ByteRepresentation*>(data_vector.data());

            // Number of bytes for an entry
            constexpr const std::size_t length_of_entry{sizeof(BaseType{})};
            // Minimum number of bytes required
            const std::size_t minimum_n_bytes_required =
                length_of_entry * data_vector.size();

            return Encode_(vector_as_bytes_ptr, minimum_n_bytes_required);
        } else {
            // Here we need to flip the order of elements, for simplicity the
            // data is copied into a temporary vector (works best with current
            // implementation but might be slower than other approaches)
            std::vector<BaseType> copy_of_matrix;
            copy_of_matrix.reserve(data_vector.rows() * data_vector.cols());

            // For readability we manipulate the Matrix using (memory safe)
            // operator() overloads
            if (row_wise) {
                for (index_t i = 0; i < data_vector.rows(); ++i) {
                    for (index_t j = 0; j < data_vector.cols(); ++j) {
                        copy_of_matrix.push_back(data_vector(i, j));
                    }
                }
            } else {
                for (index_t j = 0; j < data_vector.cols(); ++j) {
                    for (index_t i = 0; i < data_vector.rows(); ++i) {
                        copy_of_matrix.push_back(data_vector(i, j));
                    }
                }
            }
            // Use vector overload
            return Encode(copy_of_matrix);
        }
    }

    /**
     * @brief Reading a B64 string, transforming it into a vector of a
     * specific type
     *
     * @todo: In the future copies could be avoided by using string_view
     *
     * @tparam OutputType target type
     */
    template <typename OutputType>
    static std::vector<OutputType> Decode(const std::string& base64string) {
        // Safeguard
        const std::string& base64string_trimmed = trimWhitespaces(base64string);
        // Check validity of string
        GISMO_ASSERT(isValidBase64String(base64string_trimmed),
                     "Validity check failed");

        // Init return value
        const std::size_t number_of_groups{base64string_trimmed.size() / 4};
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
                buffer[i] = base64string_trimmed[buffer_index + i] != '='
                                ? char_decode_table(static_cast<unsigned>(
                                      base64string_trimmed[buffer_index + i]))
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

    /**
     * @brief Decode a string and copy into requested gismo Type
     *
     * Read a base64 string depending on a format flag (that determines base
     * type) and write it into gsMatrix<ScalarType>. Other overloads provide
     * functions to copy into a std::vector<ScalarType>
     *
     * @tparam ScalarType
     * @param base64_string
     * @param base_type_flag_
     * @param result
     */
    template <typename ScalarType>
    static void DecodeIntoGsType(const std::string& base64_string,
                                 const std::string& base_type_flag_,
                                 gsMatrix<ScalarType>& result) {
        // Format flag in this function is case sensitive
        GISMO_ASSERT(
            std::none_of(base_type_flag_.begin(), base_type_flag_.end(),
                         isupper),
            "Format flag {ascii, b64float64, ...} must be all lowercase.");

        // Perform type checks (no integral to floting point conversion)
        if (std::is_integral<ScalarType>::value ^
            (base_type_flag_.find("int") != std::string::npos)) {
            GISMO_ERROR(
                "Conversions from integral to floating type and vice-versa is "
                "not allowed!");
        }

        // Perform the actual input (using the proper encoding type)
        if (base_type_flag_ == "b64uint16") {  // Unsigned int 16
            CopyIntoGsMatrix(Decode<uint16_t>(base64_string), result);
        } else if (base_type_flag_ == "b64uint32") {  // Unsigned int 32
            CopyIntoGsMatrix(Decode<uint32_t>(base64_string), result);
        } else if (base_type_flag_ == "b64bint64") {  // Unsigned int 64
            CopyIntoGsMatrix(Decode<uint64_t>(base64_string), result);
        } else if (base_type_flag_ == "b64int16") {  // Int 16
            CopyIntoGsMatrix(Base64::Decode<int16_t>(base64_string), result);
        } else if (base_type_flag_ == "b64int32") {  // Int 32
            CopyIntoGsMatrix(Base64::Decode<int32_t>(base64_string), result);
        } else if (base_type_flag_ == "b64int64") {  // Int 64
            CopyIntoGsMatrix(Base64::Decode<int64_t>(base64_string), result);
        } else if (base_type_flag_ == "b64float32") {  // Float 32
            CopyIntoGsMatrix(Base64::Decode<float>(base64_string), result);
        } else if (base_type_flag_ == "b64float64") {  // Float 64
            CopyIntoGsMatrix(Base64::Decode<double>(base64_string), result);
        } else {
            GISMO_ERROR("Reading matrix from XML found unknown type");
        }
    }

    /**
     * @brief Decode a string and copy into requested gismo Type
     *
     * Read a base64 string depending on a format flag (that determines base
     * type) and write it into std::vector<ScalarType>. Other overloads provide
     * functions to copy into a gsMatrix<ScalarType>
     *
     * @tparam ScalarType
     * @param base64_string
     * @param base_type_flag_
     * @param result
     */
    template <typename ScalarType>
    static void DecodeIntoGsType(const std::string& base64_string,
                                 const std::string& base_type_flag_,
                                 std::vector<ScalarType>& result) {
        // Format flag in this function is case sensitive
        GISMO_ASSERT(
            std::none_of(base_type_flag_.begin(), base_type_flag_.end(),
                         isupper),
            "Format flag {ascii, b64float64, ...} must be all lowercase.");

        // Perform type checks (no integral to floting point conversion)
        if (std::is_integral<ScalarType>::value ^
            (base_type_flag_.find("int") != std::string::npos)) {
            GISMO_ERROR(
                "Conversions from integral to floating type and vice-versa is "
                "not allowed!");
        }

        // Perform the actual input (using the proper encoding type)
        if (base_type_flag_ == "b64uint16") {  // Unsigned int 16
            CopyIntoVector(Decode<uint16_t>(base64_string), result);
        } else if (base_type_flag_ == "b64uint32") {  // Unsigned int 32
            CopyIntoVector(Decode<uint32_t>(base64_string), result);
        } else if (base_type_flag_ == "b64bint64") {  // Unsigned int 64
            CopyIntoVector(Decode<uint64_t>(base64_string), result);
        } else if (base_type_flag_ == "b64int16") {  // Int 16
            CopyIntoVector(Base64::Decode<int16_t>(base64_string), result);
        } else if (base_type_flag_ == "b64int32") {  // Int 32
            CopyIntoVector(Base64::Decode<int32_t>(base64_string), result);
        } else if (base_type_flag_ == "b64int64") {  // Int 64
            CopyIntoVector(Base64::Decode<int64_t>(base64_string), result);
        } else if (base_type_flag_ == "b64float32") {  // Float 32
            CopyIntoVector(Base64::Decode<float>(base64_string), result);
        } else if (base_type_flag_ == "b64float64") {  // Float 64
            CopyIntoVector(Base64::Decode<double>(base64_string), result);
        } else {
            GISMO_ERROR("Reading matrix from XML found unknown type");
        }
    }
};
}  // namespace gismo
