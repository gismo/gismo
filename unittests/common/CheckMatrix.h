/** @file common/CheckMatrix.h

    @brief common UnitTest++ macro definitions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H. Weiner
**/

#ifndef COMMONS_CHECKMATRIX_H
#define COMMONS_CHECKMATRIX_H

#include <gismo.h>

namespace UnitTest {
  
/**
 * check that a matrix 'actual' is equal to an 'expected' matrix
 * (or at least close, with respect to the 'tolerance' value)
 */
template< typename T >
void CheckMatrixClose(TestResults& results, gismo::gsMatrix<T> const& expected, gismo::gsMatrix<T> const& actual,
                      T const& tolerance, TestDetails const& details)
{
 bool equal = true;
  const index_t rows = actual.rows();
  const index_t cols = actual.cols();
  
  equal &= (expected.rows() == rows);
  equal &= (expected.cols() == cols);
  for (index_t i=0; equal && i<rows; ++i)
  {
      for (index_t j=0; j<cols; ++j)
      {
          equal &= (gismo::math::abs( expected(i, j) - actual(i, j) )<= tolerance);
      }
  }
  
  if (!equal)
  {
      UnitTest::MemoryOutStream stream;
      
      if (expected.rows() != rows) {
          stream << "Expected rows=" << expected.rows() << " but was " << rows <<"!";
          results.OnTestFailure(details, stream.GetText());
          return;
      }
      if (expected.cols() != cols) {
          stream << "Expected cols=" << expected.cols() << " but was " << cols <<"!";
          results.OnTestFailure(details, stream.GetText());
          return;
      }

      stream << "Expected [ ";    

      for (int expectedRow = 0; expectedRow < rows; ++expectedRow)
      {
          stream << "[ ";
          for (int expectedColumn = 0; expectedColumn < cols; ++expectedColumn)
              stream << expected(expectedRow, expectedColumn) << " ";
          stream << "] ";
      }

      stream << "] +/- " << tolerance << " but was [ ";

      for (int actualRow = 0; actualRow < rows; ++actualRow)
      {
          stream << "[ ";
          for (int actualColumn = 0; actualColumn < cols; ++actualColumn)
              stream << actual(actualRow,actualColumn) << " ";
          stream << "] ";
      }

      stream << "]";

      results.OnTestFailure(details, stream.GetText());
  }
}


template< typename T >
void CheckMatrixPartitionOfUnitClose(TestResults& results, gismo::gsMatrix<T> const& actual,
                      T const& tolerance, TestDetails const& details)
{
  bool result = true;
  const index_t cols = actual.cols();
  
  for (index_t j=0; j<cols; ++j)
  {
       const real_t sum = actual.col(j).sum();
       result &= (gismo::math::abs( sum - 1.0) <= tolerance);
  }
  
  
  if (!result)
  {
      UnitTest::MemoryOutStream stream;
      
      stream << "expected for all columns that sum='" << 1.0 << "' +/- '" << tolerance<< "', but was [";
      
      for (index_t j=0; j<cols; ++j)
      {
           const real_t sum = actual.col(j).sum();
           stream << "[expected='" << 1.0 << ", actual.cols('";
           stream << j << "').sum()='" << sum << "'']";
      }

      stream << "]";

      results.OnTestFailure(details, stream.GetText());
  }
}


}

#endif // COMMONS_CHECKMATRIX_H

