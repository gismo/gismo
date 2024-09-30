/**
 * @file gsMpiTraits.h
 *
 * @brief Traits classes for mapping types onto MPI_Datatype.
 * @author C.Hofer, taken from DUNE - author Markus Blatt
 */

#pragma once

#include <utility>
#include <cstddef>

namespace gismo
{

  /**
   * @brief A traits class describing the mapping of types onto MPI_Datatypes.
   *
   * Specializations exist for the default types.
   * Specializations should provide a static method
   * \code
   * static MPI_Datatype getType();
   * \endcode
   */
  template<typename T>
  struct MPITraits
  {
  private:
    MPITraits(){}
    MPITraits(const MPITraits&){}
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;
  public:
    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(sizeof(T),MPI_BYTE,&datatype);
        MPI_Type_commit(&datatype);
      }
      return datatype;
    }

  };
  template<class T>
  MPI_Datatype MPITraits<T>::datatype = MPI_DATATYPE_NULL;

#ifndef DOXYGEN

  // A Macro for defining traits for the primitive data types
#define ComposeMPITraits(p,m) \
  template<> \
  struct MPITraits<p>{ \
    static inline MPI_Datatype getType(){ \
      return m; \
    } \
  }

  // Identify inbuild MPI types with the C types
  ComposeMPITraits(char, MPI_CHAR);
  ComposeMPITraits(unsigned char,MPI_UNSIGNED_CHAR);
  ComposeMPITraits(short,MPI_SHORT);
  ComposeMPITraits(unsigned short,MPI_UNSIGNED_SHORT);
  ComposeMPITraits(int,MPI_INT);
  ComposeMPITraits(unsigned int,MPI_UNSIGNED);
  ComposeMPITraits(long,MPI_LONG);
  ComposeMPITraits(unsigned long,MPI_UNSIGNED_LONG);
  ComposeMPITraits(float,MPI_FLOAT);
  ComposeMPITraits(double,MPI_DOUBLE);
  ComposeMPITraits(long double,MPI_LONG_DOUBLE);


#undef ComposeMPITraits

  /// Forward declaration of gsVector class
  template<class T, int _Rows, int _Options> class gsVector;

  /// Specialization for fixed-size gsVector class
  template<class T, int _Rows, int _Options>
  struct MPITraits<gsVector<T, _Rows, _Options> >
  {
  public:
    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL)
        {
          MPI_Type_contiguous(_Rows, MPITraits<T>::getType(), &vectortype);
          MPI_Type_commit(&vectortype);
          gsVector<T, _Rows, _Options> fvector;
          MPI_Aint base;
          MPI_Aint displ;
          MPI_Get_address(&fvector, &base);
          MPI_Get_address(&(fvector[0]), &displ);
          displ -= base;
          int length[1]={1};

          MPI_Type_create_struct(1, length, &displ, &vectortype, &datatype);
          MPI_Type_commit(&datatype);
        }
      return datatype;
    }

  private:
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;
  };

  template<class T, int _Rows, int _Options>
  MPI_Datatype MPITraits<gsVector<T, _Rows, _Options> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int _Rows, int _Options>
  MPI_Datatype MPITraits<gsVector<T, _Rows, _Options> >::vectortype = {MPI_DATATYPE_NULL};



  /// Forward declaration of gsMatrix class
  template<class T, int _Rows, int _Cols, int _Options> class gsMatrix;

  /// Specialization for fixed-size gsMatrix class
  template<class T, int _Rows, int _Cols, int _Options>
  struct MPITraits<gsMatrix<T, _Rows, _Cols, _Options> >
  {
  public:
    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL)
        {
          MPI_Type_contiguous(_Rows*_Cols, MPITraits<T>::getType(), &matrixtype);
          MPI_Type_commit(&matrixtype);
          gsMatrix<T, _Rows, _Cols, _Options> fmatrix;
          MPI_Aint base;
          MPI_Aint displ;
          MPI_Get_address(&fmatrix, &base);
          MPI_Get_address(&(fmatrix(0,0)), &displ);
          displ -= base;
          int length[1]={1};

          MPI_Type_create_struct(1, length, &displ, &matrixtype, &datatype);
          MPI_Type_commit(&datatype);
        }
      return datatype;
    }

  private:
    static MPI_Datatype datatype;
    static MPI_Datatype matrixtype;
  };

  template<class T, int _Rows, int _Cols, int _Options>
  MPI_Datatype MPITraits<gsMatrix<T, _Rows, _Cols, _Options> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int _Rows, int _Cols, int _Options>
  MPI_Datatype MPITraits<gsMatrix<T, _Rows, _Cols, _Options> >::matrixtype = {MPI_DATATYPE_NULL};


  /// Forward declaration gsSparseMatrix
  template<typename T, int _Options, typename _Index> class gsSparseMatrix;

  template<typename T, int _Options, typename _Index>
  struct MPITraits<gsSparseMatrix<T, _Options, _Index> >
  {
  public:
    static inline MPI_Datatype getType()
    {
      if (datatype==MPI_DATATYPE_NULL)
        {
          gsSparseMatrix<T, _Options, _Index> mat;
          MPI_Aint base;
          MPI_Get_address(mat.data, &base);
          //MPI_Aint* displ = new int[n];
          //for (int i=0; i<n; ++i)
          //  {
          //    MPI_Get_address(mat[i], &displ[i]);
          //    displ[i] -= base;
          //  }

          //MPI_Type_hindexed(n, m, displ, MPITraits<T>::getType(), &datatype);
          //MPI_Type_commit(&datatype);
        }
      return datatype;
    }

  private:
    static MPI_Datatype datatype;
    static MPI_Datatype matrixtype;
  };

  template<class T, int _Options, typename _Index>
  MPI_Datatype MPITraits<gsSparseMatrix<T, _Options, _Index> >::datatype = MPI_DATATYPE_NULL;
  template<class T, int _Options, typename _Index>
  MPI_Datatype MPITraits<gsSparseMatrix<T, _Options, _Index> >::matrixtype = {MPI_DATATYPE_NULL};

  /*
  template<int k>
  class bigunsignedint;

  template<int k>
  struct MPITraits<bigunsignedint<k> >
  {
    static MPI_Datatype datatype;
    static MPI_Datatype vectortype;

    static inline MPI_Datatype getType()
    {
      if(datatype==MPI_DATATYPE_NULL) {
        MPI_Type_contiguous(bigunsignedint<k>::n, MPITraits<std::uint16_t>::getType(),
                            &vectortype);
        //MPI_Type_commit(&vectortype);
        bigunsignedint<k> data;
        MPI_Aint base;
        MPI_Aint displ;
        MPI_Get_address(&data, &base);
        MPI_Get_address(&(data.digit), &displ);
        displ -= base;
        int length[1]={1};
        MPI_Type_create_struct(1, length, &displ, &vectortype, &datatype);
        MPI_Type_commit(&datatype);
      }
      return datatype;
    }
  };

  template<int k>
  MPI_Datatype MPITraits<bigunsignedint<k> >::datatype = MPI_DATATYPE_NULL;
  template<int k>
  MPI_Datatype MPITraits<bigunsignedint<k> >::vectortype = MPI_DATATYPE_NULL;
*/

  template<typename T1, typename T2>
  struct MPITraits<std::pair<T1,T2> >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2>
  MPI_Datatype MPITraits<std::pair<T1,T2> >::getType()
  {
    if(type==MPI_DATATYPE_NULL)
    {
        int length[2] = {1, 1};
        MPI_Aint disp[2];
        MPI_Datatype types[2] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType()};

        typedef std::pair<T1, T2> Pair;
        //static_assert(std::is_standard_layout<Pair>::value, "offsetof() is only defined for standard layout types");
        disp[0] = offsetof(Pair, first);
        disp[1] = offsetof(Pair, second);

        MPI_Datatype tmp;
        MPI_Type_create_struct(2, length, disp, types, &tmp);

        MPI_Type_create_resized(tmp, 0, sizeof(Pair), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2>
  MPI_Datatype MPITraits<std::pair<T1,T2> >::type=MPI_DATATYPE_NULL;


  template<typename T1, typename T2, typename T3>
  struct MPITraits<std::tuple<T1,T2,T3> >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2, typename T3>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3> >::getType()
  {
    if(type==MPI_DATATYPE_NULL)
    {
        int length[3] = {1, 1, 1};
        MPI_Aint disp[3];
        MPI_Datatype types[3] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType(),
                                 MPITraits<T3>::getType()};

        typedef std::tuple<T1, T2, T3> Tuple3;
        Tuple3 dummy_tuple;
        MPI_Aint base_address;

        MPI_Get_address(&dummy_tuple, &base_address);
        MPI_Get_address(&(std::get<0>(dummy_tuple)), &disp[0]);
        MPI_Get_address(&(std::get<1>(dummy_tuple)), &disp[1]);
        MPI_Get_address(&(std::get<2>(dummy_tuple)), &disp[2]);
        disp[0] = MPI_Aint_diff(disp[0], base_address);
        disp[1] = MPI_Aint_diff(disp[1], base_address);
        disp[2] = MPI_Aint_diff(disp[2], base_address);


        MPI_Datatype tmp;
        MPI_Type_create_struct(3, length, disp, types, &tmp);

        MPI_Type_create_resized(tmp, 0, sizeof(Tuple3), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2, typename T3>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3> >::type=MPI_DATATYPE_NULL;

  template<typename T1, typename T2, typename T3, typename T4>
  struct MPITraits<std::tuple<T1,T2,T3,T4> >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2, typename T3, typename T4>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4> >::getType()
  {
    if(type==MPI_DATATYPE_NULL)
    {
        int length[4] = {1, 1, 1, 1};
        MPI_Aint disp[4];
        MPI_Datatype types[4] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType(),
                                 MPITraits<T3>::getType(),
                                 MPITraits<T4>::getType()};

        typedef std::tuple<T1, T2, T3, T4> Tuple4;
        Tuple4 dummy_tuple;
        MPI_Aint base_address;

        MPI_Get_address(&dummy_tuple, &base_address);
        MPI_Get_address(&(std::get<0>(dummy_tuple)), &disp[0]);
        MPI_Get_address(&(std::get<1>(dummy_tuple)), &disp[1]);
        MPI_Get_address(&(std::get<2>(dummy_tuple)), &disp[2]);
        MPI_Get_address(&(std::get<3>(dummy_tuple)), &disp[3]);
        disp[0] = MPI_Aint_diff(disp[0], base_address);
        disp[1] = MPI_Aint_diff(disp[1], base_address);
        disp[2] = MPI_Aint_diff(disp[2], base_address);
        disp[3] = MPI_Aint_diff(disp[3], base_address);

        MPI_Datatype tmp;
        MPI_Type_create_struct(4, length, disp, types, &tmp);

        MPI_Type_create_resized(tmp, 0, sizeof(Tuple4), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2, typename T3, typename T4>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4> >::type=MPI_DATATYPE_NULL;

  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  struct MPITraits<std::tuple<T1,T2,T3,T4,T5> >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4,T5> >::getType()
  {
    if(type==MPI_DATATYPE_NULL)
    {
        int length[5] = {1, 1, 1, 1, 1};
        MPI_Aint disp[5];
        MPI_Datatype types[5] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType(),
                                 MPITraits<T3>::getType(),
                                 MPITraits<T4>::getType(),
                                 MPITraits<T5>::getType()};

        typedef std::tuple<T1, T2, T3, T4, T5> Tuple5;
        Tuple5 dummy_tuple;
        MPI_Aint base_address;

        MPI_Get_address(&dummy_tuple, &base_address);
        MPI_Get_address(&(std::get<0>(dummy_tuple)), &disp[0]);
        MPI_Get_address(&(std::get<1>(dummy_tuple)), &disp[1]);
        MPI_Get_address(&(std::get<2>(dummy_tuple)), &disp[2]);
        MPI_Get_address(&(std::get<3>(dummy_tuple)), &disp[3]);
        MPI_Get_address(&(std::get<4>(dummy_tuple)), &disp[4]);
        disp[0] = MPI_Aint_diff(disp[0], base_address);
        disp[1] = MPI_Aint_diff(disp[1], base_address);
        disp[2] = MPI_Aint_diff(disp[2], base_address);
        disp[3] = MPI_Aint_diff(disp[3], base_address);
        disp[4] = MPI_Aint_diff(disp[4], base_address);

        MPI_Datatype tmp;
        MPI_Type_create_struct(5, length, disp, types, &tmp);

        MPI_Type_create_resized(tmp, 0, sizeof(Tuple5), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4,T5> >::type=MPI_DATATYPE_NULL;

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  struct MPITraits<std::tuple<T1,T2,T3,T4,T5,T6> >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4,T5,T6> >::getType()
  {
    if(type==MPI_DATATYPE_NULL)
    {
        int length[6] = {1, 1, 1, 1, 1, 1};
        MPI_Aint disp[6];
        MPI_Datatype types[6] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType(),
                                 MPITraits<T3>::getType(),
                                 MPITraits<T4>::getType(),
                                 MPITraits<T5>::getType(),
                                 MPITraits<T6>::getType()};

        typedef std::tuple<T1, T2, T3, T4, T5, T6> Tuple6;
        Tuple6 dummy_tuple;
        MPI_Aint base_address;

        MPI_Get_address(&dummy_tuple, &base_address);
        MPI_Get_address(&(std::get<0>(dummy_tuple)), &disp[0]);
        MPI_Get_address(&(std::get<1>(dummy_tuple)), &disp[1]);
        MPI_Get_address(&(std::get<2>(dummy_tuple)), &disp[2]);
        MPI_Get_address(&(std::get<3>(dummy_tuple)), &disp[3]);
        MPI_Get_address(&(std::get<4>(dummy_tuple)), &disp[4]);
        MPI_Get_address(&(std::get<5>(dummy_tuple)), &disp[5]);
        disp[0] = MPI_Aint_diff(disp[0], base_address);
        disp[1] = MPI_Aint_diff(disp[1], base_address);
        disp[2] = MPI_Aint_diff(disp[2], base_address);
        disp[3] = MPI_Aint_diff(disp[3], base_address);
        disp[4] = MPI_Aint_diff(disp[4], base_address);
        disp[5] = MPI_Aint_diff(disp[5], base_address);

        MPI_Datatype tmp;
        MPI_Type_create_struct(6, length, disp, types, &tmp);

        MPI_Type_create_resized(tmp, 0, sizeof(Tuple6), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4,T5,T6> >::type=MPI_DATATYPE_NULL;

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  struct MPITraits<std::tuple<T1,T2,T3,T4,T5,T6,T7> >
  {
  public:
    inline static MPI_Datatype getType();
  private:
    static MPI_Datatype type;
  };
  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4,T5,T6,T7> >::getType()
  {
    if(type==MPI_DATATYPE_NULL)
    {
        int length[7] = {1, 1, 1, 1, 1, 1, 1};
        MPI_Aint disp[7];
        MPI_Datatype types[7] = {MPITraits<T1>::getType(),
                                 MPITraits<T2>::getType(),
                                 MPITraits<T3>::getType(),
                                 MPITraits<T4>::getType(),
                                 MPITraits<T5>::getType(),
                                 MPITraits<T6>::getType(),
                                 MPITraits<T7>::getType()};

        typedef std::tuple<T1, T2, T3, T4, T5, T6, T7> Tuple7;
        Tuple7 dummy_tuple;
        MPI_Aint base_address;

        MPI_Get_address(&dummy_tuple, &base_address);
        MPI_Get_address(&(std::get<0>(dummy_tuple)), &disp[0]);
        MPI_Get_address(&(std::get<1>(dummy_tuple)), &disp[1]);
        MPI_Get_address(&(std::get<2>(dummy_tuple)), &disp[2]);
        MPI_Get_address(&(std::get<3>(dummy_tuple)), &disp[3]);
        MPI_Get_address(&(std::get<4>(dummy_tuple)), &disp[4]);
        MPI_Get_address(&(std::get<5>(dummy_tuple)), &disp[5]);
        MPI_Get_address(&(std::get<6>(dummy_tuple)), &disp[6]);
        disp[0] = MPI_Aint_diff(disp[0], base_address);
        disp[1] = MPI_Aint_diff(disp[1], base_address);
        disp[2] = MPI_Aint_diff(disp[2], base_address);
        disp[3] = MPI_Aint_diff(disp[3], base_address);
        disp[4] = MPI_Aint_diff(disp[4], base_address);
        disp[5] = MPI_Aint_diff(disp[5], base_address);
        disp[6] = MPI_Aint_diff(disp[6], base_address);

        MPI_Datatype tmp;
        MPI_Type_create_struct(7, length, disp, types, &tmp);

        MPI_Type_create_resized(tmp, 0, sizeof(Tuple7), &type);
        MPI_Type_commit(&type);
        MPI_Type_free(&tmp);
    }
    return type;
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  MPI_Datatype MPITraits<std::tuple<T1,T2,T3,T4,T5,T6,T7> >::type=MPI_DATATYPE_NULL;
}

#endif
