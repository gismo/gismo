/** @file gsKnotVectors_test.cpp

    @brief test knot vectors

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
**/

#include "gismo_unittest.h"

SUITE(gsKnotVectors_test)
{
    struct data
    {
        struct TestVector
        {
            std::vector<real_t> vector;
            std::string         name;
            int                 deg;

            TestVector (const std::vector<real_t> &v, const std::string &n, int d)
                : vector(v), name(n), deg(d)
            {}
        };
        std::vector<TestVector> testVector;


        data()
        {
            {
                std::vector<real_t> temp;
                std::string         name;

                name="deg0";
                temp.resize(3);
                temp[0] = 0;
                temp[1] = 0.5;
                temp[2] = 1;
                testVector.push_back(TestVector(temp,name,0));
            }
            {
                std::vector<real_t> temp;
                std::string         name;

                name="deg1";
                temp.resize(5);
                temp[0] = 0;
                temp[1] = 0;
                temp[2] = 0.5;
                temp[3] = 1;
                temp[4] = 1;
                testVector.push_back(TestVector(temp,name,0));
            }
        }

        // It would be great to merge the two tests to one but I am off to look for the problem in gsHierarchical_box_test.cpp, as Angelos promised me 0.50 EUR for solving it.
    };


    void testConstructor (const data::TestVector &test)
    {
        gsKnotVector<real_t> unique(test.vector, test.deg);
        gsKnotVector<real_t> copied = gsKnotVector<real_t>( unique );
        gsKnotVector<real_t> assigned = copied;
        gsKnotVector<real_t> iterated( test.deg, test.vector.begin(), test.vector.end() );

        const size_t size=test.vector.size();

        CHECK( static_cast<size_t>(unique.size())   == size);
        CHECK( static_cast<size_t>(copied.size())   == size);
        CHECK( static_cast<size_t>(assigned.size()) == size);
        CHECK( static_cast<size_t>(iterated.size()) == size);

        for( size_t i = 0; i < size; ++i)
        {
            CHECK( unique[i]   == test.vector[i] );
            CHECK( copied[i]   == test.vector[i] );
            CHECK( assigned[i] == test.vector[i] );
            CHECK( iterated[i] == test.vector[i] );
        }
    }
    TEST_FIXTURE( data, constructors)
    {
        for (size_t t=0; t<testVector.size();++t)
            testConstructor(testVector[t]);
    }


    TEST_FIXTURE( data, constructorWithRegularKnots )
    {
        std::vector<real_t> uniqueKnots = testVector[0].vector;

        gsKnotVector<real_t> regular(uniqueKnots, 2, 1);
        CHECK( regular[0] == 0 );
        CHECK( regular[1] == 0 );
        CHECK( regular[2] == 0 );
        CHECK( regular[3] == 0.5 );
        CHECK( regular[4] == 1 );
        CHECK( regular[5] == 1 );
        CHECK( regular[6] == 1 );

        CHECK( regular.size() == 7 );
    }

    TEST(notClampedKnotVector)
    {
	gsKnotVector<real_t> KV(0.0, 1.0, 1, 1, 1, 1);

	std::vector<real_t> repKnots = KV.get();
	CHECK( repKnots.size() == 5 );
	CHECK_CLOSE( repKnots[0], -0.5, 1e-6 );
	CHECK_CLOSE( repKnots[1],  0.0, 1e-6 );
	CHECK_CLOSE( repKnots[2],  0.5, 1e-6 );
	CHECK_CLOSE( repKnots[3],  1.0, 1e-6 );
	CHECK_CLOSE( repKnots[4],  1.5, 1e-6 );

	const gsKnotVector<real_t>::mult_t* multSum = KV.multSumData();
	CHECK( *multSum == 1 );
	CHECK( *(multSum + 1) == 2 );
	CHECK( *(multSum + 2) == 3 );
	CHECK( *(multSum + 3) == 4 );
	CHECK( *(multSum + 4) == 5 );
    }

    TEST(gsBasis_refine)
    {
        /* This test checks that the refinement specified in the documentation of
         * gsBasis::refine() really refines the way it is described there
         * (or it has been today, 12.10.2015). */
        gsKnotVector<real_t> kv(0, 1, 4, 4, 1, 3);
        gsTensorBSplineBasis<2,real_t> tensorBasis(kv,kv);
        gsTHBSplineBasis<2> thbBasis(tensorBasis);
        gsMatrix<real_t> refinementBoxes(2,4);
        refinementBoxes << 0, 0.2, 0.8, 1, 0.4, 0.6, 0.2, 0.4;
        thbBasis.refine(refinementBoxes);

        gsMatrix<index_t> corrLowLeft(7,2), corrUppRigh(7,2);
        gsMatrix<index_t> compLowLeft, compUppRigh;
        gsVector<index_t> corrLvl(7), compLvl;
        thbBasis.tree().getBoxes(compLowLeft, compUppRigh, compLvl);

        // These matrices both contain indices of the second level (draw a picture;)).
        corrLowLeft << 8, 4,
                8, 2,
                8, 0,
                2, 0,
                0, 6,
                0, 4,
                0, 0;
        corrUppRigh << 10, 10,
                10,  4,
                10,  2,
                8, 10,
                2, 10,
                2,  6,
                2,  4;
        // Levels of the boxes.
        corrLvl << 0,
                1,
                0,
                0,
                0,
                1,
                0;
        // Check the sizes first.
        CHECK( compLowLeft.rows()  == corrLowLeft.rows()  &&
               compLowLeft.cols()  == corrLowLeft.cols()  &&
               compUppRigh.rows() == corrUppRigh.rows() &&
               compUppRigh.cols() == corrUppRigh.cols());
        CHECK( compLvl.size()      == corrLvl.size() );

        // Then check also the contents.
        CHECK( (compLvl.array()    == corrLvl.array()).all() );
        CHECK( (corrLowLeft.array() == compLowLeft.array()).all() );
        CHECK( (corrUppRigh.array() == compUppRigh.array()).all() );
    }

    template <typename KV>
    void testFindSpan (const KV &knots)
    {
        size_t deg=knots.degree();
        size_t span;
        typename KV::const_iterator spanIt;
        for (size_t k=0; k<knots.size()-deg-1;++k)
        {
            if(knots[k+1]!=knots[k])
            {
                // check for knot aligned points
                span   = knots.iFind(knots[k]) - knots.begin();
                spanIt = knots.iFind(knots[k]);
                CHECK(span==k);
                CHECK(static_cast<size_t>(spanIt-knots.begin())==static_cast<size_t>(span));

                // check for a point in the middle of the element
                span   = knots.iFind(knots[k]/2+knots[k+1]/2) - knots.begin();
                spanIt = knots.iFind(knots[k]/2+knots[k+1]/2);
                CHECK(span==k);
                CHECK(static_cast<size_t>(spanIt-knots.begin())==static_cast<size_t>(span));
            }
        }
        // check last point in the domain
        const size_t lastKnot=knots.size()-deg-1;
        span   = knots.iFind(knots[lastKnot]) - knots.begin();
        spanIt = knots.iFind(knots[lastKnot]);
        CHECK(span==lastKnot-1);
        CHECK(static_cast<size_t>(spanIt-knots.begin())==static_cast<size_t>(span));
    }

    TEST(findSpan_KV)
    {
        int deg=2;
        real_t data[]={1,1,1,1.1,1.3,1.7,1.8,1.9,2,2,2};
        std::vector<real_t> v(data,data+sizeof(data)/sizeof(real_t));
        gsKnotVector<real_t> kv1(deg,v.begin(),v.end());
        testFindSpan(kv1);
        gsKnotVector<real_t> kv2(deg,v.begin(),v.end());
        testFindSpan(kv2);
    }
}

SUITE(gsKnotVectors_test_2)
{

    typedef gsKnotVector<real_t>::uiterator uniqIter;
    typedef gsKnotVector<real_t>::iterator  knotIter;
    typedef gsKnotVector<real_t>::mult_t   mult_t;

    /// Checks whether KV has: the correct knots (corrKnots), the
    /// correct knots without repetitions (corrUKnots), the correct
    /// number of knots/ending positions (corrEndPos) and the correct
    /// degree (corrDeg)
    void compareUKVs( const gsKnotVector<real_t> & KV,
                      real_t corrKnots[],
                      real_t corrUKnots[],
                      mult_t corrEndPos[],
                      int    corrDeg)
    {
        knotIter kit( KV.begin() );
        knotIter kitEnd( KV.end() );
        mult_t i = 0;
        for(  ; kit != kitEnd; ++kit, ++i )
            CHECK_CLOSE( *kit, corrKnots[i], EPSILON );
            //CHECK( math::abs(*kit - corrKnots[i]) <= std::numeric_limits<real_t>::epsilon()*100 );

        uniqIter uit = KV.ubegin();
        uniqIter uitEnd = KV.uend();
        i = 0;
        for( ; uit != uitEnd; ++uit, ++i )
        {
            CHECK_CLOSE( *uit, corrUKnots[i], EPSILON );
            //CHECK( math::abs(*uit - corrUKnots[i])<= std::numeric_limits<real_t>::epsilon()*100 );
            CHECK( uit.multSum() == corrEndPos[i] );
        }
        CHECK(KV.degree() == corrDeg);
    }

    TEST( uFind )
    {
        {
            gsKnotVector<real_t> KV( 0, 1, 3, 3, 1 );

            CHECK( KV.uFind(0.00).uIndex() == 0 );
            CHECK( KV.uFind(0.49).uIndex() == 1 );
            CHECK( KV.uFind(0.50).uIndex() == 2 );
            CHECK( KV.uFind(0.55).uIndex() == 2 );
            CHECK( KV.uFind(1.00).uIndex() == 3 );

            CHECK( KV.uFind(0.00).multSum() == 3 );
            CHECK( KV.uFind(0.49).multSum() == 4 );
            CHECK( KV.uFind(0.50).multSum() == 5 );
            CHECK( KV.uFind(0.55).multSum() == 5 );
            CHECK( KV.uFind(1.00).multSum() == 6 );
        }
        {
            real_t data[]={-0.1, 0.0, 0.0, 0.2 ,0.3, 0.3, 0.5, 0.7, 1.0, 1.1, 1.2};
            gsKnotVector<real_t> KV( 2, data, data+11);

            CHECK( KV.numLeftGhosts() == 1 );
            CHECK( KV[0] == KV(-1) ); // unique indices of left ghosts are negative
            CHECK( KV[1] == KV(0)  );
            CHECK( KV[10] == KV(7) );

            CHECK( KV.uFind(0.00).uIndex() == 0 );
            CHECK( KV.uFind(0.49).uIndex() == 2 );
            CHECK( KV.uFind(0.50).uIndex() == 3 );
            CHECK( KV.uFind(0.55).uIndex() == 3 );
            CHECK( KV.uFind(0.99).uIndex() == 4 );
            CHECK( KV.uFind(1.00).uIndex() == 4 );

            CHECK( KV.uFind(0.00).uCardinalIndex() == 1 );
            CHECK( KV.uFind(0.49).uCardinalIndex() == 3 );
            CHECK( KV.uFind(0.50).uCardinalIndex() == 4 );
            CHECK( KV.uFind(0.55).uCardinalIndex() == 4 );
            CHECK( KV.uFind(0.99).uCardinalIndex() == 5 );
            CHECK( KV.uFind(1.00).uCardinalIndex() == 5 );

            CHECK( KV.iFind(0.00) - KV.begin() == 2 );
            CHECK( KV.iFind(0.49) - KV.begin() == 5 );
            CHECK( KV.iFind(0.50) - KV.begin() == 6 );
            CHECK( KV.iFind(0.55) - KV.begin() == 6 );
            CHECK( KV.iFind(0.99) - KV.begin() == 7 );
            CHECK( KV.iFind(1.00) - KV.begin() == 7 );
        }

    }

    TEST( uniqIter )
    {
        // Note that KV has different multiplicities than in TEST(uFind).
        gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
        uniqIter uit = KV.ubegin();
        uniqIter uitEnd = KV.uend();

        real_t correctUniqueKnots[] = {0, .25, .5, .75, 1};
        mult_t correctIndices[] = {3, 5, 7, 9, 12};
        mult_t correctUIndices[] = {0, 1, 2, 3, 4};

        CHECK( uitEnd - uit == 5 );
        for( mult_t i = 0; uit != uitEnd; ++uit, ++i )
        {
            CHECK( correctUniqueKnots[i] == *uit );
            CHECK( correctIndices[i] == uit.multSum() );
            CHECK( correctUIndices[i] == uit.uIndex() );
        }
    }
    TEST( uniformRefine )
    {
        {
            gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
            KV.uniformRefine();
            real_t corrKnots[] = {0, 0, 0, .125, .25, .25, .375, .5, .5, .625, .75, .75, .875, 1, 1, 1};
            real_t corrUKnots[] = {0, .125, .25, .375, .5, .625, .75, .875, 1};
            mult_t corrEndPos[] = {3, 4, 6, 7, 9, 10, 12, 13, 16};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
            KV.uniformRefine(1,2);
            real_t corrKnots[] = {0, 0, 0, .125, .125, .25, .25, .375, .375, .5, .5, .625, .625, .75, .75, .875, .875, 1, 1, 1};
            real_t corrUKnots[] = {0, .125, .25, .375, .5, .625, .75, .875, 1};
            mult_t corrEndPos[] = {3, 5, 7, 9, 11, 13, 15, 17, 20};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
            KV.uniformRefine(3,1);
            real_t corrKnots[] = {0, 0, 0, 0.0625, .125, .1875, .25, .25, .3125, .375, .4375, .5, .5, .5625, .625, .6875, .75, .75, .8125, .875, .9375, 1, 1, 1};
            real_t corrUKnots[] = {0, 0.0625, .125, .1875, .25, .3125, .375, .4375, .5, .5625, .625, .6875, .75, .8125, .875, .9375, 1};
            mult_t corrEndPos[] = {3, 4, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 18, 19, 20, 21, 24};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
            KV.uniformRefine(0,1);
            real_t corrKnots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {3, 5, 7, 9, 12};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
            KV.uniformRefine(2,0);
            real_t corrKnots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {3, 5, 7, 9, 12};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
    }

    TEST( iFind )
    {
        // TODO: Check really all the cases, a bug has come through!
        gsKnotVector<real_t> KV( 0, 1, 3, 3, 2 );
        CHECK( KV.iFind(0) - KV.begin() == 2 );
        CHECK( KV.iFind(0.25) - KV.begin() == 4 );
        CHECK( KV.iFind(0.3) - KV.begin()  == 4 );
        CHECK( KV.iFind(1) - KV.begin()    == 8 );
    }

    TEST( constructor )
    {
        real_t uk[] = {0, .2, .5, .7, 2};
        std::vector<real_t> knots(uk, uk + sizeof(uk)/sizeof(real_t));
        {
            gsKnotVector<real_t> KV( knots, 2, 1 );
            real_t corrKnots[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
            real_t corrUKnots[] = {0, .2, .5, .7, 2};
            mult_t corrEndPos[] = {3, 4, 5, 6, 9};
            compareUKVs(KV, corrKnots, corrUKnots, corrEndPos, 2);
        }
        {
            gsKnotVector<real_t> KV( knots, 3, 0 );
            real_t corrKnots[] = {0, 0, 0, 0, .2, .2, .2, .5, .5, .5, .7, .7, .7, 2, 2, 2, 2};
            real_t corrUKnots[] = {0, .2, .5, .7, 2};
            mult_t corrEndPos[] = {4, 7, 10, 13, 17};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 3 );
        }

        real_t ak[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
        std::vector<real_t> aknots(ak, ak + sizeof(ak)/sizeof(real_t));
        {
            gsKnotVector<real_t> KV(2, aknots.begin(), aknots.end());
            real_t corrKnots[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
            real_t corrUKnots[] = {0, .2, .5, .7, 2};
            mult_t corrEndPos[] = {3, 4, 5, 6, 9};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            gsKnotVector<real_t> KV( give(aknots) );
            real_t corrKnots[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
            real_t corrUKnots[] = {0, .2, .5, .7, 2};
            mult_t corrEndPos[] = {3, 4, 5, 6, 9};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
            real_t corrKnots[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
            real_t corrUKnots[] = {0, .2, .5, .7, 2};
            mult_t corrEndPos[] = {3, 4, 5, 6, 9};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
    }

    TEST( XML_IO )
    {
        std::string filename("unittest_knotvector.xml");
        // The ugly number is ugly intentionally in order to test the conversion real_t <-> ASCII.
        real_t ak[] = {0, 0, 0, .2, .5, .7, .7123456789123456456789, 2, 2, 2};
        gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
        gsKnotVector<real_t> ReadKV;
        gsFileData<> write;
        write << KV;
        write.dump(filename);
        gsReadFile<>(filename, ReadKV);

        CHECK_ARRAY_CLOSE(KV.begin(), ReadKV.begin(), KV.size(), EPSILON);
        //CHECK( ReadKV == KV );
    }

    TEST( degreeElevate )
    {
        {
            gsKnotVector<real_t> KV(0,1,3,3,2);
            KV.degreeElevate();
            real_t corrKnots[] = {0, 0, 0, 0, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 1, 1, 1, 1};
            real_t corrUKnots[] = {0, 0.25, 0.5, 0.75, 1};
            mult_t corrEndPos[] = {4, 7, 10, 13, 17};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 3 );
        }
        {
            gsKnotVector<real_t> KV(0,1,3,3,2);
            KV.degreeElevate(3);
            real_t corrKnots[] = {0, 0, 0, 0, 0, 0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75, 0.75, 0.75, 1, 1, 1, 1, 1, 1};
            real_t corrUKnots[] = {0, 0.25, 0.5, 0.75, 1};
            mult_t corrEndPos[] = {6, 11, 16, 21, 27};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 5 );
        }
    }

    TEST( addConstant )
    {
        real_t ak[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
        // The degree is intentionally 3.
        gsKnotVector<real_t> KV(3, ak, ak + sizeof(ak)/sizeof(real_t));
        KV.addConstant(1.3);
        real_t corrKnots[] = {1.3, 1.3, 1.3, 1.5, 1.8, 2, 3.3, 3.3, 3.3};
        real_t corrUKnots[] = {1.3, 1.5, 1.8, 2, 3.3};
        mult_t corrEndPos[] = { 3, 4, 5, 6, 9 };
        compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 3 );
    }
    TEST( findElements )
    {
        real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
        // The degree is intentionally -1.
        gsKnotVector<real_t> KV(-1, knots, knots + sizeof(knots)/sizeof(real_t));  //(0,1,3,3,2);

        std::vector<index_t> spans;
        spans.push_back(0);
        KV.refineSpans(spans, 1);
        {
            real_t corrKnots[] = {0, 0, 0, 0.125, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1};
            real_t corrUKnots[] = {0, 0.125, 0.25, 0.5, 0.75, 1};
            mult_t corrEndPos[] = {3, 4, 6, 8, 10, 13};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, -1 );
        }
        spans[0] = 3;
        spans.push_back(4);
        KV.refineSpans(spans, 3);
        {
            real_t corrKnots[] = {0, 0, 0, 0.125, 0.25, 0.25, 0.5, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.75, 0.8125, 0.875, 0.9375, 1, 1, 1};
            real_t corrUKnots[] = {0, 0.125, 0.25, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1};
            mult_t corrEndPos[] = {3, 4, 6, 8, 9, 10, 11, 13, 14, 15, 16, 19};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, -1 );
        }
    }

    TEST( knotIter )
    {
        real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
        knotIter kit( KV.begin() );
        knotIter kitEnd( KV.end() );

        real_t correctKnots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};

        CHECK( kitEnd - kit == 12 );
        for( mult_t i = 0; kit != kitEnd; ++kit, ++i )
        {
            CHECK( correctKnots[i] == *kit );
            CHECK( i == kit - KV.begin() );
        }
    }

    void insertTest( real_t value,
                     int mult,
                     real_t corrKnots[],
                     real_t corrUKnots[],
                     mult_t corrEndPos[] )
    {
        real_t knots[] = {0, 0, 0, .25, .5, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 1 );
        KV.insert( value, mult );
        compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
    }

    TEST( insert )
    {
        {
            real_t corrKnots[] = {0, 0, 0, 0, .25, .5, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {4, 5, 6, 7, 10};
            insertTest( 0, 1, corrKnots, corrUKnots, corrEndPos );
        }
        {
            real_t corrKnots[] = {0, 0, 0, .25, .25, .25, .5, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {3, 6, 7, 8, 11};
            insertTest( .25, 2, corrKnots, corrUKnots, corrEndPos );
        }
        {
            real_t corrKnots[] = {0, 0, 0, .25, .3, .3, .3, .3, .3, .3, .3, .3, .3, .5, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .3, .5, .75, 1};
            mult_t corrEndPos[] = {3, 4, 13, 14, 15, 18};
            insertTest( .3, 9, corrKnots, corrUKnots, corrEndPos );
        }
        {
            real_t corrKnots[] = {0, 0, 0, .25, .5, .75, 1, 1, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {3, 4, 5, 6, 11};
            insertTest( 1, 2, corrKnots, corrUKnots, corrEndPos );
        }
        {
            real_t corrKnots[] = { -.27, 0, 0, 0, .25, .5, .75, 1, 1, 1};
            real_t corrUKnots[] = {-.27, 0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {1, 4, 5, 6, 7, 10};
            insertTest( -.27, 1, corrKnots, corrUKnots, corrEndPos );
        }
        {
            real_t corrKnots[] = {0, 0, 0, .25, .5, .75, 1, 1, 1, 3, 3};
            real_t corrUKnots[] = {0, .25, .5, .75, 1, 3};
            mult_t corrEndPos[] = {3, 4, 5, 6, 9, 11};
            insertTest( 3, 2, corrKnots, corrUKnots, corrEndPos );
        }
        {
            gsKnotVector<real_t> KV;
            KV.insert( .5 );
            KV.insert( .5 );
            real_t corrKnots[] = {.5, .5};
            real_t corrUKnots[] = {.5};
            mult_t corrEndPos[] = {2};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, -1 );
        }
    }

    TEST( batch_insert )
    {
        {
            real_t knots[] = {0, 0, 0, .25, .5, .75, 1, 1, 1};
            gsKnotVector<real_t> KV1(2, knots, knots+sizeof(knots)/sizeof(real_t)); //( 0, 1, 3, 3, 1 );
            gsKnotVector<real_t> KV2(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 1 );

            real_t newKnots[] = {0, 0.25, .3, .3, 1, 1};
            mult_t nNKnots = sizeof(newKnots) / sizeof(real_t);
            KV1.insert( newKnots, newKnots + nNKnots );

            for( mult_t i = 0; i < nNKnots; ++i )
                KV2.insert( newKnots[i] );

            CHECK( KV1 == KV2 );
        }
        {
            real_t knots[] = {0, 0, 0, .25, .5, .75, 1, 1, 1};
            mult_t nKnots = sizeof(knots) / sizeof(real_t);

            gsKnotVector<real_t> KV1( 2 ); // We have to set the degree because of the comparison later.
            gsKnotVector<real_t> KV2( 2, knots, knots + nKnots );  //( 0, 1, 3, 3, 1 );

            KV1.insert( knots, knots + nKnots );

            CHECK( KV1 == KV2 );
        }
    }

    TEST( remove )
    {
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV.remove( KV.ubegin() + 3 );
            real_t corrKnots[] = {0, 0, 0, .25, .25, .5, .5, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {3, 5, 7, 8, 11};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV.remove( KV.ubegin(), 3 );
            real_t corrKnots[] = {.25, .25, .5, .5, .75, .75, 1, 1, 1};
            real_t corrUKnots[] = {.25, .5, .75, 1};
            mult_t corrEndPos[] = {2, 4, 6, 9};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV.remove( KV.ubegin()+1, 7 );
            real_t corrKnots[] = {0, 0, 0, .5, .5, .75, .75, 1, 1, 1};
            real_t corrUKnots[] = {0, .5, .75, 1};
            mult_t corrEndPos[] = {3, 5, 7, 10};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV.remove( KV.uend()-1, 2 );
            real_t corrKnots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1};
            real_t corrUKnots[] = {0, .25, .5, .75, 1};
            mult_t corrEndPos[] = {3, 5, 7, 9, 10};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV.remove( KV.uend()-1, 3 );
            KV.remove( KV.ubegin(), 3 );
            real_t corrKnots[] = {.25, .25, .5, .5, .75, .75};
            real_t corrUKnots[] = {.25, .5, .75};
            mult_t corrEndPos[] = {2, 4, 6};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        // Now test also the remove with value.
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV.remove( .75, 3 );
            real_t corrKnots[] = {0, 0, 0, .25, .25, .5, .5, 1, 1, 1};
            real_t corrUKnots[] = {0, .25, .5, 1};
            mult_t corrEndPos[] = {3, 5, 7, 10};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        // Removing a value that is not present should not do anything.
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV0(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            gsKnotVector<real_t> KV1(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
            KV0.remove( .3, 3 );
            CHECK( KV0 == KV1 );
        }
    }
    TEST ( operatorStdVector )
    {
        real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(2, knots, knots+sizeof(knots)/sizeof(real_t));  //( 0, 1, 3, 3, 2 );
        std::vector<real_t> corrVec( knots, knots + sizeof(knots) / sizeof(real_t) );
        std::vector<real_t> vec = KV;
        CHECK( vec == corrVec );
    }
    TEST( reverse )
    {
        {
            real_t ak[] = {0, 0, 0, .2, .5, .7, 2, 2, 2};
            gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
            KV.reverse();
            real_t corrKnots[] = {0, 0, 0, 1.3, 1.5, 1.8, 2, 2, 2};
            real_t corrUKnots[] = {0, 1.3, 1.5, 1.8, 2};
            mult_t corrEndPos[] = { 3, 4, 5, 6, 9 };
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
        {
            real_t ak[] = {0, 0, 0, .2, .2, .5, .7, .7, .7, 2, 2, 2};
            gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
            KV.reverse();
            real_t corrKnots[] = {0, 0, 0, 1.3, 1.3, 1.3, 1.5, 1.8, 1.8, 2, 2, 2};
            real_t corrUKnots[] = {0, 1.3, 1.5, 1.8, 2};
            mult_t corrEndPos[] = {3, 6, 7, 9, 12};
            compareUKVs( KV, corrKnots, corrUKnots, corrEndPos, 2 );
        }
    }
    // TODO: Unittest for affineTransformTo.
    TEST( multiplicity )
    {
        {
            real_t ak[] = {0, 0, 0, .2, .5, .7, .7, 1, 1, 1};
            gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
            mult_t corrMult0[] = { 0, 3,  1,  0,  1,  2, 3, 0 };
            real_t uinputs[]   = {-1, 0, .2, .3, .5, .7, 1, 1.1 };
            for( size_t i = 0; i < sizeof(corrMult0)/sizeof(mult_t); ++i )
                CHECK( corrMult0[i] == KV.multiplicity(uinputs[i]) );

            mult_t corrMult1[] = {3, 3, 3, 1, 1, 2, 2, 3, 3, 3};
            for( size_t i = 0; i < sizeof(corrMult1)/sizeof(mult_t); ++i )
                CHECK( corrMult1[i] == KV.multiplicityIndex(i) );
        }
        {
            real_t knots[] = {0, 0, 0, .2, .2, .4, .4, .6, .6, .8, .8, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));  //( 0, 1, 4, 3, 2 );
            mult_t corrMult[] = { 3, 2, 2, 2, 2, 3 };
            mult_t i = 0;
            for( uniqIter uit = KV.ubegin(); uit != KV.uend(); ++uit, ++i )
                CHECK( uit.multiplicity() == corrMult[i] );
        }
    }
    TEST( index )
    {
        real_t ak[] = {0, 0, 0, .2, .5, .7, .7, 1, 1, 1};
        gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
        mult_t corrFirst[] = { 0, 3, 4, 5, 7 };
        mult_t corrLast[]  = { 3, 4, 5, 7, 10 };
        mult_t i = 0;
        for( uniqIter uit = KV.ubegin(); uit != KV.uend(); ++uit, ++i )
        {
            CHECK( uit.firstAppearance() == corrFirst[i] );
            CHECK( uit.multSum()  == corrLast[i] );
        }
    }
    TEST( isConsistent )
    {
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t)); //(0,1,3,3,2);
            CHECK( KV.check() );
        }
        {
            real_t ak[] = {0, 0, 0, .2, .5, .7, .7, 2, 2, 2};
            gsKnotVector<real_t> KV(2, ak, ak + sizeof(ak)/sizeof(real_t));
            CHECK( KV.check() );
        }
        /*{
            // Non-sorted knots. Constructor throws in debug mode.
            // In release mode, isConsistent() has to return false.
            real_t ak[] = {0, 0, 0, .2, .5, .7, .6, 2, 2, 2};
            std::vector<real_t> knots(ak, ak + sizeof(ak)/sizeof(real_t));
            #ifndef NDEBUG
                CHECK_THROW( gsKnotVector<real_t> KV(-1,knots),
                             std::runtime_error);
            #else
            gsKnotVector<real_t> KV(ak, ak + sizeof(ak)/sizeof(real_t));
                CHECK(KV.check()==false);
            #endif
        }*/
        {
            real_t ak[] = {0, 0, 0, .2, .5, .7, .7, 2, 2, 2};
            std::vector<real_t> aknots(ak, ak + sizeof(ak)/sizeof(real_t));
            mult_t am[] = {2,1,4,6,9};
            std::vector<mult_t> amults(am, am + sizeof(am)/sizeof(mult_t));
            CHECK( gsKnotVector<real_t>::isConsistent( aknots,
                                                              amults) == false );
        }
    }
    TEST( sizes )
    {
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));  //(0,1,3,3,2);
            CHECK( KV.size() == 12 );
            CHECK( KV.uSize() == 5 );
            CHECK( KV.numElements() == 4 );
        }
        {
            real_t knots[] = {0, .2, .4, .6, .8, 1};
            gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));
            CHECK( KV.size() == 6 );
            CHECK( KV.uSize() == 6 );
            CHECK( KV.numElements() == 1 );
        }
    }
    TEST( swap ) // Checks the copy-constructor as well by the way.
    {
        real_t knotsa[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
        real_t knotsb[] = {0, 0, .2, .4, .6, .8, 1, 1};
        gsKnotVector<real_t> KVa(2, knotsa, knotsa + sizeof(knotsa) / sizeof(real_t));  //(0,1,3,3,2);
        gsKnotVector<real_t> KVb(1, knotsb, knotsb + sizeof(knotsb) / sizeof(real_t));//(0,1,4,2,1);
        gsKnotVector<real_t> KVc(KVa);
        gsKnotVector<real_t> KVd(KVb);
        CHECK( KVc == KVa );
        CHECK( KVd == KVb );

        KVa.swap(KVb);
        CHECK( KVa == KVd );
        CHECK( KVb == KVc );
    }

    TEST(iteratorAssignement)
    {
        gsKnotVector<real_t> KV(-3,1,2,1,2);
        gsKnotVector<real_t>::iterator testIt;
        testIt=KV.begin();
        CHECK(*testIt==-3);
        gsKnotVector<real_t>::uiterator UtestIt;
        UtestIt=KV.ubegin();
        CHECK(*UtestIt==-3);

        internal::gsKnotIterator<real_t> aa(KV,0);
    }

    TEST( beginAt_and_endAt )
    {
        typedef gsKnotVector<real_t>::iterator knotIter;

        real_t knots[] = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        knotIter corrBegin, corrEnd;
        uniqIter uit;

        // Test 0: *uit == 0
        uit = KV.ubegin();
        corrBegin = KV.begin();
        corrEnd   = KV.begin();
        corrEnd += 4;
        CHECK( KV.beginAt(uit.uIndex()) == corrBegin );
        CHECK( KV.endAt  (uit.uIndex()) == corrEnd   );

        // Test 1: *uit == .25
        corrBegin = corrEnd;
        corrEnd += 3;
        ++uit;
        CHECK( KV.beginAt(uit.uIndex()) == corrBegin );
        CHECK( KV.endAt  (uit.uIndex()) == corrEnd   );

        // Test 2: *uit == .5
        corrBegin = corrEnd;
        corrEnd += 1;
        ++uit;
        CHECK( KV.beginAt(uit.uIndex()) == corrBegin );
        CHECK( KV.endAt  (uit.uIndex()) == corrEnd   );

        // Test 3: *uit == .75
        corrBegin = corrEnd;
        corrEnd += 2;
        ++uit;
        CHECK( KV.beginAt(uit.uIndex()) == corrBegin );
        CHECK( KV.endAt  (uit.uIndex()) == corrEnd   );

        // Test 4: *uit == 1
        corrBegin = corrEnd;
        corrEnd += 3;
        ++uit;
        CHECK( KV.beginAt(uit.uIndex()) == corrBegin );
        CHECK( KV.endAt  (uit.uIndex()) == corrEnd   );
        CHECK( KV.endAt  (uit.uIndex()) == KV.end()  );
    }

}

SUITE(gsKnotVectors_test_uniqIter)
{
    typedef gsKnotVector<real_t>::uiterator uniqIter;
    TEST(firstAppearance_and_multSum)
    {
        real_t knots[] = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        uniqIter uit = KV.ubegin();
        // *uit == 0
        CHECK( uit.firstAppearance() == 0 );
        CHECK( uit.multSum()   == 4 );

        ++uit;
        // *uit == .25
        CHECK( uit.firstAppearance() == 4 );
        CHECK( uit.multSum()   == 7 );

        uit += 2;
        // *uit == .75
        CHECK( uit.firstAppearance() == 8 );
        CHECK( uit.multSum()   == 10 );

        uit = KV.uend();
        --uit;
        // *uit == 1
        CHECK( uit.firstAppearance() == 10 );
        CHECK( uit.multSum()   == 13 );

        uit -= 2;
        // *uit == .5
        CHECK( uit.firstAppearance() == 7 );
        CHECK( uit.multSum()   == 8 );
    }

    // Note that the multiplicity is tested at the same spot as
    // gsKnotVector<T>::multiplicity(...).

    TEST( increment_and_decrement )
    {
        real_t knots[] = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));
        real_t uknots[] = {0, .25, .5, .75, 1};

        unsigned i = 0;
        uniqIter uit;

        //=================//
        // ++uit and --uit //
        //=================//

        for( uit = KV.ubegin(); uit != KV.uend(); ++uit, ++i )
            CHECK( *uit == uknots[i] );

        --i;
        for( uit = KV.uend() - 1; uit != KV.ubegin(); --uit, --i )
            CHECK( *uit == uknots[i] );
        // Zero in the beginning is checked separately.
        CHECK( *uit == knots[i]);

        //=================//
        // uit++ and uit-- //
        //=================//

        for( uit = KV.ubegin(); uit != KV.uend(); uit++, ++i )
            CHECK( *uit == uknots[i] );

        --i;
        for( uit = KV.uend() - 1; uit != KV.ubegin(); uit--, --i )
            CHECK( *uit == uknots[i] );
        // Zero again separately.
        CHECK( *uit == knots[i]);

        //============//
        // += and -=  //
        //============//

        uit = KV.ubegin();
        uit += 2;
        CHECK( *uit == .5 );
        uit += 2;
        CHECK( *uit == 1 );
        uit -= 3;
        CHECK( *uit == .25 );
        uit -= 1;
        CHECK( *uit == 0 );
        CHECK( uit == KV.ubegin() );
        uit += 5;
        CHECK( uit == KV.uend() );

        //=========//
        // + and - //
        //=========//

        uit = KV.ubegin();
        CHECK( *(uit +  0 ) ==  0  );
        CHECK( *(uit +  1 ) == .25 );
        CHECK( *(uit +  2 ) == .5  );
        CHECK( *(uit +  3 ) == .75 );
        CHECK( *(uit +  4 ) ==  1  );

        CHECK( *(uit -  0 ) ==  0  );
        CHECK( *(uit - -1 ) == .25 );
        CHECK( *(uit - -2 ) == .5  );
        CHECK( *(uit - -3 ) == .75 );
        CHECK( *(uit - -4 ) ==  1  );

        uit = KV.uend();
        CHECK( *(uit -  1 ) ==  1  );
        CHECK( *(uit -  2 ) == .75 );
        CHECK( *(uit -  3 ) == .5  );
        CHECK( *(uit -  4 ) == .25 );
        CHECK( *(uit -  5 ) ==  0  );

        CHECK( *(uit + -1 ) ==  1  );
        CHECK( *(uit + -2 ) == .75 );
        CHECK( *(uit + -3 ) == .5  );
        CHECK( *(uit + -4 ) == .25 );
        CHECK( *(uit + -5 ) ==  0  );

        //=====//
        // [a] //
        //=====//

        uit = KV.ubegin();
        CHECK(uit[ 0] ==  0  );
        CHECK(uit[ 1] == .25 );
        CHECK(uit[ 2] == .5  );
        CHECK(uit[ 3] == .75 );
        CHECK(uit[ 4] ==  1  );

        ++uit;
        CHECK(uit[-1] ==  0  );
        CHECK(uit[ 0] == .25 );
        CHECK(uit[ 1] == .5  );
        CHECK(uit[ 2] == .75 );
        CHECK(uit[ 3] ==  1  );

        uit = KV.uend();
        CHECK(uit[-1] ==  1  );
        CHECK(uit[-2] == .75 );
        CHECK(uit[-3] == .5  );
        CHECK(uit[-4] == .25 );
        CHECK(uit[-5] ==  0  );

        //====//
        // -> //
        //====//

        uit = KV.ubegin();
        CHECK( *(uit.operator->()) ==  0  );
        ++uit;
        CHECK( *(uit.operator->()) == .25 );
        uit = KV.uend() - 1;
        CHECK( *(uit.operator->()) ==  1  );
    }

    // Checks that smaller is really smaller than bigger.
    void uniqIterComparison_test( uniqIter smaller, uniqIter bigger )
    {
        CHECK( smaller == smaller );
        CHECK( smaller != bigger  );
        CHECK( smaller <= bigger  );
        CHECK( smaller <  bigger  );
        CHECK( bigger  >= smaller );
        CHECK( bigger  >  smaller );
        CHECK( bigger  != smaller );
        CHECK( bigger  == bigger  );
    }

    TEST( comparisons )
    {
        real_t knots[] = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        uniqIter uit0 = KV.ubegin();
        ++uit0;
        uniqIter uit1 = KV.uend();
        --uit1;

        uniqIterComparison_test( KV.ubegin(), KV.uend() );
        uniqIterComparison_test( KV.ubegin(), uit0 );
        uniqIterComparison_test( KV.ubegin(), uit1 );
        uniqIterComparison_test( uit0, uit1 );
        uniqIterComparison_test( uit0, KV.uend() );
        uniqIterComparison_test( uit1, KV.uend() );
    }

    TEST( uIndex )
    {
        typedef gsKnotVector<real_t>::mult_t mult_t;

        real_t knots[] = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        mult_t i = 0;
        for( uniqIter uit = KV.ubegin(); uit != KV.uend(); ++uit, ++i )
        {
            CHECK( uit.uIndex() == i );
        }
    }

}

SUITE(gsKnotVectors_test_gsKnotIterator)
{
    typedef gsKnotVector<real_t>::smart_iterator iter;
    typedef gsKnotVector<real_t>::mult_t mult_t;
    TEST( functionality )
    {
        real_t corrUIndices[]         = {0, 0, 0, 0,  1,   1,   1,   2,  3,   3,  4,  4,  4 };
        real_t corrMults[]            = {4, 4, 4, 4,  3,   3,   3,   1,  2,   2,  3,  3,  3 };
        real_t corrFirstAppearances[] = {0, 0, 0, 0,  4,   4,   4,   7,  8,   8,  10, 10, 10};
        real_t corrLastAppearances[]  = {3, 3, 3, 3,  6,   6,   6,   7,  9,   9,  12, 12, 12};
        real_t knots[]                = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1,  1,  1 };
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        mult_t i = 0;
        iter it = KV.sbegin();
        for( ; it != KV.send(); ++it, ++i )
        {
            CHECK( it.index()           == i                       );
            CHECK( it.uIndex()          == corrUIndices[i]         );
            CHECK( it.value()           == knots[i]                );
            CHECK( it.multiplicity()    == corrMults[i]            );
            CHECK( it.firstAppearance() == corrFirstAppearances[i] );
            CHECK( it.lastAppearance()  == corrLastAppearances[i]  );
        }

        it = KV.send() - 1;
        i = sizeof(knots)/sizeof(real_t) - 1;
        for( ; it > KV.sbegin(); --it, --i )
        {
            CHECK( it.index()           == i                       );
            CHECK( it.uIndex()          == corrUIndices[i]         );
            CHECK( it.value()           == knots[i]                );
            CHECK( it.multiplicity()    == corrMults[i]            );
            CHECK( it.firstAppearance() == corrFirstAppearances[i] );
            CHECK( it.lastAppearance()  == corrLastAppearances[i]  );
        }
        CHECK( it.index()           == 0                       );
        CHECK( it.uIndex()          == corrUIndices[0]         );
        CHECK( it.value()           == knots[0]                );
        CHECK( it.multiplicity()    == corrMults[0]            );
        CHECK( it.firstAppearance() == corrFirstAppearances[0] );
        CHECK( it.lastAppearance()  == corrLastAppearances[0]  );
    }

    // remains: +-, comparisons and backToFirstAppearance, uIndex+-.
    TEST( uPrev_and_uNext )
    {
        real_t uKnots[] = {0, .25, .5, .75, 1};
        real_t knots[]  = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        size_t i = 0;
        iter it;
        for( it = KV.sbegin(); it != KV.send(); it.uNext(), ++i )
            CHECK( *it == uKnots[i] );

        --i;
        for( it = KV.send() - 1; it != KV.sbegin(); it.uPrev(), --i )
            CHECK( *it == uKnots[i] );
        CHECK( *it == uKnots[0]);
    }

    TEST( increments_and_decrements )
    {
        mult_t corrUIndices[] = {0, 0, 0, 0,  1,   1,   1,   2,  3,   3,   4,  5, 5, 5};
        real_t knots[]        = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, .85, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        // ++it
        iter it = KV.sbegin();
        size_t i = 0;
        for( ; it < KV.send(); ++it, ++i )
        {
            CHECK( *it == knots[i] );
            CHECK( it.uIndex() == corrUIndices[i] );
        }

        // --it
        it = KV.send() - 1;
        i  = sizeof(knots)/sizeof(real_t) - 1;
        for( ; it > KV.sbegin(); --it, --i )
        {
            CHECK( *it == knots[i] );
            CHECK( it.uIndex() == corrUIndices[i] );
        }
        CHECK( *it == knots[0] );
        CHECK( it.uIndex() == corrUIndices[0] );

        // it++
        it = KV.sbegin();
        i  = 0;
        for( ; it < KV.send(); it++, ++i )
        {
            CHECK( *it == knots[i] );
            CHECK( it.uIndex() == corrUIndices[i] );
        }

        // it--
        it = KV.send() - 1;
        i  = sizeof(knots)/sizeof(real_t) - 1;
        for( ; it > KV.sbegin(); it--, --i )
        {
            CHECK( *it == knots[i] );
            CHECK( it.uIndex() == corrUIndices[i] );
        }
        CHECK( *it == knots[0] );
        CHECK( it.uIndex() == corrUIndices[0] );

        // +=
        it = KV.sbegin();
        i  = 0; //test assumes an even number of knots
        for( ; it < KV.send(); it += 2, i+= 2 )
        {
            CHECK( *it == knots[i] );
            CHECK( it.uIndex() == corrUIndices[i] );
        }

        // -=
        it = KV.send()-1;
        i  = sizeof(knots)/sizeof(real_t) - 1;
        for( ; it > KV.sbegin()+1; it -= 2, i -= 2 )
        {
            CHECK( *it == knots[i] );
            CHECK( it.uIndex() == corrUIndices[i] );
        }
    }

    TEST( plus_and_minus )
    {
        real_t knots[]  = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));
        ptrdiff_t size = KV.size();

        // If you are confused about the for cycles,
        // cf. uniqIter::increment_and_decrement, where they are done explicitely.

        iter it = KV.sbegin();
        ptrdiff_t i = 0;

        // + with a positive argument
        for( ; i < size; ++i )
            CHECK( *(it +  i ) ==  knots[i]  );

        // + with a negative argument
        it = KV.send() - 1;
        i = 0;
        for( ; i > -size; --i )
            CHECK( *(it + i) == knots[size+i-1] );

        // - with a positive argument
        it = KV.send() - 1;
        i = size  - 1;
        for( ; i >= 0; --i )
            CHECK( *(it -  i ) ==  knots[size-i-1] );

        // - with a negative argument
        it = KV.sbegin();
        i = 0;
        for( ; i > -size; --i )
            CHECK( *(it - i ) == knots[-i] );

        // + and - from the middle
        ptrdiff_t halfSize = size / 2;
        it = KV.sbegin() + halfSize;
        i = - halfSize;
        for( ; i + halfSize < size; ++i )
            CHECK( *(it + i ) == knots[i+halfSize] );

        // [] from the beginning
        it = KV.sbegin();
        i = 0;
        for( ; i < size; ++i )
            CHECK( it[i] == knots[i] );

        // [] from inside
        ++it;
        i = -1;
        for( ; i < size - 1; ++i )
            CHECK( it[i] == knots[i+1] );

        // [] from the end
        it = KV.send();
        i = size;
        for( ; i > 0; --i )
            CHECK( it[-i] == knots[size-i]);

        // TODO Improve by also checking throwing if asked for it[0] or it[-size-1].
    }

    TEST( ptr_difference )
    {
        real_t knots[]  = {0, 0, 0, 0, .25, .25, .25, .5, .75, .75, 1, 1, 1};
        gsKnotVector<real_t> KV(3, knots, knots+sizeof(knots)/sizeof(real_t));

        iter it0 = KV.sbegin();
        iter it1 = KV.sbegin();
        CHECK( it1 - it0 ==  0 );
        ++it0;
        CHECK( it1 - it0 == -1 );
        CHECK( it0 - it1 ==  1 );
        it1 += 3;
        CHECK( it1 - it0 ==  2 );
        CHECK( it0 - it1 == -2 );
        it0 = KV.send();
        CHECK( it1 - it0 == -10 );
        CHECK( it0 - it1 ==  10 );
        it1 = KV.send();
        CHECK( it1 - it0 ==  0 );
        CHECK( it0 - it1 ==  0 );
        it0= KV.sbegin();
        CHECK( it1 - it0 ==  13 );
        CHECK( it0 - it1 == -13 );
    }

    TEST( uniqueAndBreaks )
    {
        {
            real_t knots[] = {0, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1};
            gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));
            std::vector<real_t> unique = KV.unique();
            std::vector<real_t> breaks = KV.breaks();
            real_t corrUnique[] = {0, .25, .5, .75, 1};
            real_t corrBreaks[] = {0, .25, .5, .75, 1};

            CHECK( unique.size() == 5 );
            CHECK( breaks.size() == 5 );
            for( size_t i = 0; i < 5; ++i )
            {
                CHECK( unique[i] == corrUnique[i] );
                CHECK( breaks[i] == corrBreaks[i] );
            }
        }
        {
            real_t knots[] = {0, .2, .4, .6, .8, 1};
            gsKnotVector<real_t> KV(2, knots, knots + sizeof(knots)/sizeof(real_t));
            std::vector<real_t> unique = KV.unique();
            std::vector<real_t> breaks = KV.breaks();
            real_t corrUnique[] = {0, .2, .4, .6, .8, 1};

            CHECK( unique.size() == 6 );
            CHECK( breaks.size() == 2 );
            for( size_t i = 0; i < 6; ++i )
                CHECK( unique[i] == corrUnique[i] );
        }
    }
}
