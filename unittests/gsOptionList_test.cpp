/** @file gsOptionList_test.cpp

    @brief test gsIO/gsOptionList

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris,  H. Weiner, J. Vogl
 **/

#include "gismo_unittest.h"

const std::string LABEL_STR_1 = "LABEL_STR_1";
const std::string LABEL_STR_2 = "LABEL_STR_2";
const std::string LABEL_INT_1 = "LABEL_INT_1";
const std::string LABEL_INT_2 = "LABEL_INT_2";
const std::string LABEL_REAL_1 = "LABEL_REAL_1";
const std::string LABEL_REAL_2 = "LABEL_REAL_2";
const std::string LABEL_BOOL_1 = "LABEL_BOOL_1";
const std::string LABEL_BOOL_2 = "LABEL_BOOL_2";
const std::string DESC_1 = "DESC_1";
const std::string DESC_2 = "DESC_1";
const std::string STR_ANY = "ANY";
const std::string STR_1 = "STR_1";
const std::string STR_2 = "STR_2";
const int INT_ANY = -1;
const int INT_1 = 1;
const int INT_2 = 2;
const real_t REAL_ANY = -1.0;
const real_t REAL_1 = 1.0;
const real_t REAL_2 = 2.0;
const bool BOOL_ANY = false;
const bool BOOL_1 = true;
const bool BOOL_2 = false;

gsOptionList loadFromFile(std::string& path)
{
    std::string input(path);
    gsOptionList result;
    gsReadFile<>(input, result);
    return result;
}

void saveToFile(std::string& path, gsOptionList& list)
{
    std::string fullPath(path);
    gsFileData<> write;
    write << list;
    write.dump(fullPath);
}

void checkAssemblerOptions(gsOptionList& myList)
{
    CHECK_EQUAL(8u, myList.size());
    CHECK_EQUAL(11, myList.getInt("DirichletStrategy"));
    CHECK_EQUAL(101, myList.getInt("DirichletValues"));
    CHECK_EQUAL(1, myList.getInt("InterfaceStrategy"));
    CHECK_EQUAL(1, myList.getInt("bdB"));
    CHECK_EQUAL(1, myList.getInt("quB"));
    CHECK_EQUAL(2.0, myList.getReal("bdA"));
    #ifndef GISMO_WITH_GMP
        CHECK_EQUAL(0.333, myList.getReal("bdO"));
    #else
        CHECK_CLOSE(0.333, myList.getReal("bdO"), 1e-3);
    #endif
    CHECK_EQUAL(1.0, myList.getReal("quA"));
}

SUITE(gsOptionList_test)
{

    /***
     * test what we get if the gsOptionList is empty
     */
    TEST(get_items_from_empty_list)
    {
        gsOptionList myList;
        CHECK_EQUAL(0u, myList.size());
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_1));
        CHECK_EQUAL(INT_1, myList.askInt(LABEL_INT_1, INT_1));
        CHECK_EQUAL(REAL_1, myList.askReal(LABEL_REAL_1, REAL_1));
        CHECK_EQUAL(BOOL_1, myList.askSwitch(LABEL_BOOL_1, BOOL_1));
        CHECK_EQUAL(BOOL_1, myList.askSwitch(LABEL_BOOL_1, BOOL_1));
        CHECK_EQUAL(STR_2, myList.askString(LABEL_STR_1, STR_2));
        CHECK_EQUAL(STR_2, myList.askString(LABEL_STR_1, STR_2));
    }

    /***
     * test what we get if the gsOptionList is empty and getString is invoked.
     * should throw an exception!
     */
    TEST(get_items_from_empty_list_and_throw_exception)
    {
        gsOptionList myList;
        CHECK_EQUAL(0u, myList.size());
        UnitTest::deactivate_output();
        CHECK_THROW_IN_DEBUG(myList.getString(LABEL_STR_1), std::exception);
        UnitTest::reactivate_output();
    }

    /***
     * test adding and asking/getting strings
     */
    TEST(add_and_get_something)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        CHECK_EQUAL(1u, myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_ANY));
        CHECK_EQUAL(STR_2, myList.askString(LABEL_STR_2, STR_2));
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_ANY));
        CHECK_EQUAL(STR_2, myList.askString(LABEL_STR_2, STR_2));
        CHECK_EQUAL(size_t(1), myList.size());
    }

    /***
     * testing add/ask/get for multiple types
     */
    TEST(add_more_and_get_something)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        myList.addString(LABEL_STR_2, DESC_2, STR_2);
        myList.addInt(LABEL_INT_1, DESC_1, INT_1);
        myList.addInt(LABEL_INT_2, DESC_2, INT_2);
        myList.addReal(LABEL_REAL_1, DESC_1, REAL_1);
        myList.addReal(LABEL_REAL_2, DESC_2, REAL_2);
        myList.addSwitch(LABEL_BOOL_1, DESC_1, BOOL_1);
        myList.addSwitch(LABEL_BOOL_2, DESC_2, BOOL_2);
        CHECK_EQUAL(size_t(8), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_ANY));
        CHECK_EQUAL(size_t(8), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_ANY));
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_2));
        CHECK_EQUAL(STR_2, myList.askString(LABEL_STR_2, STR_ANY));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(INT_1, myList.askInt(LABEL_INT_1, INT_ANY));
        CHECK_EQUAL(INT_2, myList.getInt(LABEL_INT_2));
        CHECK_EQUAL(INT_2, myList.askInt(LABEL_INT_2, INT_ANY));
        CHECK_EQUAL(REAL_1, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(REAL_1, myList.askReal(LABEL_REAL_1, REAL_ANY));
        CHECK_EQUAL(REAL_2, myList.getReal(LABEL_REAL_2));
        CHECK_EQUAL(REAL_2, myList.askReal(LABEL_REAL_2, REAL_ANY));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_1));
        CHECK_EQUAL(BOOL_1, myList.askSwitch(LABEL_BOOL_1, BOOL_ANY));
        CHECK_EQUAL(BOOL_2, myList.getSwitch(LABEL_BOOL_2));
        CHECK_EQUAL(BOOL_2, myList.askSwitch(LABEL_BOOL_2, BOOL_ANY));
        CHECK_EQUAL(size_t(8), myList.size());
    }

    /***
     * test addString and then invoke askSwitch for the same label
     */
    TEST(add_and_get_of_different_type)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        CHECK_EQUAL(size_t(1), myList.size());
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_ANY));
        CHECK_EQUAL(INT_1, myList.askInt(LABEL_STR_1, INT_1));
        CHECK_EQUAL(REAL_1, myList.askReal(LABEL_STR_1, REAL_1));
        CHECK_EQUAL(BOOL_1, myList.askSwitch(LABEL_STR_1, BOOL_1));
        CHECK_EQUAL(STR_1, myList.askString(LABEL_STR_1, STR_ANY));
        CHECK_EQUAL(size_t(1), myList.size());
    }

    /***
     * test what happens when we invoke addString for the same label twice
     */
    TEST(add_and_add)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        myList.addInt(LABEL_INT_1, DESC_1, INT_1);
        myList.addReal(LABEL_REAL_1, DESC_1, REAL_1);
        myList.addSwitch(LABEL_BOOL_1, DESC_1, BOOL_1);
        CHECK_EQUAL(size_t(4), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(REAL_1, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_1));
        myList.addString(LABEL_STR_1, DESC_1, STR_2);
        myList.addInt(LABEL_INT_1, DESC_1, INT_2);
        myList.addReal(LABEL_REAL_1, DESC_1, REAL_2);
        myList.addSwitch(LABEL_BOOL_1, DESC_1, BOOL_2);
        CHECK_EQUAL(size_t(4), myList.size());
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_1));
        UnitTest::deactivate_output();
        CHECK_THROW(myList.addString(LABEL_INT_1, DESC_1, STR_2),
                std::runtime_error);
        CHECK_THROW(myList.addInt(LABEL_REAL_1, DESC_1, INT_1), std::runtime_error);
        CHECK_THROW(myList.addReal(LABEL_BOOL_1, DESC_1, REAL_1),
                std::runtime_error);
        CHECK_THROW(myList.addSwitch(LABEL_STR_1, DESC_1, BOOL_1),
                std::runtime_error);
        UnitTest::reactivate_output();
        CHECK_EQUAL(size_t(4), myList.size());
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(size_t(4), myList.size());
    }

    /***
     * test add/get/set/get
     */
    TEST(add_get_set_get)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        myList.addString(LABEL_STR_2, DESC_2, STR_2);
        myList.addInt(LABEL_INT_1, DESC_1, INT_1);
        myList.addInt(LABEL_INT_2, DESC_2, INT_2);
        myList.addReal(LABEL_REAL_1, DESC_1, REAL_1);
        myList.addReal(LABEL_REAL_2, DESC_2, REAL_2);
        myList.addSwitch(LABEL_BOOL_1, DESC_1, BOOL_1);
        myList.addSwitch(LABEL_BOOL_2, DESC_2, BOOL_2);
        CHECK_EQUAL(size_t(8), myList.size());
        myList.setString(LABEL_STR_1, STR_1);
        myList.setInt(LABEL_INT_1, INT_1);
        myList.setReal(LABEL_REAL_1, REAL_1);
        myList.setSwitch(LABEL_BOOL_1, BOOL_1);
        CHECK_EQUAL(size_t(8), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_2));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(INT_2, myList.getInt(LABEL_INT_2));
        CHECK_EQUAL(REAL_1, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(REAL_2, myList.getReal(LABEL_REAL_2));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_1));
        CHECK_EQUAL(BOOL_2, myList.getSwitch(LABEL_BOOL_2));

        myList.setString(LABEL_STR_1, STR_2);
        myList.setInt(LABEL_INT_2, INT_1);
        myList.setReal(LABEL_REAL_1, REAL_2);
        myList.setSwitch(LABEL_BOOL_2, BOOL_1);
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_2));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_2));
        CHECK_EQUAL(REAL_2, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(REAL_2, myList.getReal(LABEL_REAL_2));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_1));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_2));

        myList.setString(LABEL_STR_2, STR_1);
        myList.setInt(LABEL_INT_1, INT_2);
        myList.setReal(LABEL_REAL_2, REAL_1);
        myList.setSwitch(LABEL_BOOL_1, BOOL_2);
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_2));
        CHECK_EQUAL(INT_2, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_2));
        CHECK_EQUAL(REAL_2, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(REAL_1, myList.getReal(LABEL_REAL_2));
        CHECK_EQUAL(BOOL_2, myList.getSwitch(LABEL_BOOL_1));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_2));

        myList.setString(LABEL_STR_1, STR_1);
        myList.setString(LABEL_STR_2, STR_2);
        myList.setInt(LABEL_INT_1, INT_1);
        myList.setInt(LABEL_INT_2, INT_2);
        myList.setReal(LABEL_REAL_1, REAL_1);
        myList.setReal(LABEL_REAL_2, REAL_2);
        myList.setSwitch(LABEL_BOOL_1, BOOL_1);
        myList.setSwitch(LABEL_BOOL_2, BOOL_2);
        CHECK_EQUAL(size_t(8), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(STR_2, myList.getString(LABEL_STR_2));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(INT_2, myList.getInt(LABEL_INT_2));
        CHECK_EQUAL(REAL_1, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(REAL_2, myList.getReal(LABEL_REAL_2));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_1));
        CHECK_EQUAL(BOOL_2, myList.getSwitch(LABEL_BOOL_2));
    }

    /***
     * Test what happens if we invoke set before add
     */
    TEST(set_before_add)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        myList.addInt(LABEL_INT_1, DESC_1, INT_1);
        myList.addReal(LABEL_REAL_1, DESC_1, REAL_1);
        myList.addSwitch(LABEL_BOOL_1, DESC_1, BOOL_1);
        CHECK_EQUAL(size_t(4), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        UnitTest::deactivate_output();
        CHECK_THROW(myList.setString(LABEL_STR_2, STR_2), std::runtime_error);
        CHECK_THROW(myList.setInt(LABEL_INT_2, INT_2), std::runtime_error);
        CHECK_THROW(myList.setReal(LABEL_REAL_2, REAL_2), std::runtime_error);
        CHECK_THROW(myList.setSwitch(LABEL_BOOL_2, BOOL_2), std::runtime_error);
        UnitTest::reactivate_output();
        CHECK_EQUAL(size_t(4), myList.size());
        CHECK_EQUAL(STR_1, myList.getString(LABEL_STR_1));
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_INT_1));
        CHECK_EQUAL(REAL_1, myList.getReal(LABEL_REAL_1));
        CHECK_EQUAL(BOOL_1, myList.getSwitch(LABEL_BOOL_1));
    }

    /***
     * Test getMultiInt
     */
    TEST(getMultiInt)
    {
        gsOptionList myList;
        int data[] = {5, 7, 4};
        myList.addInt("0", "", 5);
        myList.addInt("1", "", 7);
        myList.addInt("2", "", 4);
        myList.addInt("Size", "", 3);
        myList = myList.wrapIntoGroup("VEC");

        std::vector<index_t> vec = myList.getMultiInt("VEC");

        CHECK_EQUAL((size_t)3, vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            CHECK_EQUAL(data[i], vec[i]);
        }
    }

    /***
     * Test addMultiInt, based on getMultiInt
     */
    TEST(addMultiInt)
    {
        gsOptionList myList1, myList2;
        int data[] = {5, 7, 4};
        myList1.addInt("0", "", 5);
        myList1.addInt("1", "", 7);
        myList1.addInt("2", "", 4);
        myList1.addInt("Size", "", 3);
        myList1 = myList1.wrapIntoGroup("VEC");

        std::vector<index_t> vec(data, data + sizeof(data)/ sizeof(int));
        myList2.addMultiInt("VEC", "", vec);

        std::vector<index_t> vec1 = myList1.getMultiInt("VEC");
        std::vector<index_t> vec2 = myList2.getMultiInt("VEC");

        CHECK_EQUAL((size_t)3, vec.size());
        CHECK_EQUAL((size_t)3, vec1.size());
        CHECK_EQUAL((size_t)3, vec2.size());

        for (size_t i = 0; i < vec1.size(); ++i) {
            CHECK_EQUAL(data[i], vec[i]);
            CHECK_EQUAL(data[i], vec1[i]);
            CHECK_EQUAL(data[i], vec2[i]);
        }
    }

    /***
     * Test getMultiReal
     */
    TEST(getMultiReal)
    {
        gsOptionList myList;
        real_t data[] = {5.4, 7.3, 4.2};
        myList.addReal("0", "", 5.4);
        myList.addReal("1", "", 7.3);
        myList.addReal("2", "", 4.2);
        myList.addInt("Size", "", 3);
        myList = myList.wrapIntoGroup("VEC");

        std::vector<real_t> vec = myList.getMultiReal("VEC");

        CHECK_EQUAL((size_t)3, vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            CHECK_EQUAL(data[i], vec[i]);
        }
    }

    /***
     * test default options
     */
    TEST(test_default_options)
    {
        gsOptionList myList = gsAssembler<real_t>::defaultOptions();
        checkAssemblerOptions(myList);
    }

    /***
     * test loading from assembler_options.xml
     */
    TEST(load_from_assembler_options_xml)
    {
        std::string path = gsFileManager::findInDataDir( "options/assembler_options.xml" );
        gsOptionList myList = loadFromFile(path);
        checkAssemblerOptions(myList);
    }

    /***
     * test loading from optionlist.xml
     */
    TEST(load_from_optionlist_xml)
    {
        std::string path = gsFileManager::findInDataDir( "options/optionlist.xml" );
        gsOptionList myList = loadFromFile(path);
        CHECK_EQUAL(size_t(4), myList.size());

        CHECK_EQUAL(true, myList.getSwitch("boundary"));
        CHECK_EQUAL(11, myList.getInt("strategy"));
        CHECK_CLOSE((real_t)1/100000, myList.getReal("tolerance"),EPSILON);
        //CHECK_EQUAL((real_t)1/100000, myList.getReal("tolerance"));
        CHECK_EQUAL("/path/to/square.xml", myList.getString("geometry"));
    }

    /***
     * test saving/loading to/from test.xml
     */
    TEST(save_load_from_test_xml)
    {
        std::string path = gsFileManager::getTempPath()
            + gsFileManager::getNativePathSeparator() + util::to_string("test.xml");
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        myList.addInt(LABEL_INT_1, DESC_1, INT_1);
        myList.addReal(LABEL_REAL_1, DESC_1, REAL_1);
        myList.addSwitch(LABEL_BOOL_1, DESC_1, BOOL_1);
        saveToFile(path, myList);
        myList.setString(LABEL_STR_1, STR_2);
        myList.setInt(LABEL_INT_1, INT_2);
        myList.setReal(LABEL_REAL_1, REAL_2);
        myList.setSwitch(LABEL_BOOL_1, BOOL_2);

        gsOptionList myList2 = loadFromFile(path);
        CHECK_EQUAL(size_t(4), myList2.size());
        CHECK_EQUAL(STR_1, myList2.getString(LABEL_STR_1));
        CHECK_EQUAL(INT_1, myList2.getInt(LABEL_INT_1));
        CHECK_EQUAL(REAL_1, myList2.getReal(LABEL_REAL_1));
        CHECK_EQUAL(BOOL_1, myList2.getSwitch(LABEL_BOOL_1));

        // clean-up afterwards
        CHECK_EQUAL(0, remove(path.c_str()));
    }

    /***
     * test adding a value with the same label but a different type
     * -> should throw an exception
     */
    TEST(same_label_different_type)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        UnitTest::deactivate_output();
        CHECK_THROW(myList.addInt(LABEL_STR_1, DESC_1, INT_1), std::runtime_error);
        UnitTest::reactivate_output();
    }

    /***
     * test loading/saving/loading from/to/from assembler_options.xml/
     * assembler_options2.xml/assembler_options2.xml
     */
    TEST(loading_saving_loading_assembler_options)
    {
        // loading
        std::string path = gsFileManager::findInDataDir( "options/assembler_options.xml" );
        gsOptionList myList = loadFromFile(path);
        // saving
        std::string actual = gsFileManager::getTempPath()
            + gsFileManager::getNativePathSeparator() + util::to_string("assembler_options2.xml");
        saveToFile(actual, myList);
        // loading
        gsOptionList myList2 = loadFromFile(actual);
        checkAssemblerOptions(myList2);

        // clean-up afterwards
        CHECK_EQUAL(0, remove(actual.c_str()));
    }

    /***
     * test remove
     *
     */

    TEST(test_remove)
    {
        gsOptionList myList;
        myList.addString(LABEL_STR_1, DESC_1, STR_1);
        CHECK_EQUAL(size_t(1), myList.size());
        myList.remove(LABEL_STR_1);
        CHECK_EQUAL(size_t(0), myList.size());
        myList.addInt(LABEL_STR_1, DESC_1, INT_1);
        CHECK_EQUAL(size_t(1), myList.size());
        CHECK_EQUAL(INT_1, myList.getInt(LABEL_STR_1));
    }

    /***
     *  test groups
     *
     */

    TEST(test_groups)
    {
        gsOptionList myList;
        myList.addString("LABEL_STR_1", DESC_1, STR_1);
        myList.addInt("LABEL_INT_1", DESC_1, INT_1);
        myList.addReal("LABEL_REAL_1", DESC_1, REAL_1);
        myList.addSwitch("LABEL_BOOL_1", DESC_1, BOOL_1);
        CHECK_EQUAL(true,myList.hasGlobals());
        CHECK_EQUAL(false,myList.hasGroup("GR"));

        gsOptionList myList2 = myList.wrapIntoGroup("GR");
        CHECK_EQUAL(false, myList2.hasGlobals());
        CHECK_EQUAL(true, myList2.hasGroup("GR"));
        CHECK_EQUAL(size_t(4), myList2.size());
        myList2.addString("LABEL_STR_1", DESC_1, STR_2);
        myList2.addString("LABEL_STR_2", DESC_1, STR_2);
        CHECK_EQUAL(true, myList2.hasGlobals());
        CHECK_EQUAL(true, myList2.hasGroup("GR"));
        CHECK_EQUAL(size_t(6), myList2.size());
        CHECK_EQUAL(STR_2, myList2.getString("LABEL_STR_1"));
        CHECK_EQUAL(STR_2, myList2.getString("LABEL_STR_2"));
        CHECK_EQUAL(STR_1, myList2.getString("GR.LABEL_STR_1"));
        CHECK_EQUAL(INT_1, myList2.getInt("GR.LABEL_INT_1"));
        CHECK_EQUAL(REAL_1, myList2.getReal("GR.LABEL_REAL_1"));
        CHECK_EQUAL(BOOL_1, myList2.getSwitch("GR.LABEL_BOOL_1"));

        gsOptionList myList3 = myList2.getGroup("GR");
        CHECK_EQUAL(true, myList3.hasGlobals());
        CHECK_EQUAL(false, myList3.hasGroup("GR"));
        CHECK_EQUAL(size_t(4), myList3.size());
        CHECK_EQUAL(STR_1, myList3.getString("LABEL_STR_1"));
        CHECK_EQUAL(INT_1, myList3.getInt("LABEL_INT_1"));
        CHECK_EQUAL(REAL_1, myList3.getReal("LABEL_REAL_1"));
        CHECK_EQUAL(BOOL_1, myList3.getSwitch("LABEL_BOOL_1"));
        CHECK_EQUAL("undefined", myList3.askString("LABEL_STR_2","undefined"));

    }


}
