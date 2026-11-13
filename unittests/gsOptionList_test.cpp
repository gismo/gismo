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
    #ifndef gsGmp_ENABLED
        CHECK_EQUAL(0.333, myList.getReal("bdO"));
    #else
        CHECK_CLOSE(0.333, myList.getReal("bdO"), 1e-3);
    #endif
    CHECK_EQUAL(1.0, myList.getReal("quA"));
}

// Target list of the update() tests: entries in all four tables, with two ints
// so that an unknown/mismatching int label sorts between valid int labels.
static gsOptionList makeUpdateTarget()
{
    gsOptionList T;
    T.addString("S_a", "dS_a", "alpha");
    T.addInt   ("N_a", "dN_a", 1);
    T.addInt   ("N_z", "dN_z", 3);
    T.addReal  ("N_m", "dN_m", 2.25);
    T.addReal  ("R_a", "dR_a", 0.5);
    T.addSwitch("B_a", "dB_a", false);
    return T;
}

// Changes matching makeUpdateTarget() in label and type, in all four tables;
// all values and descriptions differ from the target's.
static gsOptionList makeValidChanges()
{
    gsOptionList V;
    V.addString("S_a", "other", "beta");
    V.addInt   ("N_a", "other", 10);
    V.addInt   ("N_z", "other", 30);
    V.addReal  ("R_a", "other", 0.75);
    V.addSwitch("B_a", "other", true);
    return V;
}

static bool hasLabel(const gsOptionList & list, const std::string & label)
{
    const std::vector<gsOptionList::OptionListEntry> e = list.getAllEntries();
    for (size_t i = 0; i != e.size(); ++i)
        if (e[i].label == label)
            return true;
    return false;
}

// Entry of the list with the given label (the label must exist).
static gsOptionList::OptionListEntry entryOf(const gsOptionList & list, const std::string & label)
{
    const std::vector<gsOptionList::OptionListEntry> e = list.getAllEntries();
    for (size_t i = 0; i != e.size(); ++i)
        if (e[i].label == label)
            return e[i];
    CHECK(false);
    return gsOptionList::OptionListEntry();
}

// Checks that two lists agree in size and, entry by entry, in type, label,
// description and exact value. Values are read through the typed getters, since
// getAllEntries() keeps only a truncated string of them.
static void checkListsEqual(const gsOptionList & a, const gsOptionList & b)
{
    CHECK_EQUAL(a.size(), b.size());
    const std::vector<gsOptionList::OptionListEntry> ea = a.getAllEntries();
    const std::vector<gsOptionList::OptionListEntry> eb = b.getAllEntries();
    CHECK_EQUAL(ea.size(), eb.size());
    for (size_t i = 0; i < ea.size() && i < eb.size(); ++i)
    {
        CHECK_EQUAL(ea[i].type,  eb[i].type);
        CHECK_EQUAL(ea[i].label, eb[i].label);
        CHECK_EQUAL(ea[i].desc,  eb[i].desc);
        const std::string & l = ea[i].label;
        if (ea[i].type == "string")
            CHECK_EQUAL(a.getString(l), b.getString(l));
        else if (ea[i].type == "int")
            CHECK_EQUAL(a.getInt(l), b.getInt(l));
        else if (ea[i].type == "real")
            CHECK_EQUAL(a.getReal(l), b.getReal(l));
        else if (ea[i].type == "bool")
            CHECK_EQUAL(a.getSwitch(l), b.getSwitch(l));
        else
            CHECK(false);
    }
}

// Runs T.update(other, mode) with the test output silenced; returns whether a
// runtime_error was thrown and stores its what() in msg.
static bool updateThrows(gsOptionList & T, const gsOptionList & other,
                         gsOptionList::updateType mode, std::string & msg)
{
    bool thrown = false;
    msg.clear();
    UnitTest::deactivate_output();
    try { T.update(other, mode); }
    catch (const std::runtime_error & e) { thrown = true; msg = e.what(); }
    UnitTest::reactivate_output();
    return thrown;
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

    /***
     * update() with ignoreIfUnknown (explicit and default argument): unknown
     * labels are skipped, known labels of the same type are applied
     */
    TEST(update_ignore_skips_unknown_label)
    {
        gsOptionList other;
        other.addString("S_a",   "other", "beta");
        other.addString("X_new", "other", "x");

        gsOptionList T1 = makeUpdateTarget();
        gsOptionList T2 = makeUpdateTarget();
        T1.update(other);
        T2.update(other, gsOptionList::ignoreIfUnknown);

        CHECK(!hasLabel(T1, "X_new"));
        CHECK(!hasLabel(T2, "X_new"));
        CHECK_EQUAL(size_t(6), T1.size());
        CHECK_EQUAL(size_t(6), T2.size());
        CHECK_EQUAL("beta", T1.getString("S_a"));
        CHECK_EQUAL("beta", T2.getString("S_a"));
    }

    /***
     * update() with ignoreIfUnknown: a label that exists under a different type
     * is not skipped but throws
     */
    TEST(update_ignore_throws_on_type_mismatch)
    {
        gsOptionList T = makeUpdateTarget();
        gsOptionList other;
        other.addInt("N_m", "other", 7);
        UnitTest::deactivate_output();
        CHECK_THROW(T.update(other, gsOptionList::ignoreIfUnknown), std::runtime_error);
        UnitTest::reactivate_output();
    }

    /***
     * update() with addIfUnknown: unknown labels are added together with their
     * description, known labels are applied
     */
    TEST(update_add_adds_unknown_label)
    {
        gsOptionList other;
        other.addReal  ("R_new", "dR_new", 1.5);
        other.addString("S_a",   "other",  "beta");

        gsOptionList T = makeUpdateTarget();
        T.update(other, gsOptionList::addIfUnknown);

        CHECK_EQUAL(1.5, T.getReal("R_new"));
        CHECK(hasLabel(T, "R_new"));
        if (hasLabel(T, "R_new"))
        {
            const gsOptionList::OptionListEntry e = entryOf(T, "R_new");
            CHECK_EQUAL("dR_new", e.desc);
            CHECK_EQUAL("real", e.type);
        }
        CHECK_EQUAL(size_t(7), T.size());
        CHECK_EQUAL("beta", T.getString("S_a"));
    }

    /***
     * update() with errorIfUnknown: an unknown label throws a runtime_error
     * whose message names the label
     */
    TEST(update_error_unknown_label_throws_naming_label)
    {
        gsOptionList T = makeUpdateTarget();
        gsOptionList other = makeValidChanges();
        other.addInt("N_q", "other", 5);

        bool thrown = false;
        std::string msg;
        UnitTest::deactivate_output();
        try { T.update(other, gsOptionList::errorIfUnknown); }
        catch (const std::runtime_error & e) { thrown = true; msg = e.what(); }
        UnitTest::reactivate_output();

        CHECK(thrown);
        CHECK(msg.find("N_q") != std::string::npos);
    }

    /***
     * update() with errorIfUnknown: a label present under another type throws
     * a runtime_error whose message names the label and both types
     */
    TEST(update_error_type_mismatch_throws_naming_label)
    {
        gsOptionList T = makeUpdateTarget();
        gsOptionList other = makeValidChanges();
        other.addInt("N_m", "other", 7);

        bool thrown = false;
        std::string msg;
        UnitTest::deactivate_output();
        try { T.update(other, gsOptionList::errorIfUnknown); }
        catch (const std::runtime_error & e) { thrown = true; msg = e.what(); }
        UnitTest::reactivate_output();

        CHECK(thrown);
        CHECK(msg.find("N_m")  != std::string::npos);
        CHECK(msg.find("real") != std::string::npos);
        CHECK(msg.find("int")  != std::string::npos);
    }

    /***
     * update() with errorIfUnknown: when an unknown label makes it throw, no
     * valid change has been applied
     */
    TEST(update_error_unknown_label_leaves_list_unchanged)
    {
        gsOptionList T = makeUpdateTarget();
        gsOptionList other = makeValidChanges();
        other.addInt("N_q", "other", 5);
        const gsOptionList before = T;

        bool thrown = false;
        UnitTest::deactivate_output();
        try { T.update(other, gsOptionList::errorIfUnknown); }
        catch (const std::runtime_error &) { thrown = true; }
        UnitTest::reactivate_output();

        CHECK(thrown);
        checkListsEqual(before, T);
    }

    /***
     * update() with errorIfUnknown: when a type mismatch makes it throw, no
     * valid change has been applied
     */
    TEST(update_error_type_mismatch_leaves_list_unchanged)
    {
        gsOptionList T = makeUpdateTarget();
        gsOptionList other = makeValidChanges();
        other.addInt("N_m", "other", 7);
        const gsOptionList before = T;

        bool thrown = false;
        UnitTest::deactivate_output();
        try { T.update(other, gsOptionList::errorIfUnknown); }
        catch (const std::runtime_error &) { thrown = true; }
        UnitTest::reactivate_output();

        CHECK(thrown);
        checkListsEqual(before, T);
    }

    /***
     * update() with errorIfUnknown: a bad entry in any of the four tables
     * (unknown string/real/switch label, or a real label given as int) throws
     * naming the label, and leaves the target unchanged although `other`
     * also carries valid changes in every table
     */
    TEST(update_error_rejects_bad_entries_in_every_table)
    {
        // Runs one case: the update must throw a runtime_error whose message
        // contains every word of `needles`, and must not modify the target.
        auto expectRejected = [](const gsOptionList & other,
                                 const std::vector<std::string> & needles)
        {
            gsOptionList T = makeUpdateTarget();
            const gsOptionList snapshot = T;

            UnitTest::deactivate_output();
            CHECK_THROW(T.update(other, gsOptionList::errorIfUnknown), std::runtime_error);

            std::string msg;
            try { T.update(other, gsOptionList::errorIfUnknown); }
            catch (const std::runtime_error & e) { msg = e.what(); }
            UnitTest::reactivate_output();

            for (size_t i = 0; i != needles.size(); ++i)
                CHECK(msg.find(needles[i]) != std::string::npos);
            checkListsEqual(T, snapshot);
        };

        {   // unknown string label
            gsOptionList other = makeValidChanges();
            other.addString("S_zz", "other", "x");
            expectRejected(other, {"S_zz"});
        }
        {   // unknown real label
            gsOptionList other = makeValidChanges();
            other.addReal("R_zz", "other", 1.5);
            expectRejected(other, {"R_zz"});
        }
        {   // unknown switch label
            gsOptionList other = makeValidChanges();
            other.addSwitch("B_zz", "other", true);
            expectRejected(other, {"B_zz"});
        }
        {   // R_a is a real in the target, an int in other (a label is unique
            // across the tables of one list, so other cannot start from
            // makeValidChanges(), which holds R_a as a real)
            gsOptionList other;
            other.addString("S_a", "other", "beta");
            other.addInt   ("N_a", "other", 10);
            other.addInt   ("R_a", "other", 7);
            other.addReal  ("N_m", "other", 3.5);
            other.addSwitch("B_a", "other", true);
            expectRejected(other, {"R_a", "real", "int"});
        }
    }

    /***
     * update() with errorIfUnknown: a fully valid list applies every value of
     * every type and leaves descriptions and size of the target untouched
     */
    TEST(update_error_valid_update_applies_all)
    {
        gsOptionList T = makeUpdateTarget();
        gsOptionList other = makeValidChanges();
        other.addReal("N_m", "other", 3.5);
        const gsOptionList before = T;

        // Every value of other differs from the target's
        CHECK(before.getString("S_a") != other.getString("S_a"));
        CHECK(before.getInt("N_a")    != other.getInt("N_a"));
        CHECK(before.getInt("N_z")    != other.getInt("N_z"));
        CHECK(before.getReal("N_m")   != other.getReal("N_m"));
        CHECK(before.getReal("R_a")   != other.getReal("R_a"));
        CHECK(before.getSwitch("B_a") != other.getSwitch("B_a"));

        T.update(other, gsOptionList::errorIfUnknown);

        CHECK_EQUAL("beta", T.getString("S_a"));
        CHECK_EQUAL(10, T.getInt("N_a"));
        CHECK_EQUAL(30, T.getInt("N_z"));
        CHECK_EQUAL(3.5,  T.getReal("N_m"));
        CHECK_EQUAL(0.75, T.getReal("R_a"));
        CHECK_EQUAL(true, T.getSwitch("B_a"));

        CHECK_EQUAL(before.size(), T.size());
        const std::vector<gsOptionList::OptionListEntry> eb = before.getAllEntries();
        const std::vector<gsOptionList::OptionListEntry> et = T.getAllEntries();
        CHECK_EQUAL(eb.size(), et.size());
        for (size_t i = 0; i < eb.size() && i < et.size(); ++i)
        {
            CHECK_EQUAL(eb[i].type,  et[i].type);
            CHECK_EQUAL(eb[i].label, et[i].label);
            CHECK_EQUAL(eb[i].desc,  et[i].desc);
        }
    }

    /***
     * update() from an empty list changes nothing, in every mode
     */
    TEST(update_empty_other_all_modes)
    {
        const gsOptionList other;
        const gsOptionList::updateType modes[3] =
            { gsOptionList::ignoreIfUnknown, gsOptionList::addIfUnknown,
              gsOptionList::errorIfUnknown };
        for (int i = 0; i < 3; ++i)
        {
            gsOptionList T = makeUpdateTarget();
            T.update(other, modes[i]);
            checkListsEqual(T, makeUpdateTarget());
        }
    }

    /***
     * update() into an empty list with ignoreIfUnknown (explicit and default)
     * adds nothing
     */
    TEST(update_into_empty_target_ignore)
    {
        const gsOptionList other = makeValidChanges();
        gsOptionList T;
        T.update(other, gsOptionList::ignoreIfUnknown);
        CHECK_EQUAL(size_t(0), T.size());
        T.update(other);
        CHECK_EQUAL(size_t(0), T.size());
    }

    /***
     * update() into an empty list with addIfUnknown copies the other list,
     * descriptions included
     */
    TEST(update_into_empty_target_add)
    {
        const gsOptionList other = makeValidChanges();
        gsOptionList T;
        T.update(other, gsOptionList::addIfUnknown);
        checkListsEqual(T, other);
    }

    /***
     * update() into an empty list with errorIfUnknown throws on the first
     * validated label (the only string of other) and leaves the target empty
     */
    TEST(update_into_empty_target_error)
    {
        const gsOptionList other = makeValidChanges();
        gsOptionList T;
        std::string msg;
        CHECK(updateThrows(T, other, gsOptionList::errorIfUnknown, msg));
        CHECK(msg.find("unknown option") != std::string::npos);
        CHECK(msg.find("S_a") != std::string::npos);
        CHECK_EQUAL(size_t(0), T.size());
    }

    /***
     * Updating a list from itself is safe and a no-op in every mode: all labels
     * of other exist in the target, so nothing is inserted while iterating
     */
    TEST(update_self_all_modes)
    {
        const gsOptionList::updateType modes[3] =
            { gsOptionList::ignoreIfUnknown, gsOptionList::addIfUnknown,
              gsOptionList::errorIfUnknown };
        for (int i = 0; i < 3; ++i)
        {
            gsOptionList T = makeUpdateTarget();
            T.update(T, modes[i]);
            checkListsEqual(T, makeUpdateTarget());
        }
    }

    /***
     * update() with grouped (dotted) labels in all three modes
     */
    TEST(update_group_labels_all_modes)
    {
        gsOptionList base;
        base.addInt   ("Name", "dName", 1);
        base.addSwitch("Flag", "dFlag", false);
        const gsOptionList target = base.wrapIntoGroup("Group");

        gsOptionList o;
        o.addInt ("Name",  "other", 5);
        o.addReal("Extra", "other", 2.5);
        const gsOptionList other = o.wrapIntoGroup("Group");
        CHECK(hasLabel(other, "Group.Extra"));

        {   // ignoreIfUnknown
            gsOptionList T = target;
            T.update(other, gsOptionList::ignoreIfUnknown);
            CHECK_EQUAL(5, T.getInt("Group.Name"));
            CHECK(!hasLabel(T, "Group.Extra"));
            CHECK_EQUAL(size_t(2), T.size());
        }
        {   // addIfUnknown
            gsOptionList T = target;
            T.update(other, gsOptionList::addIfUnknown);
            CHECK_EQUAL(5, T.getInt("Group.Name"));
            CHECK_EQUAL(2.5, T.getReal("Group.Extra"));
            CHECK_EQUAL(size_t(3), T.size());
            CHECK(T.hasGroup("Group"));
            if (T.hasGroup("Group"))
                CHECK_EQUAL(5, T.getGroup("Group").getInt("Name"));
        }
        {   // errorIfUnknown with an unknown dotted label
            gsOptionList T = target;
            std::string msg;
            CHECK(updateThrows(T, other, gsOptionList::errorIfUnknown, msg));
            CHECK(msg.find("Group.Extra") != std::string::npos);
            checkListsEqual(T, target);
        }
        {   // errorIfUnknown with a valid list applies
            gsOptionList v;
            v.addInt("Name", "other", 5);
            const gsOptionList valid = v.wrapIntoGroup("Group");
            gsOptionList T = target;
            T.update(valid, gsOptionList::errorIfUnknown);
            CHECK_EQUAL(5, T.getInt("Group.Name"));
        }
    }

    /***
     * update() without a mode behaves as ignoreIfUnknown: unknown labels are
     * neither added nor rejected
     */
    TEST(update_default_mode_is_ignore)
    {
        gsOptionList other = makeValidChanges();
        other.addInt("N_q", "other", 5);

        gsOptionList A = makeUpdateTarget();
        gsOptionList B = makeUpdateTarget();
        A.update(other);
        B.update(other, gsOptionList::ignoreIfUnknown);
        CHECK(!hasLabel(A, "N_q"));
        checkListsEqual(A, B);
    }

    /***
     * errorIfUnknown: a switch of the target given as a string in other is
     * rejected with a message naming the label and both types ("switch")
     */
    TEST(update_string_vs_switch_mismatch)
    {
        const gsOptionList target = makeUpdateTarget();
        gsOptionList other;
        other.addString("B_a", "other", "yes");
        other.addInt   ("N_a", "other", 10);

        std::string msg;
        {
            gsOptionList T = target;
            CHECK(updateThrows(T, other, gsOptionList::ignoreIfUnknown, msg));
        }
        gsOptionList T = target;
        CHECK(updateThrows(T, other, gsOptionList::errorIfUnknown, msg));
        CHECK(msg.find("B_a") != std::string::npos);
        CHECK(msg.find("has type switch, not string") != std::string::npos);
        CHECK(msg.find("bool") == std::string::npos);
        checkListsEqual(T, target);
    }

    /***
     * errorIfUnknown: a string of the target given as a switch in other is
     * rejected with a message naming the label and both types ("switch")
     */
    TEST(update_switch_vs_string_mismatch)
    {
        const gsOptionList target = makeUpdateTarget();
        gsOptionList other;
        other.addSwitch("S_a", "other", true);
        other.addInt   ("N_a", "other", 10);

        std::string msg;
        {
            gsOptionList T = target;
            CHECK(updateThrows(T, other, gsOptionList::ignoreIfUnknown, msg));
        }
        gsOptionList T = target;
        CHECK(updateThrows(T, other, gsOptionList::errorIfUnknown, msg));
        CHECK(msg.find("S_a") != std::string::npos);
        CHECK(msg.find("has type string, not switch") != std::string::npos);
        CHECK(msg.find("bool") == std::string::npos);
        checkListsEqual(T, target);
    }

}
