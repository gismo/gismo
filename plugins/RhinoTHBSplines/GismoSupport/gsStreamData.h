/** @file gsStreamData.h
    @brief Provides stream interface for reading and writing XML data

*/
#pragma once

/** @brief stream xml data
*/
class gsStreamData :
    public gsFileData<double>
{
    // Note: the reason this class exists is because gsFileData<>::readGismoXmlStream is protected
    //       and not accessible by users of gsFileData<>. To circumvent this shortcoming this class
    //       offers reading and writing of XML to an istream and ostream, respectively.
public:
    gsStreamData();
    ~gsStreamData();

    /// Read XML from a stream
    bool ReadXmlStream(std::istream&);

    /// Write XML to a stream
    bool WriteXmlStream(std::ostream&);

};

