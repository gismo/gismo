#include "stdafx.h"
#include "gsStreamData.h"


gsStreamData::gsStreamData()
{
}


gsStreamData::~gsStreamData()
{
}

bool gsStreamData::ReadXmlStream(std::istream& is)
{
    return this->readGismoXmlStream(is);
}

bool gsStreamData::WriteXmlStream(std::ostream& os)
{
    os << *this;
    return true;
}
