

QList<int> getDegrees(QDomElement e)
{
    QList<int> deg ;
    QDomNodeList nodelist = e.elementsByTagName("order") ;
    for(unsigned int i = 0 ; i < nodelist.length() ; i++)  {
        QDomElement element = nodelist.item(i).toElement() ;
        QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
        foreach(QString s, list) deg << s.toInt() - 1;
    }
    return deg;
}


void getKnots(QDomElement e, int i, gsKnotVector<> & result);
{ 
    // Note: degree NOT set inside this fuction
    QDomNodeList nodelist = e.elementsByTagName("knots") ;
    QDomElement element = nodelist.item(i).toElement() ;
    QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
    int knotsSize = list.size();

    result.clear();
    result.reserve(knotsSize);
    foreach(QString s, list)
    {
        result.push_back( s.toDouble() );
    }
}

void getControlPoints(QDomElement e, int i, int dim, gsMatrix<> & result)
{
    QDomNodeList nodelist = e.elementsByTagName("points") ;
    QDomElement element = nodelist.item(i).toElement() ;
    QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
    int coeffsSize = list.size();
    result.resize(dim, coeffsSize / dim);
    double * coeffs = result.data();
    for(QStringList::iterator itr = list.begin() ; itr != list.end() ; itr++) 
    {
        *coeffs = itr->toFloat();
        ++coeffs;
    }

    result.transposeInPlace();
}
