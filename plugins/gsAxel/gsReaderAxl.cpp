
#include "gsReaderAxl.h"
#include "gsGeometryData.h"
#include "gsMultiPatchData.h"
#include "gsBasisData.h"

#include <dtkCoreSupport/dtkAbstractData.h>
#include <dtkCoreSupport/dtkAbstractDataFactory.h>

//to read fields applied on this object
#include <axlCore/axlAbstractField.h>
#include <axlCore/axlAbstractFieldParametricSurface.h>
#include <axlCore/axlFieldReadersFactory.h>

namespace {

int getDimension(QDomElement e)
{
    int d = 0 ;
    QDomNodeList nodelist = e.elementsByTagName("dimension") ;
    for(int i = 0 ; i < nodelist.length() ; i++)  {
        QDomElement element = nodelist.item(i).toElement() ;
        d = element.text().simplified().toInt() ;
    }
    return d ;
}

QList<int> getNumbers(QDomElement e)
{
    QList<int> numbers ;
    QDomNodeList nodelist = e.elementsByTagName("number") ;
    for(int i = 0 ; i < nodelist.length() ; i++)  {
        QDomElement element = nodelist.item(i).toElement() ;
        QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
        foreach(QString s, list) numbers << s.toInt() ;
    }
    return numbers ;
}

QList<int> getOrders(QDomElement e)
{
    QList<int> orders ;
    QDomNodeList nodelist = e.elementsByTagName("order") ;
    for(int i = 0 ; i < nodelist.length() ; i++)  {
        QDomElement element = nodelist.item(i).toElement() ;
        QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
        foreach(QString s, list) orders << s.toInt() ;
    }
    return orders ;
}

double *getKnots(QDomElement e, int i)
{
    QDomNodeList nodelist = e.elementsByTagName("knots") ;
    QDomElement element = nodelist.item(i).toElement() ;
    QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
    int knotsSize = list.size();
    double *knots = new double[knotsSize] ;
    int j=0;
    foreach(QString s, list)
    {
        knots[j]= s.toDouble() ;
        j++;
    }
    return knots ;
}

double *getControlPoints(QDomElement e,int i, int dim)
{
    QDomNodeList nodelist = e.elementsByTagName("points") ;
    QDomElement element = nodelist.item(i).toElement() ;
    QStringList list = element.text().simplified().split(QRegExp("\\s+")) ;
    int coeffsSize = list.size();
    //qDebug()<<"-- nb coef "<<coeffsSize;
    //qDebug()<<list;
    double *coeffs = new double[coeffsSize] ;
    int j=0;
    for(QStringList::iterator itr = list.begin() ; itr != list.end() ; itr+=1) {
        coeffs[j]=itr->toFloat();
        j++;
        /*
        if(dim > 1){
            coeffs[j+1]=(itr+1)->toFloat();
        }
        if(dim > 2)
        {
            coeffs[j+2]=(itr+2)->toFloat();
            j++;
            itr++;
        }
        if(dim> 1){
            j+=2;
            itr+=1;
        }else{
            j++;
        } */
    }
   
    return coeffs ;
}

}


gsReaderAxl::gsReaderAxl(void)
{

}

gsReaderAxl::~gsReaderAxl(void)
{

}

QString gsReaderAxl::identifier(void) const
{
    return "gsReaderAxl";
}

QString gsReaderAxl::description(void) const
{
    return "gsReaderAxl";
}

QStringList gsReaderAxl::handled(void) const
{
    return QStringList() << "SplineSurface";
}

bool gsReaderAxl::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataReaderType("gsReaderAxl", QStringList(), creategsReaderAxl);
}

bool gsReaderAxl::accept(const QDomNode& node)
{
    QDomElement element = node.toElement();

    if(element.tagName() != "surface")
    // or   element.tagName() != "curve" 
        return false;

    if(element.attribute("type") != "bspline")
        return false;

    if(!hasChildNode(element, "dimension"))
        return false;

    if(!hasChildNode(element, "number"))
        return false;

    if(!hasChildNode(element, "order"))
        return false;

    if(!hasChildNode(element, "knots"))
        return false;

    if(!hasChildNode(element, "points"))
        return false;

    return true;
}

bool gsReaderAxl::reject(const QDomNode& node)
{
    return !this->accept(node);
}

axlAbstractData *gsReaderAxl::read(const QDomNode& node)
{
    QDomElement element = node.toElement();

    axlAbstractSurfaceBSpline *bsplineSurface = dynamic_cast<axlAbstractSurfaceBSpline *>(gsAxelPlugin::dataFactSingleton->create("SplineSurface"));
    //if there are some field, read them thanks to the factory.
    QDomNodeList nodeListField = element.elementsByTagName("field");
    if(!nodeListField.isEmpty()){
        for(int i =0; i < nodeListField.count(); i++){
            QDomElement fieldElement = nodeListField.at(i).toElement();
            QString fieldType = fieldElement.attribute("type");
            if(!fieldType.isEmpty()){
                axlAbstractDataReader *field_reader = dynamic_cast<axlAbstractDataReader *>(axlFieldReadersFactory::instance()->create(fieldType));
                axlAbstractField * fieldToAdd = dynamic_cast<axlAbstractField *>(field_reader->read(fieldElement));
                if(fieldToAdd){
                    dtkWarn()<< fieldToAdd->identifier();
                    bsplineSurface->addField(fieldToAdd);
                    if(axlAbstractFieldParametricSurface *paramField = dynamic_cast<axlAbstractFieldParametricSurface *>(fieldToAdd)){
                        paramField->setSurface(bsplineSurface);
                    }
                }
            }
        }

        for(int j = 0;j < nodeListField.count(); j++){
            QDomNode child = nodeListField.at(j);
            element.removeChild(child);
        }
    }

    QString r = element.attribute("rational");
    bool rational = (!r.isEmpty()) ;

    int dimension = getDimension(element) ;
    QList<int> orders = getOrders(element) ;
    QList<int> numbers = getNumbers(element) ;
    double *knots_u = getKnots(element, 0) ;
    double *knots_v = getKnots(element, 1) ;
    double *points = getControlPoints(element, 0, dimension+rational) ;

    //qDebug()<<"test";
    //axlAbstractSurfaceBSpline *bsplineSurface = dynamic_cast<axlAbstractSurfaceBSpline *>(gsAxelPlugin::dataFactSingleton->create("SplineSurface"));

    bsplineSurface->setSurface(numbers[0], numbers[1], orders[0], orders[1], dimension, knots_u, knots_v, points, rational);

    delete knots_u;
    delete knots_v;
    delete points;

    QString name = element.attribute("name");
    if(!name.isEmpty())
    {
        bsplineSurface->setObjectName(name);
    }

    QString color = element.attribute("color");
    if(!color.isEmpty())
    {
        QStringList colorList = color.split(" ");
        if(colorList.size() == 3) // rgb components
            bsplineSurface->setColor(QColor(colorList.at(0).toInt(), colorList.at(1).toInt(), colorList.at(2).toInt()));

    }

    QString size = element.attribute("size");
    if(!size.isEmpty())
    {
        bsplineSurface->setSize(size.toDouble());
    }

    QString shader = element.attribute("shader");
    QString dirShader;
    if(!shader.isEmpty())
    {
        // try to read from axelShader.qrc
        dirShader = ":axlShader/shader/"+shader;
        if(!QFile::exists(dirShader))
        {// if shader wasn't found, we try to look at the shader setting path
            QSettings settings("inria", "dtk");
            QString defaultPath;
            settings.beginGroup("shader");
            dirShader = settings.value("path", defaultPath).toString();
            qDebug()<< "gsReaderAxl:: dirShader" <<dirShader;

            //            dirShader = this->file().left(this->file().lastIndexOf("axel-data") + 9);
            dirShader.append("/Shader/"+shader);
        }
        bsplineSurface->setShader(dirShader);
    }

    return bsplineSurface;
}

dtkAbstractDataReader *creategsReaderAxl(void)
{
    return new gsReaderAxl;
}
