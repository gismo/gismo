/* gsNativeAxelReader.cpp ---
 *
 */

#include "gsAxelPlugin.h"

#include "gsNativeAxelReader.h"

#include "gsGeometryData.h"
#include "gsBasisData.h"
#include "gsMultiPatchData.h"

#include <axlCore/axlAbstractField.h>

#include <dtkCoreSupport/dtkAbstractData.h>
#include <dtkCoreSupport/dtkAbstractDataFactory.h>

//to read fields applied on thi object
#include <axlCore/axlAbstractField.h>
#include <axlCore/axlAbstractFieldParametricSurface.h>
#include <axlCore/axlFieldReadersFactory.h>

// /////////////////////////////////////////////////////////////////
// gsNativeAxelReader
// /////////////////////////////////////////////////////////////////

gsNativeAxelReader::gsNativeAxelReader(void)
{

}

gsNativeAxelReader::~gsNativeAxelReader(void)
{

}

QString gsNativeAxelReader::identifier(void) const
{
    return "gsNativeAxelReader";
}

QString gsNativeAxelReader::description(void) const
{
    return "gsNativeAxelReader";
}

QStringList gsNativeAxelReader::handled(void) const
{
    return QStringList() << "goSurfaceBSpline";
}

bool gsNativeAxelReader::registered(void)
{
    return gsAxelPlugin::dataFactSingleton->registerDataReaderType("goSurfaceBSplineReader", QStringList(), creategsNativeAxelReader);
}

bool gsNativeAxelReader::accept(const QDomNode& node)
{
    QDomElement element = node.toElement();

    if(element.tagName() == "surface")
      if(element.attribute("type") == "bspline")
        return true;

    if(element.tagName() == "curve")
      if(element.attribute("type") == "bspline")
        return true;

    return false;
}

bool gsNativeAxelReader::reject(const QDomNode& node)
{
    return !this->accept(node);
}

axlAbstractData *gsNativeAxelReader::read(const QDomNode& node)
{

  return NULL;
}

dtkAbstractDataReader *creategsNativeAxelReader(void)
{
    return new gsNativeAxelReader;
}
