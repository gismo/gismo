#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsHSplinesParaview.h>
#include <gsHSplines/gsHSplinesParaview.hpp>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{


TEMPLATE_INST
void gsWriteParaview(const gsHBox<2,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void gsWriteParaview(const gsHBox<3,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void gsWriteParaview(const gsHBoxContainer<2,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void gsWriteParaview(const gsHBoxContainer<3,T> & hbox, std::string const & fn, short_t mode);

TEMPLATE_INST
void writeSingleHBox(const gsHBox<2,T> & box, std::string const & fn);

TEMPLATE_INST
void writeSingleHBox(const gsHBox<3,T> & box, std::string const & fn);

} // namespace gismo

#undef T
#undef uZ
#undef Z
