#include <iostream>
#include <string>
#include <gismo.h>
#include <autodiff/reverse/var.hpp>

using namespace gismo;

// forward declaration of some utility functions
void print(const gsBSplineBasis<autodiff::var>& bsb, const std::string& name);

int main(int argc, char* argv[])
{
    gsInfo << "--- Test for autodiff::var in gsBSplineBasis ---\\n";

    autodiff::var a = 0;
    autodiff::var b = 1;
    index_t interior = 4;
    index_t multEnd = 3;
    int degree = multEnd - 1;

    gsKnotVector<autodiff::var> kv(a, b, interior, multEnd);
    gsBSplineBasis<autodiff::var> bsb(kv);
    print(bsb, "bsb");

    return 0;
}

void print(const gsBSplineBasis<autodiff::var>& bsb,
           const std::string& name)
{
    gsInfo << name << ": \\n";
    // bsb.print(gsInfo); // print is not supported for autodiff::var
    gsInfo << "Size: " << bsb.size() << std::endl;
    gsInfo << "Degree: " << bsb.degree() << std::endl;
    gsInfo << "\\n\\n";
}