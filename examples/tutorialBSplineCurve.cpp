#include <iostream>

#include <gismo.h>


using namespace gismo;
using std::cout;

int main(int argc, char *argv[])
{
    bool plot; // If set to true, paraview file is generated and launched on exit
    try 
    {
	gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
	gsArgSwitch ap("", "plot", "Plot result in ParaView format", cmd);
	cmd.parse(argc,argv);
	plot       = ap.getValue();
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << "\n"; return -1; }


  // Make a BSpline curve
  gsKnotVector<> kv(0,1,1,3);//start,end,interior knots, start/end multiplicites of knots1
  gsMatrix<> coefs(4,3);
  coefs << 0,0,0,  1,2,3, 2,1,4, 4,4,4 ;
  
  gsBSpline<> curve( kv, give(coefs));
  
  // Print the Bspline curve
  cout<< "I am a "<< curve <<"\n";

  if (plot) 
  {
      // Output a paraview file
      gsWriteParaview( curve, "bsplinecurve", 100);
      return system("bsplinecurve.pvd");
  }

  return 0;
}
