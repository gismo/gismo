#include <gismo.h>
#include <gsModeling/gsSurfaceReparameterization.h>

using namespace gismo;

int main(int argc, char *argv[]) {
  // Handle command-line arguments for the input file
  std::string INPUT_FILE = (argc > 1) ? argv[1] : "surfaces/crazySurf.xml";

  // Load XML file - multi-patch domain
  if (!gsFileManager::fileExists(INPUT_FILE)) {
	gsWarn << "The file cannot be found!\n";
	return EXIT_FAILURE;
  }

  // MultiPatch reader
  gsInfo << "Reading file: " << INPUT_FILE << "\n";
  gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(INPUT_FILE);
  gsInfo << "Loaded geometry: " << *mp << "\n";

  gsWriteParaview(*mp, "input_surface", 1000);

  // Create the surface reparametrization object
  SurfaceReparameterization<real_t> reparam(*mp);

  // Generate the final optimized geometry as a B-Spline surface
  gsMultiPatch<real_t> optSurface = reparam.solve();

  // Output the resulting geometry to a Paraview file
  gsWriteParaview(optSurface, "optimized_surface", 1000);

  return EXIT_SUCCESS;
}
