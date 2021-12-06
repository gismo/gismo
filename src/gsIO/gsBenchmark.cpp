/** @file gsBenchmark.cpp

    @brief Provides implemementation of generic benchmarking framework.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#include <gsIO/gsBenchmark.h>
#include <gsCore/gsJITCompiler.h>
#include <gsCore/gsSysInfo.h>
#include <cstring>

namespace gismo
{

  std::ostream &gsBenchmark::gsBenchmarkResultSet::print(std::ostream &os) const
  {
    os << "\\pgfplotstableread[row sep=\\\\,col sep=&]{\n"
       << "threads & " << label << " \\\\\n";

    for (auto it=results.cbegin(); it!=results.cend(); ++it)
        os << it->at(0) << "&" << it->at(2) << "\\\\\n";

    os << "}\\data" << label << "\n";

    return os;
  }

  std::ostream &gsBenchmark::gsBenchmarkSet::print(std::ostream &os) const
  {
    for (auto it=results.cbegin(); it!=results.cend(); ++it)
      (*it)->print(os);

    os << "\\begin{tikzpicture}\n"
       << "\\begin{semilogyaxis}[\n"
       << "name=MyAxis,\n"
       << "width=\\textwidth,\n"
       << "height=.5\\textwidth,\n"
       << "legend pos=outer north east,\n"        
       << "symbolic x coords={";

    //std::vector<std::array<double, 4> >::const_iterator
    auto it  = results.front()->get().cbegin();
    auto ite = results.front()->get().cend();
    for (;it!=ite; ++it)
        os << it->at(0) << (it!=ite-1 ? "," : "");
    os << "},\n"
       << "xlabel={OpenMP threads},\n";

    it = results.front()->get().cbegin();
    switch( (int)it->at(3) )
    {
    case metric::bandwidth_kb_sec:        
      os << "ylabel={Bandwidth in KB/s},\n";
      break;
    case metric::bandwidth_mb_sec:        
      os << "ylabel={Bandwidth in MB/s},\n";
      break;
    case metric::bandwidth_gb_sec:        
      os << "ylabel={Bandwidth in GB/s},\n";
      break;
    case metric::bandwidth_tb_sec:        
      os << "ylabel={Bandwidth in TB/s},\n";
      break;
    case metric::perf_kflop_sec:        
      os << "ylabel={Berformance in kFLOP/s},\n";
      break;
    case metric::perf_mflop_sec:        
      os << "ylabel={Berformance in mFLOP/s},\n";
      break;
    case metric::perf_gflop_sec:        
      os << "ylabel={Berformance in gFLOP/s},\n";
      break;
    case metric::perf_tflop_sec:        
      os << "ylabel={Berformance in tFLOP/s},\n";
      break;
    case metric::runtime_sec:
      os << "ylabel={Runtime in seconds},\n";
      break;
    default:
      throw std::runtime_error("Unsupported metric");
    }
      
    os << "title={" << title << "},\n"
       << "]";

    for (auto rit=results.cbegin(); rit!=results.cend(); ++rit)
      os << "\\addplot table[x=threads,y="
         << (*rit)->get_label()
         << "]{\\data"
         << (*rit)->get_label()
         << "};\n";

    os << "\\legend{";
    for (auto rit=results.cbegin(); rit!=results.cend(); ++rit)
      os << (*rit)->get_title() << (rit!=results.cend()-1 ? "," : "");
    os << "}\n"
        
       << "\\end{semilogyaxis}\n"

       << "\\path let \\p1=(MyAxis.west), \\p2=(MyAxis.east) in "
       << "node[below right, align=left, text=black, text width=\\x2-\\x1]\n"
       << "at ($(MyAxis.south west)+(0,-30pt)$) {%\n"
       << "G+Smo " << gsSysInfo::getGismoVersion()
       << ", Eigen " << gsSysInfo::getEigenVersion()
       << " (" << gsSysInfo::getCompilerVersion()
       << ", " << gsSysInfo::getCppVersion()
       << ", " << gsSysInfo::getStdLibVersion()
       << (gsSysInfo::getExtraLibsVersion().empty()
           ? "), \n"
           : gsSysInfo::getExtraLibsVersion()+"), \n")

       << "CPU " << gsSysInfo::getCpuInfo() << ", "
       << "Memory " << gsSysInfo::getMemoryInfo()  << ", ";

    gsJITCompilerConfig jit; jit.load(GISMO_CONFIG_DIR "jit.xml");
    std::string flags = jit.getFlags();
    os << "Compiler flags ";
      
    for (auto token=strtok(&flags[0], " "); token!=NULL; token=strtok(NULL, " ")) {
      if (token[0]=='-') {
        if (token[1]=='I' || token[1]=='L' || token[1]=='l' || token[1]=='W')
          continue;
        os << "\\verb!" << token << "! ";
      }
    }
      
    os << "};\n"
       << "\\end{tikzpicture}\n";

    return os;
  }

  std::ostream &gsBenchmark::print(std::ostream &os) const
  {
    os << "\\documentclass[tikz]{standalone}\n"
       << "\\usepackage{pgfplots}\n"
       << "\\usepackage{verbatim}\n"
       << "\\begin{document}\n"
       << "\\usetikzlibrary{calc}\n";
    
    for (auto it=benchmarks.cbegin(); it!=benchmarks.cend(); ++it)
      (*it)->print(os);
    
    os << "\\end{document}\n";
    return os;
  }
  
} // namespace gismo
