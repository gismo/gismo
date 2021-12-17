/** @file gsBenchmark.cpp

    @brief Provides implemementation of generic benchmarking framework.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Moller
*/

#include <gsCore/gsJITCompiler.h>
#include <gsCore/gsSysInfo.h>
#include <gsIO/gsBenchmark.h>

#include <cstring>

namespace gismo
{

  std::ostream &gsBenchmark::gsBenchmarkResultSet::print(std::ostream &os) const
  {
    os << "\\pgfplotstableread[col sep=space]{\n"
       << label << "\n";

    for (auto it=results.cbegin(); it!=results.cend(); ++it)
        os << it->at(2) << "\n";

    os << "}\\data" << label << "\n";

    return os;
  }

  std::ostream &gsBenchmark::gsBenchmarkSet::print(std::ostream &os) const
  {
    for (auto it=results.cbegin(); it!=results.cend(); ++it)
      (*it)->print(os);

    os << "\\begin{tikzpicture}\n"
       << "\\begin{axis}[\n"
       << "name=MyAxis,\n"
       << "width=2\\textwidth,\n"
       << "height=.8\\textwidth,\n"
       << "legend pos=outer north east,\n"
       << "ybar = 0.05cm,\n"
       << "bar width = 3pt,\n"
       << "ymajorgrids=true,\n"
       << "xticklabel style={rotate=45,anchor=east},\n"
       << "xticklabels={";

    for (auto rit=results.cbegin(); rit!=results.cend(); ++rit)
      os << (*rit)->get_title() << (rit!=results.cend()-1 ? "," : "");

    os << "},\n"
       << "xtick=data,\n";

    auto it = results.front()->get().cbegin();
    if ((metric)it->at(3) & metric::speedup) {
      switch( (int)it->at(3) & ~metric::speedup ) {
      case metric::bandwidth_kb_sec:
      case metric::bandwidth_mb_sec:
      case metric::bandwidth_gb_sec:
      case metric::bandwidth_tb_sec:
        os << "ylabel={Bandwidth [speedup]},\n";
        break;
      case metric::perf_kflop_sec:
      case metric::perf_mflop_sec:
      case metric::perf_gflop_sec:
      case metric::perf_tflop_sec:
        os << "ylabel={Performance [speedup]},\n";
        break;
      case metric::runtime_sec:
        os << "ylabel={Runtime [speedup]},\n";
        break;
      default:
        GISMO_ERROR("Unsupported metric");
      }
    } else {
      switch( (int)it->at(3) & ~metric::speedup ) {
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
        os << "ylabel={Performance in kFLOP/s},\n";
        break;
      case metric::perf_mflop_sec:
        os << "ylabel={Performance in mFLOP/s},\n";
        break;
      case metric::perf_gflop_sec:
        os << "ylabel={Performance in gFLOP/s},\n";
        break;
      case metric::perf_tflop_sec:
        os << "ylabel={Performance in tFLOP/s},\n";
        break;
      case metric::runtime_sec:
        os << "ylabel={Runtime in seconds},\n";
        break;
      default:
        GISMO_ERROR("Unsupported metric");
      }
    }

    os << "title={" << title
       << " [real\\_t:" << util::type<real_t>::name()
           << ", index\\_t:" << util::type<index_t>::name()
           << ", short\\_t:" << util::type<short_t>::name()<< "]},\n"
       << "]\n";

    for (auto rit=results.cbegin()+1; rit!=results.cend(); ++rit)
      os << "\\pgfplotstablecreatecol[copy column from "
         << "table={\\data"
         << (*rit)->get_label()
         << "}{[index] 0}] {"
         << (*rit)->get_label()
         << "} {\\data"
         << (*results.cbegin())->get_label()
         << "}\n";

    os << "\\pgfplotstabletranspose[rows/threads/.style={string type}]\\mytable{"
       << "\\data"
       << (*results.cbegin())->get_label()
       << "}\n";

    for (std::size_t i=1; i<=results.front()->get().size(); ++i)
      os << "\\addplot table[x expr=\\coordindex, y index="
         << util::to_string(i) << "]{\\mytable};\n";

    os << "\\legend{";
    it  = results.front()->get().cbegin();
    auto ite = results.front()->get().cend();
    for (;it!=ite; ++it)
      os << "Threads=" << it->at(0) << (it!=ite-1 ? "," : "");
    os << "}\n"

       << "\\end{axis}\n"

       << "\\gettikzxy{(MyAxis.south west)}{\\ax}{\\ay}\n"
       << "\\gettikzxy{(MyAxis.outer south east)}{\\bx}{\\by}\n"

       << "\\path let \\p1=(MyAxis.west), \\p2=(MyAxis.east) in "
       << "node[draw,below right, align=left, text=black, text width=\\x2-\\x1-10pt, minimum width=\\x2-\\x1]\n"
       << "at ($(\\ax, \\by-10pt)$) {%\n"
       << "G+Smo " << gsSysInfo::getGismoVersion()
       << ", Eigen " << gsSysInfo::getEigenVersion()
       << " (" << gsSysInfo::getCompilerVersion()
       << ", " << gsSysInfo::getCppVersion()
       << ", " << gsSysInfo::getStdLibVersion()
       << (gsSysInfo::getExtraLibsVersion().empty()
           ? "), \n"
           : ", "+gsSysInfo::getExtraLibsVersion()+"), \n")

       << "CPU " << gsSysInfo::getCpuInfo() << ", "
       << "Memory " << gsSysInfo::getMemoryInfo()  << "\\\\\n";

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
       << "\\usepackage{pgfplotstable}\n"
       << "\\usepackage{verbatim}\n"
       << "\\pgfplotsset{compat=1.18}\n"
       << "\\makeatletter\n"
       << "\\newcommand{\\gettikzxy}[3]{%\n"
       << "\\tikz@scan@one@point\\pgfutil@firstofone#1\\relax\n"
       << "\\edef#2{\\the\\pgf@x}%\n"
       << "\\edef#3{\\the\\pgf@y}%\n"
       << "}\n"
       << "\\makeatother\n"
       << "\\begin{document}\n"
       << "\\usetikzlibrary{calc}\n";

    for (auto it=benchmarks.cbegin(); it!=benchmarks.cend(); ++it)
      (*it)->print(os);

    os << "\\end{document}\n";
    return os;
  }

} // namespace gismo
