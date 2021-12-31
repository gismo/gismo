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
  std::ostream &gsBenchmarkResultSet::to_tikz(std::ostream &os) const
  {
    os << "\\pgfplotstableread[col sep=space]{\n"
       << label << "\n";

    for (auto it=results.cbegin(); it!=results.cend(); ++it)
        os << it->value << "\n";

    os << "}\\data" << label << "\n";

    return os;
  }

  std::ostream &gsBenchmarkSet::to_tikz(std::ostream &os) const
  {
    for (auto it=results.cbegin(); it!=results.cend(); ++it)
      (*it)->to_tikz(os);

    os << "\\begin{tikzpicture}\n"
       << "\\begin{axis}[\n"
       << "name=MyAxis,\n"
       << "width=2\\textwidth,\n"
       << "height=.8\\textwidth,\n"
       << "legend pos=outer north east,\n"
       << "ybar=0.05cm,\n"
       << "bar width=3pt,\n"
       << "ymajorgrids=true,\n"
       << "xticklabel style={rotate=45,anchor=east},\n"
       << "xticklabels={";

    for (auto rit=results.cbegin(); rit!=results.cend(); ++rit)
      os << (*rit)->get_descr() << (rit!=results.cend()-1 ? "," : "");

    os << "},\n"
       << "xtick=data,\n";

    auto it = results.front()->get().cbegin();
    if (it->metric & gismo::metric::speedup) {
      switch(it->metric & ~gismo::metric::speedup) {
      case gismo::metric::bandwidth_kb_sec:
      case gismo::metric::bandwidth_mb_sec:
      case gismo::metric::bandwidth_gb_sec:
      case gismo::metric::bandwidth_tb_sec:
        os << "ylabel={Bandwidth [speedup]},\n";
        break;
      case gismo::metric::perf_kflop_sec:
      case gismo::metric::perf_mflop_sec:
      case gismo::metric::perf_gflop_sec:
      case gismo::metric::perf_tflop_sec:
        os << "ylabel={Performance [speedup]},\n";
        break;
      case gismo::metric::runtime_sec:
        os << "ylabel={Runtime [speedup]},\n";
        break;
      default:
        GISMO_ERROR("Unsupported metric");
      }
    } else {
      switch(it->metric & ~gismo::metric::speedup) {
      case gismo::metric::bandwidth_kb_sec:
        os << "ylabel={Bandwidth in KB/s},\n";
        break;
      case gismo::metric::bandwidth_mb_sec:
        os << "ylabel={Bandwidth in MB/s},\n";
        break;
      case gismo::metric::bandwidth_gb_sec:
        os << "ylabel={Bandwidth in GB/s},\n";
        break;
      case gismo::metric::bandwidth_tb_sec:
        os << "ylabel={Bandwidth in TB/s},\n";
        break;
      case gismo::metric::perf_kflop_sec:
        os << "ylabel={Performance in kFLOP/s},\n";
        break;
      case gismo::metric::perf_mflop_sec:
        os << "ylabel={Performance in mFLOP/s},\n";
        break;
      case gismo::metric::perf_gflop_sec:
        os << "ylabel={Performance in gFLOP/s},\n";
        break;
      case gismo::metric::perf_tflop_sec:
        os << "ylabel={Performance in tFLOP/s},\n";
        break;
      case gismo::metric::runtime_sec:
        os << "ylabel={Runtime in seconds},\n";
        break;
      default:
        GISMO_ERROR("Unsupported metric");
      }
    }

    os << "title={" << descr
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
      os << "Threads=" << it->threads << (it!=ite-1 ? "," : "");
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

  std::ostream &gsBenchmark::to_tikz(std::ostream &os) const
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
      (*it)->to_tikz(os);

    os << "\\end{document}\n";
    return os;
  }

  std::ostream &gsBenchmarkResultSet::print(std::ostream &os) const
  {
    os << std::setw(8) << descr << " | ";
    for (auto it=results.cbegin(); it!=results.cend(); ++it)
      os << std::setw(4) << it->threads << " : "
         << std::setw(6) << std::scientific << std::setprecision(2) << it->value;
    os << "\n";
    return os;
  }

  std::ostream &gsBenchmarkSet::print(std::ostream &os) const
  {
    os << "[" << label << "] " << descr << "\n"
       << std::setw(8) << "memsize"
       << " | "
       << util::to_string(results.front()->get().size())
       << "x (#Threads : ";

    if (results.front()->get().cbegin()->metric & gismo::metric::speedup) {
      switch(results.front()->get().cbegin()->metric & ~gismo::metric::speedup) {
      case gismo::metric::bandwidth_kb_sec:
      case gismo::metric::bandwidth_mb_sec:
      case gismo::metric::bandwidth_gb_sec:
      case gismo::metric::bandwidth_tb_sec:
        os << "Bandwidth [speedup])\n";
        break;
      case gismo::metric::perf_kflop_sec:
      case gismo::metric::perf_mflop_sec:
      case gismo::metric::perf_gflop_sec:
      case gismo::metric::perf_tflop_sec:
        os << "Performance [speedup])\n";
        break;
      case gismo::metric::runtime_sec:
        os << "Runtime [speedup])\n";
        break;
      default:
        GISMO_ERROR("Unsupported metric");
      }
    } else {
      switch(results.front()->get().cbegin()->metric & ~gismo::metric::speedup) {
      case gismo::metric::bandwidth_kb_sec:
        os << "Bandwidth in KB/s)\n";
        break;
      case gismo::metric::bandwidth_mb_sec:
        os << "Bandwidth in MB/s)\n";
        break;
      case gismo::metric::bandwidth_gb_sec:
        os << "Bandwidth in GB/s)\n";
        break;
      case gismo::metric::bandwidth_tb_sec:
        os << "Bandwidth in TB/s)\n";
        break;
      case gismo::metric::perf_kflop_sec:
        os << "Performance in kFLOP/s)\n";
        break;
      case gismo::metric::perf_mflop_sec:
        os << "Performance in mFLOP/s)\n";
        break;
      case gismo::metric::perf_gflop_sec:
        os << "Performance in gFLOP/s)\n";
        break;
      case gismo::metric::perf_tflop_sec:
        os << "Performance in tFLOP/s)\n";
        break;
      case gismo::metric::runtime_sec:
        os << "Runtime in seconds)\n";
        break;
      default:
        GISMO_ERROR("Unsupported metric");
      }
    }
    
    for (auto it=results.cbegin(); it!=results.cend(); ++it)
      (*it)->print(os);
    return os;
  }

  std::ostream &gsBenchmark::print(std::ostream &os) const
  {
    for (auto it=benchmarks.cbegin(); it!=benchmarks.cend(); ++it)
      (*it)->print(os);
    return os;
  }
  
} // namespace gismo
