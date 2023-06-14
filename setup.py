import os
import re
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

import pybind11

# taken from pybind cmake example
# https://github.com/pybind/cmake_example/blob/master/setup.py

# Convert distutils Windows platform specifiers to CMake -A arguments
PLAT_TO_CMAKE = {
    "win32": "Win32",
    "win-amd64": "x64",
    "win-arm32": "ARM",
    "win-arm64": "ARM64",
}


# A CMakeExtension needs a sourcedir instead of a file list.
# The name must be the _single_ output extension from the CMake build.
# If you need multiple extensions, see scikit-build.
class CMakeExtension(Extension):
    def __init__(
        self, name: str, sourcedir: str = "", extra_args: dict = None
    ) -> None:
        super().__init__(name, sources=[])
        self.sourcedir = os.fspath(Path(sourcedir).resolve())
        if extra_args is not None:
            self.extra_args = extra_args
            self.debug = extra_args.get("--debug", False)


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        # Must be in this form due to bug in .resolve()
        # only fixed in Python 3.10+
        ext_fullpath = Path.cwd() / self.get_ext_fullpath(
            ext.name
        )  # type: ignore[no-untyped-call]
        extdir = ext_fullpath.parent.resolve()

        # Using this requires trailing slash for auto-detection & inclusion of
        # auxiliary "native" libs

        debug = (
            int(os.environ.get("DEBUG", 0))
            if self.debug is None
            else self.debug
        )
        # overwrite if ext.debug exists
        debug = ext.debug if hasattr(ext, "debug") else debug
        cfg = "Debug" if debug else "Release"

        # CMake lets you override the generator - we need to check this.
        # Can be set with Conda-Build, for example.
        cmake_generator = os.environ.get("CMAKE_GENERATOR", "")

        # Set Python_EXECUTABLE instead if you use PYBIND11_FINDPYTHON
        # EXAMPLE_VERSION_INFO shows you how to pass a value into the C++ code
        # from Python.
        print(
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
        )
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}{os.sep}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",  # not used on MSVC, but no harm
            # option flags for gismo
            f"-DGISMO_WITH_PYBIND11=ON",
            f"-DGISMO_BUILD_EXAMPLES=OFF",
            # find_package(pybind11) hint
            f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
        ]
        # extra cmake args
        # cmake_args.extend(ext.extra_args["cmake_args"])

        build_args = []
        # Adding CMake arguments set as environment variable
        # (needed e.g. to build for ARM OSx on conda-forge)
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [
                item for item in os.environ["CMAKE_ARGS"].split(" ") if item
            ]

        if self.compiler.compiler_type != "msvc":
            # Using Ninja-build since it a) is available as a wheel and b)
            # multithreads automatically. MSVC would require all variables be
            # exported for Ninja to pick it up, which is a little tricky to do.
            # Users can override the generator with CMAKE_GENERATOR in CMake
            # 3.15+.
            if not cmake_generator or cmake_generator == "Ninja":
                try:
                    import ninja  # noqa: F401

                    ninja_executable_path = Path(ninja.BIN_DIR) / "ninja"
                    cmake_args += [
                        "-GNinja",
                        "-DCMAKE_MAKE_PROGRAM"
                        f":FILEPATH={ninja_executable_path}",
                    ]
                except ImportError:
                    pass

        else:
            # Single config generators are handled "normally"
            single_config = any(
                x in cmake_generator for x in {"NMake", "Ninja"}
            )

            # CMake allows an arch-in-generator style
            # for backward compatibility
            contains_arch = any(x in cmake_generator for x in {"ARM", "Win64"})

            # Specify the arch if using MSVC generator, but only if it doesn't
            # contain a backward-compatibility arch spec already in the
            # generator name.
            if not single_config and not contains_arch:
                cmake_args += ["-A", PLAT_TO_CMAKE[self.plat_name]]

            # Multi-config generators have a different way to specify configs
            if not single_config:
                cmake_args += [
                    "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_"
                    f"{cfg.upper()}={extdir}"
                ]
                build_args += ["--config", cfg]

        if sys.platform.startswith("darwin"):
            # Cross-compile support for macOS - respect ARCHFLAGS if set
            archs = re.findall(r"-arch (\S+)", os.environ.get("ARCHFLAGS", ""))
            if archs:
                cmake_args += [
                    "-DCMAKE_OSX_ARCHITECTURES={}".format(";".join(archs))
                ]

        # Set CMAKE_BUILD_PARALLEL_LEVEL to control the parallel build level
        # across all generators.
        if "CMAKE_BUILD_PARALLEL_LEVEL" not in os.environ:
            # self.parallel is a Python 3 only way to set parallel jobs by hand
            # using -j in the build_ext call,
            # not supported by pip or PyPA-build.
            if hasattr(self, "parallel") and self.parallel:
                # CMake 3.12+ only.
                build_args += [f"-j{self.parallel}"]

            # if len(ext.extra_args["build_args"]) != 0:
            #    build_args.extend(ext.extra_args["build_args"])

        build_temp = Path(self.build_temp) / ext.name
        if not build_temp.exists():
            build_temp.mkdir(parents=True)

        # provide standard include path for third party libraries
        gs_external = Path().cwd() / "external"
        eigen_sym_path = Path(gs_external / "Eigen")
        eigen_sym_path.symlink_to(gs_external / "gsEigen")

        subprocess.run(
            ["cmake", ext.sourcedir] + cmake_args, cwd=build_temp, check=True
        )
        subprocess.run(
            ["cmake", "--build", "."] + build_args, cwd=build_temp, check=True
        )

        # delete eigen_sym_path
        eigen_sym_path.unlink()


# use README as long description
with open("README.md") as readme:
    long_description = readme.read()

setup(
    name="pygismo",
    version="0.0.0",
    author="Angelos Mantzaflaris",
    author_email="angelos.mantzaflaris@inria.fr",
    description="G+Smo (Geometry + Simulation Modules)",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gismo/gismo",
    install_requires=[
        "numpy",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)"
        "Programming Language :: Python",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
    ],
    ext_modules=[CMakeExtension("pygismo")],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)
