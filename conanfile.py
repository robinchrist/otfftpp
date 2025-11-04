from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.files import copy
import os

class otfftppRecipe(ConanFile):
    name = "otfftpp"
    version = "0.0.6"
    package_type = "library"

    license = "MIT"
    url = "https://github.com/robinchrist/otfftpp"
    description = "OTFFT is a high-speed FFT library using the Stockham's algorithm and SIMD"
    topics = ("FFT", "SIMD")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "with_openmp": [True, False],
        "build_testing": [True, False]
    }
    default_options = {
        "shared": False,
        "fPIC": True,
        "with_openmp": False,
        "build_testing": False
    }

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "tests/*"

    @property
    def _openmp_flags(self):
        if self.settings.compiler == "clang":
            return ["-fopenmp=libomp"]
        elif self.settings.compiler == "apple-clang":
            return ["-Xclang", "-fopenmp"]
        elif self.settings.compiler == "gcc":
            return ["-fopenmp"]
        elif self.settings.compiler == "intel-cc":
            return ["-Qopenmp"]
        elif self.settings.compiler == "sun-cc":
            return ["-xopenmp"]
        elif is_msvc(self):
            return ["-openmp"]
        return None

    def config_options(self):
        if self.settings.os == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        if self.options.shared:
            self.options.rm_safe("fPIC")

    def requirements(self):
        self.requires("simde/0.8.2")

        if self.options.with_openmp and self.settings.compiler in ["clang", "apple-clang"]:
            self.requires("llvm-openmp/17.0.6", transitive_headers=True, transitive_libs=True)

    def build_requirements(self):
        self.test_requires("boost/[^1.82.0]")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()

        tc = CMakeToolchain(self)

        tc.variables["OTFFTPP_WITH_OPENMP"] = self.options.with_openmp
        tc.variables["OTFFTPP_BUILD_TESTS"] = self.options.build_testing

        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        copy(self, "LICENSE", src=self.source_folder, dst=os.path.join(self.package_folder, "licenses"))
        
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["otfftpp"]

        if self.options.with_openmp:
            if self.settings.compiler in ["clang", "apple-clang"]:
                self.cpp_info.requires.append("llvm-openmp::llvm-openmp")
            openmp_flags = self._openmp_flags
            self.cpp_info.cflags = openmp_flags
            self.cpp_info.cxxflags = openmp_flags
            self.cpp_info.sharedlinkflags = openmp_flags
            self.cpp_info.exelinkflags = openmp_flags
