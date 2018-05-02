from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import os
import glob
import re
import sys

# if sys.platform == 'darwin':
#    from distutils import sysconfig
#    vars = sysconfig.get_config_vars()
#    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-dynamiclib')


if 'CC' in os.environ:
    print("Using CC={}".format(os.environ['CC']))
else:
    os.environ["CC"] = "gcc"
    print("Using CC={} (set by setup.py)".format(os.environ['CC']))

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
print(MODULE_DIR)
SRC_DIR = os.path.join(MODULE_DIR, "fidimag")
SUNDIALS_DIR = os.path.join(SRC_DIR, "common", "sundials")
NEB_DIR = os.path.join(SRC_DIR, "common", "neb")
NEBM_DIR = os.path.join(SRC_DIR, "common", "neb_method")
ATOM_DIR = os.path.join(SRC_DIR, "atomistic", "lib")
COMMON_DIR = os.path.join(SRC_DIR, "common", "lib")
MICRO_DIR = os.path.join(SRC_DIR, "micro", "lib")
BARYAKHTAR_DIR = os.path.join(MICRO_DIR, "baryakhtar")
DEMAG_DIR = os.path.join(SRC_DIR, "common", "dipolar")
FMM_DIR = os.path.join(SRC_DIR, "atomistic", "lib", "fmmlib")
LOCAL_DIR = os.path.join(MODULE_DIR, "local")
INCLUDE_DIR = os.path.join(LOCAL_DIR, "include")
LIB_DIR = os.path.join(LOCAL_DIR, "lib")
print("LIB_DIR={}".format(LIB_DIR))


pkg_init_path = os.path.join(
    os.path.dirname(__file__), 'fidimag', '__init__.py')


def get_version():
    with open(pkg_init_path) as f:
        for line in f:
            m = re.match(r'''__version__\s*=\s*(['"])(.+)\1''', line.strip())
            if m:
                return m.group(2)
    raise Exception("Couldn't find __version__ in %s" % pkg_init_path)


version = get_version()


def glob_cfiles(path, excludes, extension="*.c"):
    cfiles = []
    for cfile in glob.glob(os.path.join(path, extension)):
        filename = os.path.basename(cfile)
        if not filename in tuple(excludes):
            cfiles.append(cfile)
    return cfiles


sources = []
sources.append(os.path.join(ATOM_DIR, 'clib.pyx'))
sources += glob_cfiles(ATOM_DIR, excludes=["clib.c"])

common_sources = []
common_sources.append(os.path.join(COMMON_DIR, 'common_clib.pyx'))
common_sources += glob_cfiles(COMMON_DIR, excludes=["common_clib.c"])

cvode_sources = []
cvode_sources.append(os.path.join(SUNDIALS_DIR, 'cvode.pyx'))
cvode_sources += glob_cfiles(SUNDIALS_DIR, excludes=["cvode.c"])

baryakhtar_sources = []
baryakhtar_sources.append(os.path.join(BARYAKHTAR_DIR, 'baryakhtar_clib.pyx'))
baryakhtar_sources += glob_cfiles(BARYAKHTAR_DIR,
                                  excludes=["baryakhtar_clib.c"])

micro_sources = []
micro_sources.append(os.path.join(MICRO_DIR, 'micro_clib.pyx'))
micro_sources += glob_cfiles(MICRO_DIR, excludes=["micro_clib.c"])

# NEB Method ------------------------------------------------------------------

nebm_sources = []
nebm_sources.append(os.path.join(NEBM_DIR,
                                 'nebm_clib.pyx'))
nebm_sources += glob_cfiles(NEBM_DIR,
                            excludes=["nebm_clib.c"])

nebm_spherical_sources = []
nebm_spherical_sources.append(os.path.join(NEBM_DIR,
                                           'nebm_spherical_clib.pyx'))
nebm_spherical_sources += glob_cfiles(NEBM_DIR,
                                      excludes=["nebm_spherical_clib.c"])

nebm_geodesic_sources = []
nebm_geodesic_sources.append(os.path.join(NEBM_DIR,
                                          'nebm_geodesic_clib.pyx'))
nebm_geodesic_sources += glob_cfiles(NEBM_DIR,
                                     excludes=["nebm_geodesic_clib.c"])

nebm_cartesian_sources = []
nebm_cartesian_sources.append(os.path.join(NEBM_DIR,
                                           'nebm_cartesian_clib.pyx'))
nebm_cartesian_sources += glob_cfiles(NEBM_DIR,
                                      excludes=["nebm_cartesian_clib.c"])

# -----------------------------------------------------------------------------

dipolar_sources = []
dipolar_sources.append(os.path.join(DEMAG_DIR, 'dipolar.pyx'))
dipolar_sources += glob_cfiles(DEMAG_DIR, excludes=["dipolar.c"])

fmm_sources = []
fmm_sources.append(os.path.join(FMM_DIR, 'fmmlib.pyx'))
fmm_sources += glob_cfiles(FMM_DIR, excludes=[], extension="*.cpp")

print(fmm_sources)

com_libs = ['m', 'fftw3_omp', 'fftw3', 'sundials_cvodes',
            'sundials_nvecserial', 'sundials_nvecopenmp', 'blas', 'lapack']

com_args = ['-O3']
# rpath is the path relative to the compiled shared object files (e.g. clib.so, etc)
# which the dynamic linker looks for the linked libraries (e.g. libsundials_*.so) in.
# We need to set it relatively in order for it to be preserved if the parent directory is moved
# hence why it is a 'relative'(r) path. Here the relative path is with respect to
# the fidimag/fidimag/extensions directory.
RPATH = '../../local/lib'
com_link = ['-Wl,-rpath,{}'.format(LIB_DIR)]
lib_paths = [LIB_DIR]


if 'icc' in os.environ['CC']:
    com_args.append('-openmp')
    com_link.append('-openmp')
else:
    com_args.append('-fopenmp')
    com_link.append('-fopenmp')


com_inc = [numpy.get_include(), INCLUDE_DIR]

if 'SUNDIALS_DIR' in os.environ:
    lib_paths.append(os.environ['SUNDIALS_DIR'])
    com_inc.append(os.environ['SUNDIALS_INC'])

if 'FFTW_DIR' in os.environ:
    lib_paths.append(os.environ['FFTW_DIR'])
    com_inc.append(os.environ['FFTW_INC'])

com_args_cpp = com_args.copy()
com_args_cpp.append('-std=c++14')

com_args.append('-std=c99')
com_inc_cpp = com_inc.copy()
com_link_cpp = com_link.copy()
com_libs_cpp = com_libs.copy()
lib_paths_cpp = lib_paths.copy()

ext_modules = [
    Extension("fidimag.extensions.clib",
              sources=sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.common_clib",
              sources=common_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.cvode",
              sources=cvode_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.micro_clib",
              sources=micro_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.nebm_clib",
              sources=nebm_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.nebm_spherical_clib",
              sources=nebm_spherical_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.nebm_geodesic_clib",
              sources=nebm_geodesic_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.nebm_cartesian_clib",
              sources=nebm_cartesian_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.cvode",
              sources=cvode_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.baryakhtar_clib",
              sources=baryakhtar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.micro_clib",
              sources=micro_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.dipolar",
              sources=dipolar_sources,
              include_dirs=com_inc,
              libraries=com_libs,
              library_dirs=lib_paths, runtime_library_dirs=lib_paths,
              extra_compile_args=com_args,
              extra_link_args=com_link,
              ),
    Extension("fidimag.extensions.fmm",
              sources=fmm_sources,
              include_dirs=com_inc_cpp,
              libraries=com_libs_cpp,
              library_dirs=lib_paths_cpp, runtime_library_dirs=lib_paths_cpp,
              extra_compile_args=com_args_cpp,
              extra_link_args=com_link_cpp,
              language="c++",
              )
]

setup(
    name='fidimag',
    version=version,
    description='Finite difference micromagnetic code',
    packages=['fidimag',
              'fidimag.atomistic',
              'fidimag.micro',
              'fidimag.extensions',
              'fidimag.common',
              ],
    ext_modules=cythonize(ext_modules),
)
