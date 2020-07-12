from distutils.core import setup
from Cython.Build import cythonize


setup(
    name = 'dse',
    ext_modules = cythonize("DSEQin.pyx"),
)

