from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('vcalc_cy.pyx'))
