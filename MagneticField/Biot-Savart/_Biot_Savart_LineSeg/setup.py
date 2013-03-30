from distutils.core import setup, Extension
import os.path

module1 = Extension('BiotSavartLineSeg',
                    define_macros=[('MAJOR_VERSION','1'),
                                   ('MINOR_VERSION','0')],
                    include_dirs=[('/usr/local/include'),
                                  ('/usr/include'),
                                  (os.path.abspath('~/Projects/Utilities/clibs/include'))],
                    libraries=[('m'),('gsl'),('blas')],
                    library_dirs=[('/usr/local/lib64'),
                                  ('/usr/local/lib'),
                                  ('/usr/lib64'),
                                  (os.path.abspath('~/Projects/Utilities/clibs/lib64'))],
                    sources=[('BiotSavart_LineSeg.c'),('Vector.c'),('LineSeg.c')])

setup (name='BiotSavart Line Segment Helpers',
       version='1.0',
       description='This is the base c package for the Kerenel of the Biot Savart Line Seg computation',
       ext_modules=[module1]
)
