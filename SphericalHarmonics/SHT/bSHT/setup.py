from distutils.core import setup, Extension

#
#   It is probably a bad idea but I installed fftw3 from source
#   but somehow that didn't work. So I made a symbolic link
#   from the /usr/lib64/libfftw3.so.3.0.1 to /usr/local/lib64/libfftw3.a and libfftw3.la
#
#    I could at least complie then!
module1 = Extension('bSHT',
                    define_macros=[('MAJOR_VERSION','1'),
                                   ('MINOR_VERSION','0')],
                    include_dirs=[('/home/edlerk/Local/include')],
                    libraries=[('fftw3'),
                               ('m')],
                    library_dirs=[('/home/edlerk/Local/lib')],
                    runtime_library_dirs=[('/home/edlerk/Local/lib')],
                    sources=[('SHTmoduleBACKEND.c'),
                             ('SHTmodule.c'),
                             ('FST_semi_fly.c'),
                             ('naive_synthesis.c'),
                             ('pmls.c'),
                             ('cospmls.c'),
                             ('primitive.c'),
                             ('seminaive.c'),
                             ('makeweights.c')])

setup (name='Spherical Harmonic Transform',
       version='1.0',
       description='This is the base c package for SHT',
       ext_modules=[module1]
)
