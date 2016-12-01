from setuptools import setup, Extension
from Cython.Build import cythonize
from glob import glob

#-------------Transform Module------------------
filelist1 = ["cython/event.pyx", "src/utility.cpp"]
extensions = [
    Extension(
        'event',
                filelist1,
                language="c++",
        extra_compile_args=["-std=c++11"],
                libraries=["m"])
]

setup(
    ext_modules=cythonize(extensions)
        )

