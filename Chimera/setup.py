from setuptools import setup, find_packages
from Cython.Build import cythonize

setup(
    name="Chimera",
    version="0.1",
    description="A magma ocean box modeling tool.",
    author="Scott D. Hull",
    license="BSD 3-Clause License",
    packages=find_packages(),
    ext_modules=cythonize("conduction.pyx"),
)