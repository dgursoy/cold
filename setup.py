#!/usr/bin/env python3

from setuptools import setup, Extension

coldc = Extension(
    name='libcold',
    extra_compile_args=['-std=c99', '-DWIN32'],
    sources=['source/core.c'], 
    include_dirs=['include'])

setup(
    name='cold',
    version='1.0',
    description='Coded Laue Diffraction Imaging',
    author='Doga Gursoy',
    author_email='dgursoy@anl.gov',
    packages=['cold'],
    ext_modules=[coldc],
    license='BSD-3'
    )