# coding: utf8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='COVERnant',
    version='0.3.2',
    packages=['covernantlib'],
    author='Konrad U. Förstner',
    author_email='konrad@foerstner.org',
    description='A tool to generate and manipulate coverage plots of '
    'high-throughput sequencing data.',
    url='https://github.com/konrad/COVERnant',
    install_requires=[
        "pybedtools >= 0.7.10",
        "matplotlib >= 2.2.2",
        "pandas >= 0.22.0",
        "pysam >=  0.14.1",
        "numpy >= 1.14.2"
    ],
    scripts=['bin/covernant'],
    license='ISC License (ISCL)',
    long_description=open('README.rst').read(),
    classifiers=[
        'License :: OSI Approved :: ISC License (ISCL)',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
