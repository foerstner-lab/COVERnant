# coding: utf8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='COVERnant',
    version='0.1dev',
    packages=['covernantlib'],
    author='Konrad U. Förstner',
    author_email='konrad@foerstner.org',
    description='A tool for DNA/RNA coverage creation and manipulation.',
    url='',
    install_requires=[
        "matplotlib >= 1.4.3",
        "pandas >= 0.15.2",
        "pysam >= 0.8.1",
        "numpy >= 1.9.2"
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
