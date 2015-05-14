try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='COVERnant',
    version='0.1dev',
    packages=['covernant'],
    author='Konrad U. Förstner',
    author_email='konrad@foerstner.org',
    description='A tool for coverage creation and manipulation.',
    url='',
    install_requires=[
        "pysam >= 0.8.1"
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
