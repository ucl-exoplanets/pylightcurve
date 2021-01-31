from setuptools import setup, Extension
import os
import sys
import glob
import time
import codecs
import shutil
import socket
import numpy

if sys.version_info[0] > 2:
    from urllib.request import urlopen, urlretrieve
else:
    from urllib import urlopen, urlretrieve
    input = raw_input

name = 'pylightcurve'
description = 'A python package for modeling and analysing transit light-curves.'
url = 'https://https://github.com/ucl-exoplanets/pylightcurve'
install_requires = ['matplotlib', 'numpy', 'exodata', 'emcee', 'seaborn', 'astropy', 'scipy', 'sklearn', 'astroquery']

os.chdir(os.path.abspath(os.path.dirname(__file__)))

subdirs_to_include = []
for x in os.walk(name):
    if os.path.isdir(x[0]):
        if x[0] != name:
            subdirs_to_include.append(x[0])

files_to_include = []
for x in glob.glob(os.path.join(name, '*')):
    if os.path.isfile(x):
        if x.split('.')[-1] not in ['py']:
            files_to_include.append(os.path.join(name, os.path.split(x)[1]))

files_to_include.append('README.md')
files_to_include.append('LICENSE')
files_to_include.append('readme.md')
files_to_include.append('licence')

w = open('MANIFEST.in', 'w')
for i in subdirs_to_include:
    w.write('include ' + os.path.join(i, '*') + ' \n')

for i in files_to_include:
    w.write('include ' + i + ' \n')

w.close()

version = ' '
for i in open(os.path.join(name, '__init__.py')):
    if len(i.split('__version__')) > 1:
        version = i.split()[-1][1:-1]

setup(
    name=name,
    version=version,
    description=description,
    long_description='Visit https://github.com/ucl-exoplanets/pylightcurve',
    url=url,
    author='Angelos Tsiaras',
    author_email='aggelostsiaras@gmail.com',
    license='MIT',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Operating System :: MacOS :: MacOS X',
                 'Programming Language :: Python :: 3.7',
                 ],
    packages=[name],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
)
