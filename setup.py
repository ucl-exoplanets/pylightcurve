
import os
import glob

from setuptools import setup

name = 'pylightcurve'
description = 'A python package for modeling and analysing transit light-curves.'
url = 'https://https://github.com/ucl-exoplanets/pylightcurve'
install_requires = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'requirements.txt')).read().split('\n')

os.chdir(os.path.abspath(os.path.dirname(__file__)))

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
