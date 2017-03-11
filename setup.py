from setuptools import setup
import codecs
import os

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'readme.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pylightcurve',
    version='1.1.0',
    description="Modeling and analysing transit light-curves",
    long_description=long_description,
    url='https://github.com/ucl-exoplanets/pylightcurve',
    author='Angelos Tsiaras',
    author_email='aggelostsiaras@gmail.com',
    license='GPLv3',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Operating System :: MacOS :: MacOS X'
                 'Programming Language :: Python :: 2.7',
                 ],
    packages=['pylightcurve'],
    install_requires=['matplotlib', 'numpy', 'quantities'],
    include_package_data=True,
    zip_safe=False,
)
