
import os

from setuptools import setup

package = 'pylightcurve'
version = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__version__.txt')).read()
author = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__author__.txt')).read()
author_email = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__author_email__.txt')).read()
description = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__description__.txt')).read()
url = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), package, '__url__.txt')).read()

install_requires = open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'requirements.txt')).read().split('\n')
entry_point = ''

os.chdir(os.path.abspath(os.path.dirname(__file__)))

setup(
    name=package,
    version=version,
    description=description,
    long_description='Visit {0} for more details.'.format(url),
    url=url,
    author=author,
    author_email=author_email,
    license='MIT',
    classifiers=['Development Status :: 4 - Beta',
                 'Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Topic :: Scientific/Engineering :: Astronomy',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Operating System :: MacOS :: MacOS X',
                 'Programming Language :: Python :: 3.7',
                 ],
    entry_points={},
    packages=[package],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
    setup_requires=["pytest-runner"],
    tests_require=['pytest'],
)
