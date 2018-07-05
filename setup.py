from setuptools import setup
import os
import sys
import glob
import time
import codecs
import shutil
import socket

if sys.version_info[0] > 2:
    from urllib.request import urlopen, urlretrieve
else:
    from urllib import urlopen, urlretrieve
    input = raw_input

name = 'pylightcurve'
description = 'A python package for modeling and analysing transit light-curves.'
url = 'https://https://github.com/ucl-exoplanets/pylightcurve'
install_requires = ['matplotlib', 'numpy', 'exodata', 'emcee', 'seaborn', 'ephem', 'astropy', 'scipy']

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

w = open('MANIFEST.in', 'w')
for i in subdirs_to_include:
    w.write('include ' + os.path.join(i, '*') + ' \n')

for i in files_to_include:
    w.write('include ' + i + ' \n')

w.close()

with codecs.open('README.md', encoding='utf-8') as f:
    long_description = f.read()

version = ' '
for i in open(os.path.join(name, '__init__.py')):
    if len(i.split('__version__')) > 1:
        version = i.split()[-1][1:-1]

setup(
    name=name,
    version=version,
    description=description,
    long_description=long_description,
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
                 'Programming Language :: Python :: 3.6',
                 ],
    packages=[name],
    install_requires=install_requires,
    include_package_data=True,
    zip_safe=False,
)


def intsall_phoenix():
    force_update = False

    database_location = os.path.join(os.path.expanduser('~'), '.pylightcurve')
    if not os.path.isdir(database_location):
        os.mkdir(database_location)

    phoenix_database_location = os.path.join(database_location, 'phoenix_database')
    phoenix_database_last_update_file_path = os.path.join(database_location, 'phoenix_database_last_update.txt')
    phoenix_database_zip_file_path = os.path.join(database_location, 'phoenix_database.zip')
    phoenix_database_zip_file_path2 = os.path.join(database_location, 'phoenix_database2.zip')
    phoenix_database_zip_file_path3 = os.path.join(database_location, 'phoenix_database3.zip')

    phoenix_database_url = 'https://www.dropbox.com/sh/39et0eg8akk4ga8/AADcCAEVWirSGwH8bFgj2rq1a?dl=1'
    phoenix_database_url2 = 'https://www.dropbox.com/sh/pyx89lnlxa6dqvs/AADVuMRBtSGIhV09NP1XWiMOa?dl=1'
    phoenix_database_url3 = 'https://www.dropbox.com/sh/4tiqkzb60hmdewy/AAAeKprAX9LiH5phYM-RiFR2a?dl=1'

    if os.path.isdir(phoenix_database_location):
        if force_update or len(glob.glob(os.path.join(phoenix_database_location, '*'))) == 0:
            shutil.rmtree(phoenix_database_location)
        else:
            if not os.path.isfile(phoenix_database_last_update_file_path):
                shutil.rmtree(phoenix_database_location)
            elif int(open(phoenix_database_last_update_file_path).readlines()[0]) < 180428:
                shutil.rmtree(phoenix_database_location)

    if not os.path.isdir(phoenix_database_location):

        try:
            if input('Downloading phoenix database (3.9GB, in three parts)... '
                     'proceed with download now? (y/n):') == 'y':
                def reporthook(count, block_size, size):
                    progress_size = int(count * block_size / (1024 * 1024))
                    percent = int(progress_size * 100 / 1791)
                    sys.stdout.write('\rPart 1... {0}%, {1} MB'.format(percent, progress_size))
                    sys.stdout.flush()

                socket.setdefaulttimeout(500000)
                urlretrieve(phoenix_database_url, phoenix_database_zip_file_path, reporthook)
                socket.setdefaulttimeout(30)

                print('')

                def reporthook(count, block_size, size):
                    progress_size = int(count * block_size / (1024 * 1024))
                    percent = int(progress_size * 100 / 1119)
                    sys.stdout.write('\rPart 2... {0}%, {1} MB'.format(percent, progress_size))
                    sys.stdout.flush()

                socket.setdefaulttimeout(500000)
                urlretrieve(phoenix_database_url2, phoenix_database_zip_file_path2, reporthook)
                socket.setdefaulttimeout(30)
                w = open(phoenix_database_last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

                print('')

                def reporthook(count, block_size, size):
                    progress_size = int(count * block_size / (1024 * 1024))
                    percent = int(progress_size * 100 / 1081)
                    sys.stdout.write('\rPart 3... {0}%, {1} MB'.format(percent, progress_size))
                    sys.stdout.flush()

                socket.setdefaulttimeout(500000)
                urlretrieve(phoenix_database_url3, phoenix_database_zip_file_path3, reporthook)
                socket.setdefaulttimeout(30)
                w = open(phoenix_database_last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

                os.system('unzip {0} -d {1}'.format(phoenix_database_zip_file_path, phoenix_database_location))
                os.system('rm {0}'.format(phoenix_database_zip_file_path))
                os.system('unzip {0} -d {1}'.format(phoenix_database_zip_file_path2, phoenix_database_location))
                os.system('rm {0}'.format(phoenix_database_zip_file_path2))
                os.system('unzip {0} -d {1}'.format(phoenix_database_zip_file_path3, phoenix_database_location))
                os.system('rm {0}'.format(phoenix_database_zip_file_path3))

        except:
            print('Downloading phoenix database failed. A download will be attempted next time.')
            pass

    if (not os.path.isdir(phoenix_database_location) or
            len(glob.glob(os.path.join(phoenix_database_location, '*'))) == 0):
        print('Phoenix features cannot be used.')


intsall_phoenix()
