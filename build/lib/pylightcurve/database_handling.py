import os
import sys
import glob
import gzip
import time
import shutil
import socket

if sys.version_info[0] > 2:
    from urllib.request import urlopen, urlretrieve
else:
    from urllib import urlopen, urlretrieve
    input = raw_input

database_location = os.path.join(os.path.expanduser('~'), '.pylightcurve')
if not os.path.isdir(database_location):
    os.mkdir(database_location)


# clablimb database

clablimb_database_location = os.path.join(database_location, 'clablimb_database')
clablimb_database_last_update_file_path = os.path.join(database_location, 'clablimb_database_last_update.txt')
clablimb_database_zip_file_path = os.path.join(database_location, 'clablimb_database.zip')

clablimb_database_url = 'https://www.dropbox.com/sh/c1dnqgfht2ara32/AADS4yus59DJj1ifImDexYkba?dl=1'


def clablimb_database(force_update=False):

    if os.path.isdir(clablimb_database_location):
        if force_update or len(glob.glob(os.path.join(clablimb_database_location, '*'))) == 0:
            shutil.rmtree(clablimb_database_location)
        else:
            if not os.path.isfile(clablimb_database_last_update_file_path):
                shutil.rmtree(clablimb_database_location)
            elif int(open(clablimb_database_last_update_file_path).readlines()[0]) < 180428:
                shutil.rmtree(clablimb_database_location)

    if not os.path.isdir(clablimb_database_location):

        try:
            print('Downloading clablimb database...')

            urlretrieve(clablimb_database_url, clablimb_database_zip_file_path)
            w = open(clablimb_database_last_update_file_path, 'w')
            w.write(time.strftime('%y%m%d'))
            w.close()

            os.system('unzip {0} -d {1}'.format(clablimb_database_zip_file_path, clablimb_database_location))
            os.system('rm {0}'.format(clablimb_database_zip_file_path))

        except:
            print('Downloading clablimb database failed. A download will be attempted next time.')
            pass

    if (not os.path.isdir(clablimb_database_location) or
            len(glob.glob(os.path.join(clablimb_database_location, '*'))) == 0):

        print('Clablimb features cannot be used.')

    else:

        return clablimb_database_location


# phoenix database

phoenix_database_location = os.path.join(database_location, 'phoenix_database')
phoenix_database_last_update_file_path = os.path.join(database_location, 'phoenix_database_last_update.txt')
phoenix_database_zip_file_path = os.path.join(database_location, 'phoenix_database.zip')

phoenix_database_url = 'https://www.dropbox.com/sh/39et0eg8akk4ga8/AADcCAEVWirSGwH8bFgj2rq1a?dl=1'


def phoenix_database(force_update=False):

    if os.path.isdir(phoenix_database_location):
        if force_update or len(glob.glob(os.path.join(phoenix_database_location, '*'))) == 0:
            shutil.rmtree(phoenix_database_location)
        else:
            if not os.path.isfile(phoenix_database_last_update_file_path):
                shutil.rmtree(phoenix_database_location)
            elif int(open(phoenix_database_last_update_file_path).readlines()[0]) < 180428:
                shutil.rmtree(phoenix_database_location)

    if not os.path.isdir(phoenix_database_location):

        # try:
            if input('Downloading phoenix database (4.2GB)... proceed with download now? (y/n):') == 'y':

                def reporthook(count, block_size, size):
                    progress_size = int(count * block_size / (1024 * 1024))
                    percent = int(progress_size * 100 / 3992)
                    sys.stdout.write('\r... {0}%, {1} MB'.format(percent, progress_size))
                    sys.stdout.flush()

                socket.setdefaulttimeout(500000)
                urlretrieve(phoenix_database_url, phoenix_database_zip_file_path, reporthook)
                socket.setdefaulttimeout(30)
                w = open(phoenix_database_last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

                os.system('unzip {0} -d {1}'.format(phoenix_database_zip_file_path, phoenix_database_location))
                os.system('rm {0}'.format(phoenix_database_zip_file_path))

        # except:
        #     print('Downloading phoenix database failed. A download will be attempted next time.')
        #     pass

    if (not os.path.isdir(phoenix_database_location) or
            len(glob.glob(os.path.join(phoenix_database_location, '*'))) == 0):

        print('Phoenix features cannot be used.')

    else:

        return phoenix_database_location


# oec database

oec_database_location = os.path.join(database_location, 'oec_database')
if not os.path.isdir(oec_database_location):
    os.mkdir(oec_database_location)

systems_backup_file_path = os.path.join(oec_database_location, 'systems_backup.xml.gz')
systems_file_path = os.path.join(oec_database_location, 'systems.xml.gz')
systems_last_update_file_path = os.path.join(oec_database_location, 'systems_last_update.txt')

systems_backup_url = 'https://www.dropbox.com/s/ajt34peoq9u7p54/systems_backup.xml.gz?raw=1'
systems_url = 'https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz'

aaaaaa='https://www.dropbox.com/sh/08kfai6wdp2y6p0/AACS9HIxUzK4lmKW5LdJvLrea?dl=1'
def oec_database():

    if not os.path.isfile(systems_backup_file_path):

        try:
            print('Downloading backup OEC...')
            socket.setdefaulttimeout(5)
            urlretrieve(systems_backup_url, systems_backup_file_path)
            socket.setdefaulttimeout(30)
        except:
            print('Downloading backup OEC failed. A download will be attempted next time.')
            pass

    try:
        date = time.strftime('%y%m%d')
        update = False
        if not os.path.isfile(systems_last_update_file_path):
            update = True
        elif not os.path.isfile(systems_file_path):
            update = True
        elif int(open(systems_last_update_file_path).readlines()[0]) < int(date):
            update = True

        if update:

            print('Updating OEC...')

            try:
                socket.setdefaulttimeout(5)
                urlretrieve(systems_url, systems_file_path)
                socket.setdefaulttimeout(30)
                w = open(systems_last_update_file_path, 'w')
                w.write(date)
                w.close()
            except:
                print('Updating OEC failed. An update will be attempted next time.')
                pass
    except:
        pass

    if os.path.isfile(systems_file_path):

        return systems_file_path

    elif os.path.isfile(systems_backup_file_path):

        print('No updated OEC found. Backup OEC will be used.')
        return systems_last_update_file_path

    else:

        print('No OEC found. OEC features cannot be used.')
