
import os
import glob
import time
import shutil
from scipy.interpolate import interp1d

from pylightcurve.processes.files import open_dict, open_yaml, download
from pylightcurve import __version__

try:
    import zipfile
    download_zip = True
except:
    download_zip = False

databases_file = '__databases__.pickle'
package_name = 'pylightcurve4'
github_link = 'https://github.com/ucl-exoplanets/pylightcurve/raw/master/pylightcurve/__databases__.pickle?raw=true'


class PlcData:

    def __init__(self, _reset=False, _test=False):

        self.package_name = package_name
        self.version = '.'.join(__version__.split('.')[:2])

        self.build_in_databases_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), databases_file)

        self.databases_directory_path = os.path.join(os.path.abspath(os.path.expanduser('~')),
                                                     '.{0}'.format(self.package_name))

        self.databases_file_path = os.path.join(self.databases_directory_path, databases_file)
        self.databases_file_path_new = os.path.join(self.databases_directory_path, databases_file + '_new')

        # initiate databases

        if not os.path.isdir(self.databases_directory_path):
            os.mkdir(self.databases_directory_path)

        if not os.path.isfile(self.databases_file_path):
            shutil.copy(self.build_in_databases_file_path, self.databases_file_path)

        # check for updates in the databases (identified on github)

        if download(github_link, self.databases_file_path_new):
            shutil.move(self.databases_file_path_new, self.databases_file_path)
        else:
            pass

        # load databases

        self.databases = open_dict(self.databases_file_path)

        self.exotethys_loaded = self._setup_database('exotethys')
        self.ephemerides_loaded = self._setup_database('ephemerides')
        self.photometry_loaded = self._setup_database('photometry')
        self.catalogues_loaded = self._setup_database('catalogues')

        self.barycenter_data = None
        self.sun_data = None
        self.ecc_data = None
        self.all_filters_data = None

    def exotethys(self):
        return self.exotethys_loaded

    def ephemeris(self):
        return self.ephemerides_loaded

    def photometry(self):
        return self.photometry_loaded

    def catalogues(self):
        return self.catalogues_loaded

    def all_filters(self):

        if not self.all_filters_data:
            self.all_filters_data = [os.path.split(ff)[1].split('.')[0]
                                     for ff in glob.glob(os.path.join(self.photometry(), '*'))]

        return self.all_filters_data

    def ecc(self):

        if not self.ecc_data:
            stars = open_yaml(os.path.join(self.catalogues(), 'ecc_stars.yaml'))
            planets = open_yaml(os.path.join(self.catalogues(), 'ecc_planets.yaml'))

            self.ecc_data = {'stars': stars, 'planets':planets}

        return self.ecc_data

    def barycentre(self, utc_jd):

        if not self.barycenter_data:
            self.barycenter_data = open_dict(os.path.join(self.ephemeris(), 'bjd_dict.pickle'))

        bjd_dict = self.barycenter_data[int(utc_jd)]

        ssb_t = bjd_dict['t']
        ssb_ra = bjd_dict['ra']
        ssb_dec = bjd_dict['dec']
        ssb_d = bjd_dict['d']
        ssb_dt = bjd_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(utc_jd)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(utc_jd)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(utc_jd)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(utc_jd)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def heliocentre(self, utc_jd):

        if not self.sun_data:
            self.sun_data = open_dict(os.path.join(self.ephemeris(), 'hjd_dict.pickle'))

        hjd_dict = self.sun_data[int(utc_jd)]

        ssb_t = hjd_dict['t']
        ssb_ra = hjd_dict['ra']
        ssb_dec = hjd_dict['dec']
        ssb_d = hjd_dict['d']
        ssb_dt = hjd_dict['dt']

        ssb_ra = interp1d(ssb_t, ssb_ra, kind='cubic')(utc_jd)
        ssb_dec = interp1d(ssb_t, ssb_dec, kind='cubic')(utc_jd)
        ssb_d = interp1d(ssb_t, ssb_d, kind='cubic')(utc_jd)
        ssb_dt = interp1d(ssb_t, ssb_dt, kind='cubic')(utc_jd)

        return ssb_ra, ssb_dec, ssb_d, ssb_dt

    def _setup_database(self, database_name):

        # define paths

        database_directory_path = os.path.join(self.databases_directory_path, database_name)
        database_file_path = os.path.join(self.databases_directory_path, database_name + '.pickle')
        database_link_file_path = os.path.join(self.databases_directory_path, database_name + '_link.txt')
        database_file_path_new = os.path.join(self.databases_directory_path, database_name + '_new.pickle')
        database_file_path_old = os.path.join(self.databases_directory_path, database_name + '_old.pickle')
        last_update_file_path = os.path.join(self.databases_directory_path, '{0}_last_update.txt'.format(database_name))

        # define paths

        # check if everything exists, if not reset database

        if not os.path.isdir(database_directory_path) or not os.path.isfile(database_file_path) or not os.path.isfile(database_link_file_path):

            try:
                shutil.rmtree(database_directory_path)
            except:
                pass

            try:
                os.remove(database_file_path)
            except:
                pass

            try:
                os.remove(database_file_path_old)
            except:
                pass

            try:
                os.remove(database_file_path_new)
            except:
                pass

            try:
                os.remove(database_link_file_path)
            except:
                pass

            try:
                os.remove(last_update_file_path)
            except:
                pass

            os.mkdir(database_directory_path)

            if not download(self.databases[self.version][database_name], database_file_path):
                print('\n{0} features cannot be used.'.format(database_name))
                return False
            else:
                shutil.copy(database_file_path, database_file_path_old)
                w = open(database_link_file_path, 'w')
                w.write(self.databases[self.version][database_name])
                w.close()

                try:
                    new_database = open_dict(database_file_path)
                    download(new_database['zipfile'], database_directory_path + '.zip')
                    new_database = zipfile.ZipFile(database_directory_path + '.zip', 'r')
                    here = os.path.abspath('.')
                    os.chdir(self.databases_directory_path)
                    new_database.extractall()
                    os.chdir(here)
                    os.remove(database_directory_path + '.zip')
                except:
                    pass

        # check if everything exists, if not reset database

        # download database if there is an update

        if self.databases[self.version][database_name] != open(database_link_file_path).read():

            if not download(self.databases[self.version][database_name], database_file_path_new):
                pass
            else:
                shutil.move(database_file_path, database_file_path_old)
                shutil.move(database_file_path_new, database_file_path)

                w = open(database_link_file_path, 'w')
                w.write(self.databases[self.version][database_name])
                w.close()

        # download database if there is an update

        # check all files in database, remove files that need to be updated

        print('Checking {0} database...'.format(database_name))

        current_database = open_dict(database_file_path_old)
        new_database = open_dict(database_file_path)

        for dbx_file in current_database['files']:

            if dbx_file not in new_database['files']:
                try:
                    os.remove(os.path.join(self.databases_directory_path,
                                           new_database['files'][dbx_file]['local_path']))
                except:
                    pass
            elif new_database['files'][dbx_file]['link'] != current_database['files'][dbx_file]['link']:
                try:
                    os.remove(os.path.join(self.databases_directory_path,
                                           new_database['files'][dbx_file]['local_path']))
                except:
                    pass

        # check for updates, remove files that need to be updated

        # download missing files

        final_check = True

        for dbx_file in new_database['files']:
            if not os.path.isfile(os.path.join(self.databases_directory_path,
                                               new_database['files'][dbx_file]['local_path'])):
                try:
                    os.remove(last_update_file_path)
                except:
                    pass
                if not download(new_database['files'][dbx_file]['link'],
                                os.path.join(self.databases_directory_path,
                                             new_database['files'][dbx_file]['local_path'])):
                    final_check = False

        # download missing files

        # update files from external links

        frequency = new_database['frequency']
        if frequency:

            try:
                last_update_date = int(open(last_update_file_path).read())
            except:
                last_update_date = 0

            today = int(time.strftime('%y%m%d'))

            if today >= last_update_date + frequency:

                for dbx_file in new_database['files']:
                    if 'external_link' in new_database['files'][dbx_file]:
                        print('\tUpdating: ', dbx_file)
                        if not download(new_database['files'][dbx_file]['external_link'],
                                        os.path.join(self.databases_directory_path,
                                                     new_database['files'][dbx_file]['local_path'])):
                            final_check = False

                w = open(last_update_file_path, 'w')
                w.write(time.strftime('%y%m%d'))
                w.close()

        # update files from external links

        if not final_check:
            print('\n{0} features cannot be used.'.format(database_name))
            return False
        else:
            return database_directory_path


plc_data = PlcData()
