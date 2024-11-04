
import os
import time
import shutil

from pylightcurve.processes.files import open_dict, download

try:
    import zipfile
    download_zip = True
except:
    download_zip = False


def _setup_database(plc_data, database_name):

        print('pylightcurve: Checking {0} database...'.format(database_name))

        # define paths

        database_directory_path = os.path.join(plc_data.databases_directory_path, database_name)
        database_file_path = os.path.join(plc_data.databases_directory_path, database_name + '.pickle')
        database_link_file_path = os.path.join(plc_data.databases_directory_path, database_name + '_link.txt')
        database_file_path_new = os.path.join(plc_data.databases_directory_path, database_name + '_new.pickle')
        database_file_path_old = os.path.join(plc_data.databases_directory_path, database_name + '_old.pickle')
        last_update_file_path = os.path.join(plc_data.databases_directory_path, '{0}_last_update.txt'.format(database_name))

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

            if not download(plc_data.databases[plc_data.version][database_name], database_file_path):
                print('\n{0} features cannot be used.'.format(database_name))
                return False
            else:
                shutil.copy(database_file_path, database_file_path_old)
                w = open(database_link_file_path, 'w')
                w.write(plc_data.databases[plc_data.version][database_name])
                w.close()

                try:
                    new_database = open_dict(database_file_path)
                    download(new_database['zipfile'], database_directory_path + '.zip')
                    new_database = zipfile.ZipFile(database_directory_path + '.zip', 'r')
                    here = os.path.abspath('.')
                    os.chdir(plc_data.databases_directory_path)
                    new_database.extractall()
                    os.chdir(here)
                    os.remove(database_directory_path + '.zip')
                except Exception as e:
                    print('Could not download zip.', e)
                    pass

        # check if everything exists, if not reset database

        # download database if there is an update

        if plc_data.databases[plc_data.version][database_name] != open(database_link_file_path).read():

            if not download(plc_data.databases[plc_data.version][database_name], database_file_path_new):
                pass
            else:
                shutil.move(database_file_path, database_file_path_old)
                shutil.move(database_file_path_new, database_file_path)

                w = open(database_link_file_path, 'w')
                w.write(plc_data.databases[plc_data.version][database_name])
                w.close()

        # download database if there is an update

        # check all files in database, remove files that need to be updated

        current_database = open_dict(database_file_path_old)
        new_database = open_dict(database_file_path)

        for dbx_file in current_database['files']:

            if dbx_file not in new_database['files']:
                try:
                    os.remove(os.path.join(plc_data.databases_directory_path,
                                           new_database['files'][dbx_file]['local_path']))
                except:
                    pass
            elif new_database['files'][dbx_file]['link'] != current_database['files'][dbx_file]['link']:
                try:
                    os.remove(os.path.join(plc_data.databases_directory_path,
                                           new_database['files'][dbx_file]['local_path']))
                except:
                    pass

        # check for updates, remove files that need to be updated

        # download missing files

        final_check = True

        for dbx_file in new_database['files']:
            if not os.path.isfile(os.path.join(plc_data.databases_directory_path,
                                               new_database['files'][dbx_file]['local_path'])):
                try:
                    os.remove(last_update_file_path)
                except:
                    pass
                if not download(new_database['files'][dbx_file]['link'],
                                os.path.join(plc_data.databases_directory_path,
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
                                        os.path.join(plc_data.databases_directory_path,
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
            if current_database != new_database:
                shutil.copy(database_file_path, database_file_path_old)
            return database_directory_path
