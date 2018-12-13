#
#
# package_database_location = os.path.join(os.path.expanduser('~'), '.pylightcurve')
# if not os.path.isdir(package_database_location):
#     os.mkdir(package_database_location)
#
#
# def database(database_name, force_update=False):
#
#     database_info_file_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), '_0database.pickle')
#     database_location = os.path.join(package_database_location, '{0}_database'.format(database_name))
#     database_last_update_file_path = os.path.join(package_database_location,
#                                                   '{0}_database_last_update.txt'.format(database_name))
#
#     if os.path.isdir(database_location):
#         if force_update or len(glob.glob(os.path.join(database_location, '*'))) == 0:
#             shutil.rmtree(database_location)
#             os.mkdir(database_location)
#             update = True
#         else:
#             if not os.path.isfile(database_last_update_file_path):
#                 update = True
#             elif int(open(database_last_update_file_path).readlines()[0]) < 181212:
#                 update = True
#             else:
#                 update = False
#     else:
#         os.mkdir(database_location)
#         update = True
#
#     if update and database_name == 'phoenix':
#         if input('Downloading phoenix database (up to 5GB)... proceed with download now? (y/n): ') == 'y':
#             update = True
#         else:
#             update = False
#
#     if update:
#         try:
#             print('Downloading {0} database...'.format(database_name))
#
#             dbx_files = pickle.load(open(database_info_file_path, 'rb'))['{0}_database'.format(database_name)]
#
#             for i in glob.glob(os.path.join(database_location, '*')):
#                 if os.path.split(i)[1] not in dbx_files:
#                     os.remove(i)
#
#             for i in dbx_files:
#                 if not os.path.isfile(os.path.join(package_database_location, dbx_files[i]['local_path'])):
#                     print(i)
#                     urlretrieve(dbx_files[i]['link'],
#                                 os.path.join(package_database_location, dbx_files[i]['local_path']))
#
#             if database_name == 'clablimb':
#                 xx = pickle.load(open(glob.glob(os.path.join(database_location, '*'))[0], 'rb'))
#                 for i in xx:
#                     w = open(os.path.join(database_location, i), 'w')
#                     w.write(xx[i])
#                     w.close()
#
#             w = open(database_last_update_file_path, 'w')
#             w.write(time.strftime('%y%m%d'))
#             w.close()
#
#         except:
#             print('Downloading {0} database failed. A download will be attempted next time.'.format(database_name))
#             pass
#
#     if (not os.path.isdir(database_location) or
#             len(glob.glob(os.path.join(database_location, '*'))) == 0):
#
#         print('{0} features cannot be used.'.format(database_name))
#
#     else:
#
#         return database_location
#
#
# def clablimb_database(force_update=False):
#     return database('clablimb', force_update=force_update)
#
#
# def phoenix_database(force_update=False):
#     return database('phoenix', force_update=force_update)
#
#
# # oec database
#
# oec_database_location = os.path.join(package_database_location, 'oec_database')
# if not os.path.isdir(oec_database_location):
#     os.mkdir(oec_database_location)
#
# systems_backup_file_path = os.path.join(oec_database_location, 'systems_backup.xml.gz')
# systems_file_path = os.path.join(oec_database_location, 'systems.xml.gz')
# systems_last_update_file_path = os.path.join(oec_database_location, 'systems_last_update.txt')
#
# systems_backup_url = 'https://www.dropbox.com/s/ajt34peoq9u7p54/systems_backup.xml.gz?raw=1'
# systems_url = 'https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz'
#
#
# def oec_database():
#
#     if not os.path.isfile(systems_backup_file_path):
#
#         try:
#             print('Downloading backup OEC...')
#             socket.setdefaulttimeout(5)
#             urlretrieve(systems_backup_url, systems_backup_file_path)
#             socket.setdefaulttimeout(30)
#         except:
#             print('Downloading backup OEC failed. A download will be attempted next time.')
#             pass
#
#     try:
#         date = time.strftime('%y%m%d')
#         update = False
#         if not os.path.isfile(systems_last_update_file_path):
#             update = True
#         elif not os.path.isfile(systems_file_path):
#             update = True
#         elif int(open(systems_last_update_file_path).readlines()[0]) < int(date):
#             update = True
#
#         if update:
#
#             print('Updating OEC...')
#
#             try:
#                 socket.setdefaulttimeout(5)
#                 urlretrieve(systems_url, systems_file_path)
#                 socket.setdefaulttimeout(30)
#                 w = open(systems_last_update_file_path, 'w')
#                 w.write(date)
#                 w.close()
#             except:
#                 print('Updating OEC failed. An update will be attempted next time.')
#                 pass
#     except:
#         pass
#
#     if os.path.isfile(systems_file_path):
#
#         return systems_file_path
#
#     elif os.path.isfile(systems_backup_file_path):
#
#         print('No updated OEC found. Backup OEC will be used.')
#         return systems_last_update_file_path
#
#     else:
#
#         print('No OEC found. OEC features cannot be used.')
