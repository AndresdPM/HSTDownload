#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import argparse
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Catalogs
from astroquery.mast import Observations
from astropy.time import Time

def search_mast(ra, dec, search_width = 0.25, search_height = 0.25, filters = ['any'], project = ['HST'], t_exptime_min = 50, t_exptime_max = 2500):
   """
   This routine search for HST observations in MAST at a given position.
   """

   ra1 = ra - search_width / 2 + 0.056 / np.cos(np.deg2rad(dec))
   ra2 = ra + search_width / 2 - 0.056 / np.cos(np.deg2rad(dec))
   dec1 = dec - search_height / 2 + 0.056
   dec2 = dec + search_height / 2 - 0.056

   if type(filters) is not list:
      filters = [filters]

   obs_table = Observations.query_criteria(dataproduct_type=['image'], obs_collection=['HST'], s_ra=[ra1, ra2], s_dec=[dec1, dec2], instrument_name=['ACS/WFC', 'WFC3/UVIS'], filters = filters, project = project)

   data_products_by_obs = search_data_products_by_obs(obs_table)

   #Pandas is easier:
   obs_table = obs_table.to_pandas()
   data_products_by_obs = data_products_by_obs.to_pandas()

   # We are only interested in FLC and DRZ images
   data_products_by_obs = data_products_by_obs.loc[data_products_by_obs.project != 'HAP', :]
   obs_table = obs_table.merge(data_products_by_obs.loc[data_products_by_obs.productSubGroupDescription == 'FLC', :].groupby(['parent_obsid'])['parent_obsid'].count().rename_axis('obsid').rename('n_exp'), on = ['obsid'])

   obs_table['i_exptime'] = obs_table['t_exptime'] / obs_table['n_exp']

   #For convenience we add an extra column with the obstime
   obs_time = Time(obs_table['t_max'], format='mjd')
   obs_time.format = 'iso'
   obs_time.out_subfmt = 'date'
   obs_table['obs_time'] = obs_time
   obs_table['filters'] = obs_table['filters'].str.strip('; CLEAR2L CLEAR1L')

   data_products_by_obs = data_products_by_obs.merge(obs_table.loc[:, ['obsid', 'i_exptime', 'filters', 's_ra', 's_dec']].rename(columns={'obsid':'parent_obsid'}), on = ['parent_obsid'])

   #We select by individual exp time:
   obs_table = obs_table.loc[(obs_table.i_exptime > t_exptime_min) & (obs_table.i_exptime < t_exptime_max)]
   data_products_by_obs = data_products_by_obs.loc[(data_products_by_obs.i_exptime > t_exptime_min) & (data_products_by_obs.i_exptime < t_exptime_max)]

   #We add an ID
   obs_table['field_id'] = ['(%i)'%(ii+1) for ii in np.arange(len(obs_table))]

   return obs_table.astype({'obsid': 'int64'}).reset_index(drop = True), data_products_by_obs.astype({'parent_obsid': 'int64'}).reset_index(drop = True)


def search_data_products_by_obs(obs_table):
   """
   This routine search for images in MAST related to the given observations table.
   """

   data_products_by_obs = Observations.get_product_list(obs_table)

   return data_products_by_obs[((data_products_by_obs['productSubGroupDescription'] == 'FLC') | (data_products_by_obs['productSubGroupDescription'] == 'DRZ')) & (data_products_by_obs['obs_collection'] == 'HST')]


def download_HST_images(data_products_by_obs, path = './'):
   """
   This routine downloads the selected HST images from MAST.
   """

   try:
      images = Observations.download_products(Table.from_pandas(data_products_by_obs), download_dir=path)
   except:
      images = Observations.download_products(data_products_by_obs, download_dir=path)

   return images


def get_object_properties(args):
   """
   This routine will try to obtain all the required object properties from Simbad or from the user.
   """

   print('\n'+'-'*42)
   print("Commencing execution")
   print('-'*42)

   #Try to get object:
   if (args.ra is None) or (args.dec is None):
      try:
         from astroquery.simbad import Simbad
         import astropy.units as u
         from astropy.coordinates import SkyCoord

         customSimbad = Simbad()
         customSimbad.add_votable_fields('dim')

         object_table = customSimbad.query_object(args.name)

         object_name = str(object_table['MAIN_ID'][0]).replace("b'NAME ","'").replace("b' ","'")
         coo = SkyCoord(ra = object_table['RA'], dec = object_table['DEC'], unit=(u.hourangle, u.deg))

         args.ra = float(coo.ra.deg)
         args.dec = float(coo.dec.deg)

         #Try to get the search radius
         if all((args.search_radius == None, any((args.search_width == None, args.search_height == None)))):
            if (object_table['GALDIM_MAJAXIS'].mask == False):
               args.search_radius = max(np.round(float(2. * object_table['GALDIM_MAJAXIS'] / 60.), 2), 0.1)

      except:
         object_name = args.name
         if ((args.ra is None) or (args.dec is None)) and (args.quiet is False):
            print('\n')
            try:
               if (args.ra is None):
                  args.ra = float(input('R.A. not defined, please enter R.A. in degrees: '))
               if args.dec is None:
                  args.dec = float(input('Dec not defined, please enter Dec in degrees: '))
            except:
               print('No valid input. Float number required.')
               print('\nExiting now.\n')
               sys.exit(1)

         elif ((args.ra is None) or (args.dec is None)) and (args.quiet is True):
            print('HSTDownload could not find the object coordinates. Please check that the name of the object is written correctly. You can also run HSTDownload deffining explictly the coordinates using the "--ra" and "--dec" options.')
            sys.exit(1)

   else:
      object_name = args.name

   if (args.search_radius is None) and (args.quiet is False):
      print('\n')
      args.search_radius = float(input('Search radius not defined, please enter the search radius in degrees (Press enter to adopt the default value of 0.25 deg): ') or 0.25)
   elif (args.search_radius is None) and (args.quiet is True):
      args.search_radius = 0.25

   if (args.search_height is None):
      try:
         args.search_height = 2.*args.search_radius
      except:
         args.search_height = 1.0
   if (args.search_width is None):
      try:
         args.search_width = np.abs(2.*args.search_radius / np.cos(np.deg2rad(args.dec)))
      except:
         args.search_width = 1.0

   setattr(args, 'area', args.search_height * args.search_width * np.abs(np.cos(np.deg2rad(args.dec))))

   if args.hst_filters == ['any']:
      args.hst_filters = ['F555W','F606W','F625W','F775W','F814W','F850LP']

   name_coo = 'ra_%.3f_dec_%.3f_r_%.2f'%(args.ra, args.dec, args.search_radius)

   if args.name is not None:
      args.name = args.name.replace(" ", "_")
      args.base_file_name = args.name+'_'+name_coo
   else:
      args.name = name_coo
      args.base_file_name = name_coo


   #The script creates directories and set files names
   args.base_path = './%s/'%(args.name)
   args.HST_path = args.base_path+'HST/'

   args.used_HST_obs_table_filename = args.base_path + args.base_file_name+'_used_HST_images.csv'
   args.logfile = args.base_path + args.base_file_name+'.log'

   args.HST_obs_table_filename = args.HST_path + args.base_file_name+'_obs.csv'
   args.HST_data_table_products_filename = args.HST_path + args.base_file_name+'_data_products.csv'

   print('\n')
   print('-'*42)
   print('Search information')
   print('-'*42)
   print('- Object name:', object_name)
   print('- (ra, dec) = (%s, %s) deg.'%(round(args.ra, 5), round(args.dec, 5)))
   print('- Search radius = %s deg.'%args.search_radius)
   print('-'*42+'\n')

   return args


def str2bool(v):
   """
   This routine converts ascii input to boolean.
   """

   if v.lower() in ('yes', 'true', 't', 'y'):
      return True
   elif v.lower() in ('no', 'false', 'f', 'n'):
      return False
   else:
      raise argparse.ArgumentTypeError('Boolean value expected.')



def create_dir(path):
   """
   This routine creates directories.
   """

   if not os.path.isdir(path):
      try:
         tree = path.split('/')
         previous_tree = tree[0]
         for leave in tree[1:]:
            previous_tree = '%s/%s'%(previous_tree,leave)
            try:
               os.mkdir(previous_tree)
            except:
               pass
      except OSError:
         print ("Creation of the directory %s failed" % path)
      else:
         print ("Successfully created the directory %s " % path)


def hstdownload(argv):
   """
   Inputs
   """

   examples = '''Examples:

   hstdownload --name "Sculptor dSph"

   hstdownload --name "NGC 5053" --quiet

   hstdownload --ra 201.405 --dec -47.667 --search_radius 2 --hst_filters "F625W"
   '''

   parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, usage='%(prog)s [options]', description='HSTDownload downloads HST data using objetcs names or coordinates.', epilog=examples)

   # Search options
   parser.add_argument('--name', type=str, default = 'Output', help='Name for the Output table.')
   parser.add_argument('--ra', type=float, default = None, help='Central R.A.')
   parser.add_argument('--dec', type=float, default = None, help='Central Dec.')
   parser.add_argument('--search_radius', type=float, default = None, help='Radius of search in degrees.')
   parser.add_argument('--search_width', type=float, default = None, help='Width of the search rectangule in degrees.')
   parser.add_argument('--search_height', type=float, default = None, help='Height of the search rectangule in degrees.')
   parser.add_argument('--hst_filters', type=str, nargs='+', default = ['any'], help='Required filter for the HST images. Default all filters. They can be added as a list, e.g. "F814W" "F606W".')
   parser.add_argument('--hst_integration_time_min', type=float, default = 50, help='Required minimum average integration time for a set of HST images. This quantity is a limit on the average exposure time of an entire set of observations. Therefore, longer and shorter exposures may be available in an specific data set.')
   parser.add_argument('--hst_integration_time_max', type=float, default = 2000, help='Required maximum average integration time for a set of HST images. This quantity is a limit on the average exposure time of an entire set of observations. Therefore, longer and shorter exposures may be available in an specific data set. Exposures with less that 500 seconds of integration time are preferred. The default value is 2000 seconds, which is far more than the optimal value, but allow datasets with combinations of short and long exposures to be considered.')
   parser.add_argument('--project', type=str, nargs='+', default = ['HST'], help='Processing project. E.g. HST, HLA, EUVE, hlsp_legus. Default HST. They can be added as a list, e.g. "HST", "HLA".')
   parser.add_argument('--field_id', type=str, nargs='+', default = None, help='Specify the Ids of the fields to download. This is an internal id created by HSTDownload (field_id). The default value, "y", will download all the available HST observations fulfilling the required conditions. The user can also specify "n" for none, or the specific ids separated by spaces.')

   #Miscellaneus options
   parser.add_argument('--quiet', action='store_true', help='This flag deactivate the interactivity of HSTDownload. When used, HSTDownload will use all the default values without asking the user. This flag override and prevent the "--preselect_cmd" option.')

   if len(argv)==0:
      parser.print_help(sys.stderr)
      sys.exit(1)

   args = parser.parse_args(argv)
   args = get_object_properties(args)

   """
   The script creates directories and set files names
   """
   create_dir(args.base_path)
   create_dir(args.HST_path)


   """
   The script tries to load an existing HST table, otherwise it will download it from the MAST archive.
   """
   obs_table, data_products_by_obs = search_mast(args.ra, args.dec, search_width = args.search_width, search_height = args.search_height, filters = args.hst_filters, project = args.project, t_exptime_min = args.hst_integration_time_min, t_exptime_max = args.hst_integration_time_max)

   obs_table.to_csv(args.HST_obs_table_filename, index = False)
   data_products_by_obs.to_csv(args.HST_data_table_products_filename, index = False)

   if len(obs_table) > 0:

      """
      Ask whether the user wish to download the available HST images
      """

      print(obs_table.loc[:, ['obsid', 'filters', 'n_exp', 'i_exptime', 'obs_time', 'proposal_id', 's_ra', 's_dec', 'field_id']].to_string(index=False), '\n')

      if (args.quiet is True) and (args.field_id is None):
         print('HSTDownload will download the above sets of observations.\n')
         HST_obs_to_use = 'y'
      elif args.field_id is not None:
         print('HSTDownload will download the observations sets %s.'%(' '.join(str(p) for p in args.field_id) ))
         HST_obs_to_use = ' '.join([str(p) for p in args.field_id])
      else:
         print('Would you like to download the above HST observations?\n')
         print("Type 'y' or just press enter for all observations, 'n' for none. Type the id within parentheses at the right (field_id) if you wish to use that specific set of observations. You can enter several ids separated by space. \n")
         HST_obs_to_use = input('Please type your answer and press enter: ') or 'y'
         print('\n')

      try:
         HST_obs_to_use = str2bool(HST_obs_to_use)
      except:
         try:
            HST_obs_to_use = list(set([obsid for obsid in obs_table.obsid[[int(obsid)-1 for obsid in HST_obs_to_use.split()]] if np.isfinite(obsid)]))
         except:
            print('No valid input. Not downloading observations.')
            HST_obs_to_use = False

      if HST_obs_to_use is not False:
         if HST_obs_to_use is True:
            HST_obs_to_use = list(obs_table['obsid'].values)
         hst_images = download_HST_images(data_products_by_obs.loc[data_products_by_obs['parent_obsid'].isin(HST_obs_to_use), :], path = args.HST_path)
      else:
         print('\nExiting now.\n')
         sys.exit(1)

      """
      Select only flc
      """
      drz_images = data_products_by_obs[(data_products_by_obs['productSubGroupDescription'] == 'DRZ') & (data_products_by_obs['parent_obsid'].isin(HST_obs_to_use))]
      flc_images = data_products_by_obs[(data_products_by_obs['productSubGroupDescription'] == 'FLC') & (data_products_by_obs['parent_obsid'].isin(HST_obs_to_use))]


      """
      Save Gaia and HST tables
      """
      obs_table.to_csv(args.HST_obs_table_filename, index = False)
      data_products_by_obs.to_csv(args.HST_data_table_products_filename, index = False)
      flc_images.to_csv(args.used_HST_obs_table_filename, index = False)

   else:
      if args.quiet:
         print('No suitable HST observations were found. Please try with different parameters. Exiting now.')
      else:
         input('No suitable HST observations were found. Please try with different parameters.\nPress enter to exit.\n')

if __name__ == '__main__':
    hstdownload(sys.argv[1:])
    sys.exit(0)

"""
Andres del Pino Molina
"""
