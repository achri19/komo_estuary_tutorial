#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 07:28:56 2021

@author: Alchrist
This script is specifically for developing new SWORD channel make_SWORD_networks
It used a 10m water mask build in GEE from optical Sentinel-2 and/or radar Sentinel-1 images
It can be run locally or on Gattaca, but paths should match those below

"""

import sys
import os
import pandas as pd
import shutil
from datetime import datetime
from string import Template
import fnmatch
import geopandas as gpd
import platform
import psutil
import rasterio
from osgeo import gdal
from pathlib import Path
import numpy as np

if os.getcwd().split('/')[1] == 'Users':
    path = '/Users/alchrist/Documents/GitHub/ANUGA/'
    path_code = path + 'processing/code/'
    path_templates = path + 'processing/templates/'
    path_configs = path  + 'processing/configs/'
    cnes_tools_path = path + 'swot-hydroloy-toolbox/'
    path_ancillary = path + 'ancillary/'
    path_examples = '/Volumes/FortressL3/ANUGA/GlobalDeltas/examples/' #''/Volumes/FortressL3/ANUGA/GlobalDeltas/examples/'
else:
    path = '/projects/loac_hydro/alchrist/anuga/'
    path_code = path + 'code/'
    path_templates = path + 'templates/'
    path_configs = path + 'configs/'
    path_ancillary = path + 'ancillary/'
    cnes_tools_path = '/scratch_sm/loac_hydro/swot-hydrology-toolbox/'
    path_examples = path + '/examples/'

sys.path.insert(1,path_code)
print(path_code)
""" These functions were developed for the BYOM tutorials, but are used for developing """

from BYOM_Utilities_V1 import (build_directory,
                               get_extent_parameters,
                               make_polygons,
                               make_channel_networks,
                               make_watermask,
                               more_opening)

print(f"Computer network name: {platform.node()}") #Computer network name
print(f"Machine type: {platform.machine()}") #Machine type
print(f"Processor type: {platform.processor()}") #Processor type
print(f"Platform type: {platform.platform()}") #Platform type
print(f"Operating system: {platform.system()}") #Operating system
print(f"Operating system release: {platform.release()}") #Operating system release
print(f"Operating system version: {platform.version()}") #Operating system version
print(f"Number of physical cores: {psutil.cpu_count(logical=False)}") #Physical cores
print(f"Number of logical cores: {psutil.cpu_count(logical=True)}") #Logical cores
print(f"Current CPU utilization: {psutil.cpu_percent(interval=1)}") #System-wide CPU utilization
print(f"Current per-CPU utilization: {psutil.cpu_percent(interval=1, percpu=True)}") #System-wide per-CPU utilization
print(f"Total RAM installed: {round(psutil.virtual_memory().total/1000000000, 2)} GB") #Total RAM
print(f"Available RAM: {round(psutil.virtual_memory().available/1000000000, 2)} GB") #Available RAM
print(f"Used RAM: {round(psutil.virtual_memory().used/1000000000, 2)} GB") #Used RAM
print(f"RAM usage: {psutil.virtual_memory().percent}%") #RAM usage

##############################################################################
######################### Get Configuration Parameters #######################
##############################################################################
# Resolution all output files - ANUGA models are always in UTM
res = 10

configfile = '//Users/Alchrist/Documents/Github/globalcoast/configs/Deltas_SetupParameters_MASTER_V2.csv' # #'Deltas_SetupParameters_MASTER_V2.csv'path_configs + 'need_opening.csv'#'
parameter_file = pd.read_csv(configfile)
AOIs = parameter_file['AOI']

i=0
#AOIs = ['guayas']
for AOI in AOIs[:]:
    skip = False
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ANUGA Model Setup ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')
    print('################ ################# ##################')

    print('Before building the model, we need 3 pieces: \n1) model domain (shapefile) \n2) approximate tide boundary (shapefile) \n3) watermask')

    AOI = AOI.lower()
    parameters = parameter_file[parameter_file['AOI']==AOI.capitalize()].reset_index(drop=True)

    ## Set up working directory and build sub folders if needed
    print('\n')
    print('Study area is ' + AOI)
    print('Resolution of this setup is %sm' %(res))
    working_path,folders = build_directory(path_examples, AOI)

    final_watermasks = [os.path.join(dirpath,f)
        for dirpath,dirnames, files in os.walk(folders[0])
        for f in fnmatch.filter(files,'*_finalwatermask*.tif')]

    if len(final_watermasks) >0:
        print('You already have a chosen water mask:')
        print(final_watermasks)

        ref_10m,parameters = get_extent_parameters(path_ancillary,AOI,folders,res,parameters)

        print(parameters.iloc[0] )

        EPSG = parameters['EPSG'][0]
        ulx = parameters['ulx'][0]
        uly = parameters['uly'][0]
        lrx = parameters['lrx'][0]
        lry = parameters['lry'][0]

        model_domain = gpd.read_file('%s%s_modeldomain.shp' %(folders[7],AOI))
        AOI_extent = gpd.read_file('%s%s_extent_%s.shp' %(folders[7],AOI,EPSG))

        clean_with_landcover = False

        print('\n Build a reference raster (GEBCO) that is used to set CRS, extent, resolution to make sure all future files align .......\n')

        watermaskname = make_watermask(path_ancillary,
                                       AOI,
                                       folders,
                                       parameters,
                                       ref_10m,clean_with_landcover,True)
#                                       os.path.isfile('%s%s_landmask_%s.tif' %(folders[8],AOI,10)))

        how_much_opening = 3
        #more_opening(AOI,folders,watermaskname,how_much_opening,ref_10m,parameters)

        watermask = rasterio.open('%s%s_watermask_%s.tif' %(folders[8],AOI,res)).read(1)


        make_polygons(AOI,
                        folders,
                        parameters,
                        ref_10m,
                        watermaskname,
                        path_templates,True)
#                        os.path.isfile("%s%s_rivers_%sb.shp" %(folders[7],AOI,res)))


        segment_width = 20
        pixel_step = int(round(segment_width/res))
        segment_width_2 = 30
        pixel_step_2 = int(round(segment_width_2/res))
        segment_width_3 = 40
        pixel_step_3 = int(round(segment_width_3/res))
        try:
            distance,widths = make_channel_networks(folders,
                                                  AOI,
                                                  ref_10m,
                                                  parameters,
                                                  pixel_step,False)
#                                                  (os.path.isfile("%s%s_widths_%sx%sb.tif" %(folders[8],AOI,res,pixel_step))) | (os.path.isfile("%s%s_widths_%sx%s.tif" %(folders[8],AOI,res,pixel_step_2))))
        except:
            try:
                distance,widths = make_channel_networks(folders,
                                                  AOI,
                                                  ref_10m,
                                                  parameters,
                                                  pixel_step_2,False)
#                                                  os.path.isfile("%s%s_widths_%sx%sb.tif" %(folders[8],AOI,res,pixel_step_2)))
            except:
                distance,widths = make_channel_networks(folders,
                                                  AOI,
                                                  ref_10m,
                                                  parameters,
                                                  pixel_step_3,False)
#                                                  os.path.isfile("%s%s_widths_%sx%sb.tif" %(folders[8],AOI,res,pixel_step_3)))

        print('Cleaning up temporary files')
        try:shutil.rmtree(folders[1])
        except:''
        i = i+1
