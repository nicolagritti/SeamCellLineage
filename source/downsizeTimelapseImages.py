# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 15:36:52 2015

@author: gritti
"""

import os
import shutil
import glob
import numpy as np
import time
from timelapseFun import *
from skimage import io
from PIL import Image

####################################################################################################
# Give the path and the worm name to be converted
####################################################################################################

path = 'X:\\Michael\\160323_MBA21_SRF3'

# folders = [ i for i in os.listdir(path) if os.path.isdir(path+i) ]
# foldersToKeep = [ ( not i.endswith('downsized') ) & ( not i.endswith('straighten') ) for i in folders ]
# worms = [ i for idx, i in enumerate(folders) if foldersToKeep[idx] ]
# worms.sort()

worms = [ x for x in os.listdir(path) if len(x) == 3 ]
worms = [ 'C29', 'C30' ]

scaleFactor = 8

if type( worms ) == str:
    worms = list(worms)

def downsize():
    for worm in worms:
        
        print('Downsizing worm '+worm)

        ################################################################################################
        # Load the list of the images and metadata files
        ################################################################################################
        
        # define the input folder and the output folder
        inpath = os.path.join( path, worm )
        outpath = os.path.join( path, worm+'_downsized' )
        
        # read in all the images to be converted
        channels = [False,False,False]
        flist = [None,None,None]
        
        if os.path.isfile( os.path.join( inpath, 'z001_488nm.tif' ) ):
            flist[0] = glob.glob( os.path.join( inpath, 'z*488nm.tif' ) )
            flist[0].sort()
            channels[0] = True
        if os.path.isfile( os.path.join( inpath, 'z001_561nm.tif' ) ):
            flist[1] = glob.glob( os.path.join( inpath, 'z*561nm.tif') )
            flist[1].sort()
            channels[1] = True
        if os.path.isfile( os.path.join( inpath, 'z001_CoolLED.tif' ) ):
            flist[2] = glob.glob( os.path.join( inpath, 'z*CoolLED.tif') )
            flist[2].sort()
            channels[2] = True
        
        metalist = glob.glob( os.path.join( inpath, 'z*.txt' ) )
        metalist.sort()
        
        ################################################################################################
        # Create the directory and copy metadata files
        ################################################################################################
        
        # create the directory
        if not os.path.isdir(outpath):
            os.mkdir(outpath)
        
        # copy the metadataFiles
        for f in metalist:
            # print('copying ' + f)
            if not os.path.isfile(outpath+'\\'+f.split('\\')[-1]):
                shutil.copyfile(f, outpath+'\\'+f.split('\\')[-1])
        
        ################################################################################################
        # Downsize and save the images
        ################################################################################################
        
        # for each channel, if True
        for jdx, chn in enumerate( channels ):
        
            if chn:
                
                # for each image in that channel
                for idx in np.arange(len(metalist)):
                    
                    
                    f = flist[jdx][idx]
                    if not os.path.isfile(outpath+'\\'+f.split('\\')[-1]):
                
                        print('loading ' + flist[jdx][idx])

                        start = time.time()

                        imgs = Image.open( f, mode = 'r' )
                        print(imgs.seek(0))
                        
                        end1 = time.time()
                        print( end1-start )

                        # smallimgs = []
                        # for img in imgs:
                        #     Nbig = img.shape[0]
                        #     Nsmall = img.shape[0]/scaleFactor
                        #     smallimg = ( img.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1) ).astype(np.uint16)
                        #     smallimgs.append( smallimg )
                        # smallimgs = np.array(smallimgs)

                        # end2 = time.time()
                        # print(end2-end1)
                    
                        # imsave(outpath+'\\'+f.split('\\')[-1],smallimgs)

                        # print( time.time() - end2 )
        
    print('Waiting for one hour to restart!')
    time.sleep(60*60)

while True:
    downsize()