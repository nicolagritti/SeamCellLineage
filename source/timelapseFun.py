# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 08:57:15 2015

@author: gritti

NB: this file is also in Miniconda3/Lib/site-packages. This way the file doesn't need to be in the script folder.
I have a copy here just as a backup!

"""

import glob
from tifffile import *
import numpy as np
#import PIL.Image as Image
import os.path
import matplotlib.pyplot as plt
import pickle
import scipy.interpolate as ip
import re
from datetime import datetime
from skimage.filter import threshold_otsu, rank
from skimage import measure
from skimage.morphology import remove_small_objects, disk
from scipy.ndimage import morphology, filters
import pandas as pd

# from uiget_dir import *
#import psutil




##############################################################################
# LOADING STACK, IMSHOW AND PLOTTING FUNCTIONS
##############################################################################

def loadstack(filename):

    '''
    This function reads a (multi)tif file and returns the images in a numpy 
    array using tifffile package from Gohlke.
    '''

    with TiffFile(filename) as f:
        return f.asarray()



def change_space_time(increment, st, upperLimit, lowerLimit=0):

    '''
    This function changes a parameters (could be the slice or the timepoint)
    by an increment. The new value is clipped between 0 and an upper limit 
    (often being the number of slices or the number of timepoints).
    '''

    return np.clip(st + increment, lowerLimit, upperLimit)



def show_img_outline(fig, tp, sl, img, outline, spline, movetype='', clim=np.array([1., 1.]), lethState = ''):

    '''
    Shows on the fig the image with the outline and the spline.
    '''

    fig.clf()
    ax = fig.add_axes([-0., -0., 1., 1.])

    imgplot = ax.imshow( img, cmap='gray' )
    imgplot.set_clim( np.min(img) * clim[0], np.max(img) * clim[1] )
    
#    if outline.shape[1] != 0:
    ax.plot( outline[:,0], outline[:,1], '-o', color='red', ms=8, mew=1.5, alpha=.5, lw = 2 )
    ax.plot( spline[:,0], spline[:,1], '-', color='green', lw = 2 )

    if lethState != '':
        ax.text( 20, 40, lethState, color = 'red' )

    # legend of cells in the previous timepoint
    fig.canvas.set_window_title('Timepoint: %i,' % tp + ' Slice: %i,' % sl +
                                ' Mouse wheel function: %s' % movetype)

    plt.axis('tight')
    plt.axis('equal')
    plt.axis('off')
    plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)



def show_img_cells(fig, tp, sl, img, cells, movetype=None, clim=np.array([1., 1.]), i=1):

    '''    
    Shows on the fig the image with the labelled cells.    
    '''

    ax = fig.add_axes([0., .2 * i, 1., 1.])

    imgplot = ax.imshow(img, cmap='gray')
    imgplot.set_clim((np.min(img) + np.std(img) / 5.) * clim[0], np.max(img) * clim[1])
    ax.text(0, -30, 'Timepoint: %i,' % tp + ' Slice: %i,' % sl, color='blue', size='large')

    for idx, cell in cells.iterrows():

        if cell.cZpos == sl:
            clabel = str(cell.cname) + ' ' + str(cell.cside)
            ax.text( cell.cXpos, cell.cYpos + 80, clabel, color='red', size='large', alpha=.8,
                    rotation=90)

    if movetype:
        fig.canvas.set_window_title(' Mouse wheel function: %s' % movetype +
                                    '   NB: f1=save, i (u) = invert side of cells (R bottom, L top)')
    plt.axis('off')



def show_img_cells_shape(fig, tp, sl, img, cells, edges, movetype=None, clim=np.array([1., 1.]) ):

    '''  
    Shows on the fig the image with the labelled cells.  
    '''

    ax = fig.add_axes([-.225, -.225, 1.45, 1.45])

    imgplot = ax.imshow(img, cmap='gray')
    imgplot.set_clim((np.min(img) + np.std(img) / 5.) * clim[0], np.max(img) * clim[1])
    ax.text(0, -30, 'Timepoint: %i,' % tp + ' Slice: %i,' % sl, color='blue', size='large')

    for c in cells:

        if c[1] == sl:
            ax.text(c[2][0], c[2][1] + 80, str(c[0]) + ' ' + str(c[-1]), color='red', size='large', alpha=.8,
                    rotation=90)
    
    for pos in edges:
        if pos[-1] == sl:
            ax.plot(pos[0], pos[1], 'yo', linewidth = 10)
            

    if movetype:
        fig.canvas.set_window_title(' Mouse wheel function: %s' % movetype +
                                    '   NB: f1=save, i (u) = invert side of cells (R bottom, L top)')
    plt.axis('tight')
    plt.axis('equal')
    plt.axis('off')
    plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)




##############################################################################
# UTILS
##############################################################################

def ortho(v):
    
    '''
    Returns the orthogonal vector as numpy array.
    '''

    return np.array([-v[1], v[0]])



def norm(v):
    
    '''
    Returns the norm of the vector as a float.    
    '''

    return np.sqrt(v[0] ** 2 + v[1] ** 2)



def find_closest_point(ref,pos):
    
    '''
    Returns the index of the closest position in pos to ref
    '''
    
    dist = []
    for p in pos:
        dist.append( (ref[0]-p[0])**2 + (ref[1]-p[1])**2 )
    return dist.index(np.min(dist))



def create_worm( path, wormnumber, magnification, n, firsttp):

    worm = {
        'rowtype': np.nan,
        'wormname': np.nan,
        'magnification': np.nan,
        'tidx': np.nan,
        'times': np.nan,
        'outline': np.nan,
        'spline': np.nan,
        'length': np.nan,
        'cname': np.nan,
        'cXpos': np.nan,
        'cYpos': np.nan,
        'cZpos': np.nan,
        'cside': np.nan,
        'cedges': np.nan,
        'csize': np.nan
        }
        
    df = pd.DataFrame( worm, index = [0] )

    # print(path+wormnumber+'_downsized\\*.tif')
    flist = glob.glob(path+wormnumber+'_downsized\\*.tif')
    size = 2048. / loadstack(flist[0])[0].shape[0]

    # insert parameters
    df.ix[0,'rowtype'] = 'param'
    df.ix[0,'wormname'] = path+'\\'+wormnumber
    df.ix[0,'magnification'] = magnification
    df.ix[0,'compression'] = 2048. / size
    
    # insert body data
    times = extract_times(path+'\\'+wormnumber,firsttp)[firsttp-1:]
    df1 = pd.DataFrame( { 'rowtype': 'body',
                          'tidx': np.arange(len(times)),
                          'times': times } )
    
    df = pd.concat( [ df, df1 ] ).reset_index(drop=True)
        
    df.cedges = df.cedges.astype(object)
    df.cname = df.cname.astype(object)
    df.cside = df.cside.astype(object)
    df.outline = df.outline.astype(object)
    df.spline = df.spline.astype(object)
    for i in np.arange(1,len(df.outline.values)):    
        df.outline.values[i] = np.array([[np.nan,np.nan,np.nan]])
        df.spline.values[i] = np.array([[np.nan,np.nan,np.nan]])
    df.tidx = df.tidx.astype(float)
    df.magnification = df.magnification.astype(float)
    
    pickle.dump( df, open(path+'\\worm'+wormnumber+'.pickle','wb'), protocol=2 )




def extract_times(path, zero=1):

    '''
    Reads all the metadata files in the specified folder and returns the 
    timesteps (in minutes) from an xero timepoint. Times are returned as a list
    '''
    zero = int(zero)
    flist = glob.glob(path + '\\z*.txt')
    flist.sort()
    # print(flist,path)
    times = []
    ftimezero = flist[zero-1]

    with open(ftimezero, 'r') as f:

        line = ''
        while 'Date/Time' not in line:
            line = f.readline()

        timezero = datetime.strptime(line.strip().split(': ')[1:][0], '%Y-%m-%d %H:%M:%S')

    for fname in flist:

        with open(fname, 'r') as f:

            line = ''
            while 'Date/Time' not in line:
                line = f.readline()

            date_time = datetime.strptime(line.strip().split(': ')[1:][0], '%Y-%m-%d %H:%M:%S')

            times.append(( date_time - timezero ).total_seconds() / 60. / 60.)

    return times



def closer_cell(pos, cdf, ftp):

    '''
    Returns the index and the cell in the cell dataframe of the closer cell to the position pos.
    Very useful when labeling cells. In that case, [x,y] is the position of the mouse on the figure.    
    '''
    
    tmppos = np.copy( pos )
    # read the z-step from txt file
    npxlz = 0
    with open(ftp, 'r') as f:

        line = ''
        while 'From' not in line:
            line = f.readline()
        npxlz = -float(line.strip().split(' ')[1])

        line = ''
        while 'To' not in line:
            line = f.readline()
        npxlz += float(line.strip().split(' ')[1])

        line = ''
        while 'Steps' not in line:
            line = f.readline()
        npxlz = npxlz / ( float(line.strip().split(' ')[1][:-1]) * 0.1 )

    c = np.nan
    idx = np.nan
    dist = 100000
    tmppos[-1] *= npxlz

    for j, cell in cdf.iterrows():
        cpos = np.array( cell[['cXpos','cYpos','cZpos']].values )
        cpos[-1] *= npxlz

        if dist > np.sqrt( np.sum( (cpos-tmppos)**2 ) ):
            dist = np.sqrt( np.sum( (cpos-tmppos)**2 ) )
            c = cell
            idx = j

    return (idx, c)






def automatic_outline_michalis( imgs, o, sl, tp, scalefactor = 1 ):
    
    '''
    outline[0] = head
    outline[1] = tail
    '''
    
    minsl = 0#np.clip( np.min( [ outline[0][-1], outline[-1][-1] ] ), 0, len(imgs) )
    maxsl = 19#np.clip( np.max( [ outline[0][-1], outline[-1][-1] ] ), 0, len(imgs) )
    
    # compute the max over the stack
    img = np.max( imgs[minsl:maxsl+1], 0 )

    # resize the image
    Nbig = img.shape[0]
    Nsmall = img.shape[0]/scalefactor
    smallimg = img.reshape([Nsmall, Nbig/Nsmall, Nsmall, Nbig/Nsmall]).mean(3).mean(1)
    
    # compute local variance with generic filter function
    varimg = filters.generic_filter( smallimg, np.var, size=2 ).astype('uint16')
    
    thr = threshold_otsu( varimg )    
    labels = measure.label( morphology.binary_dilation( varimg>thr, disk(2) ) )
    labels = remove_small_objects( labels , 12 )

    pos = [ [ scalefactor*i[ 'centroid' ][1], scalefactor*i[ 'centroid' ][0] ] for i in measure.regionprops(labels) ]

    sortedpos = np.array([o[0]])
    for i in np.arange(len(pos)):
#            print(len(pos))
        idx = find_closest_point(sortedpos[-1],pos)
        newsl = minsl + list(imgs[minsl:,pos[idx][1],pos[idx][0]]).index(np.max(imgs[minsl:,pos[idx][1],pos[idx][0]]))
        pos[idx].append(newsl)
        sortedpos = np.append( sortedpos, np.array( [ pos[idx] ] ), axis=0 )
        pos.pop(idx)
    sortedpos = np.append( sortedpos, np.array( [o[-1]] ), axis=0 )
    
    return sortedpos



def update_lethargus(tp, oldleth, event):
    
    '''
    updates the lethargus stages
    '''
    
    stages = ['L1','L2','L3','L4','egg']
    newleth = list(oldleth)
    
    if event.key == 'f3':  # add a developmental timing
    
        if newleth[tp-1] != '':
            
            newleth[tp] = newleth[tp-1]
            
        elif all( [ i=='' for i in newleth[:tp] ] ):
            
            newleth[tp] = 'L1'
            
        else:
            
            lastLeth = newleth[ len(newleth[:tp]) - [i!='' for i in newleth[:tp]][::-1].index(True) - 1 ]
            idxNewLeth = stages.index( lastLeth ) + 1
            newleth[tp] = stages[ idxNewLeth ]
            
    else:
        newleth[tp] = ''
    
    if newleth[tp] == 'egg':

        newleth[tp:] = [ 'egg' for i in newleth[tp:] ]

    return newleth




def interp_spline(outline):

    '''
    From the outline, this function performs a 1d spline interpolation with evenly spaced points (ds = 1 pixel).
    It first computes a 1d cubic interpolation. This gives rise to not equally spaced points.
    That's why a linear interpolation is than needed to compute the final spline curve.    
    '''

    spline = [[], [], []]

    outline2d = outline[:,:2]  # remove the z component

    if len(outline2d) < 4:
        spline[0] = []
        spline[1] = []

    if len(outline2d) > 3:
        L3d = np.array(outline).T
        L = L3d[:-1]  # remove the z component

        # spline with cubic interpolation
        dist = np.cumsum(np.insert(np.sqrt(np.sum(np.diff(L) ** 2, 0)), 0, 0))
        dnew = np.linspace(0, dist[-1], dist[-1])
        sf = ip.interp1d(dist, L, kind='cubic')
        splinetmp = sf(dnew)

        # linear spline for evenly spaced points
        splinedist = np.cumsum(np.insert(np.sqrt(np.sum(np.diff(splinetmp) ** 2, 0)), 0, 0))
        splinednew = np.linspace(0, splinedist[-1], splinedist[-1])

        spline[0] = np.interp(splinednew, splinedist, splinetmp[0])
        spline[1] = np.interp(splinednew, splinedist, splinetmp[1])

        # linear interpolation of the slices
        spline[2] = np.interp(splinednew, dist, L3d[-1])

    return np.array( spline ).T


def straighten(path, flist, spline3d, width=100):
    
    '''
    Takes the images in the specified folder. Using the spline curve, 
    it performs a 2d interpolation and creates a strighten z-stack 
    image. The stack is than saved in a different	multitif file (called 'straight%time_488nm.tif') in the same folder.
    '''

    if not os.path.exists(path + '_straighten'):
        os.mkdir(path+'_straighten')

    for idx in np.arange(len(flist)):
        time = idx+1
        channel = flist[idx].split('_')[-1]
#        print(flist[time],channel)

        spline = interp_spline(spline3d[idx])[:,:-1]
        print( flist[idx] )

        if len(spline[:,0]) > 1 and not os.path.isfile(path + '_straighten\\straight%.3d_%s'% (time, channel)):

            imgs = loadstack(flist[idx])
            straightenimgs = np.zeros(( len(imgs), 2 * width + 1, len(spline[:,0]) - 2 ))

            # first compute orthogonal vectors and coordinates
           # print('\t computing coordinates tp: %.3d' % time)
            coord = np.zeros(( len(spline[:,0]) - 2, 2 * width + 1, 2 ))

            for i, pos in enumerate(spline[1:-1,:]):
                # for each position on the spline, compute tangent and orthogonal vector (normalized)
                tang = spline[i] - spline[i+2]#np.array([spline[i,0] - spline[i+2,0], spline[i,1] - spline[i+2,1]])
                vort = ortho(tang) / norm(tang)

                # build the line with the coordinates for the interpolation
                coord[i, :, :] = np.array([np.arange(-width, width + 1), np.arange(-width, width + 1)]).T * vort + spline[i + 1,:]
                                 

            # then compute, for each slice, the interp2d and extrapolate the value for each pixel in the coord array
            for s, img in enumerate(imgs):

                print('tp: %.3d\t\t straightening slice%3d' % (time,s))

                # interpolation function of the image
                xo = np.arange(len(img[0]))
                yo = np.arange(len(img[0]))
                f = ip.interp2d(xo, yo, img, kind='cubic')

                # interpolation of each vertical line of the straighten image
                for i, c in enumerate(coord):
                    # this takes into account the case when the coords are decreasing: interp2d can only interpolate
                    # on a grid with sorted values (increasing)
                    sx = 2 * ( c[-1, 0] > c[0, 0] ) - 1
                    sy = 2 * ( c[-1, 1] > c[0, 1] ) - 1

                    # interpolate the line: interpolate the grid and take the appropriate diagonal
                    straightenimgs[s, :, i] = np.clip(
                        np.diagonal(f(np.sort(c[:, 0]), np.sort(c[:, 1]))[::sy, ::sx]).astype('uint16'),
                         0, 2 ** 16 - 1)

            imsave(path + '_straighten\\straight%.3d_%s'% (time, channel), straightenimgs.astype('uint16'))

def straighten_1(path, flist, spline3d, width=100):
    
    '''
    Takes the images in the specified folder. Using the spline curve, 
    it performs a 2d interpolation and creates a strighten z-stack 
    image. The stack is than saved in a different	multitif file (called 'straight%time_488nm.tif') in the same folder.
    '''

    if not os.path.exists(path + '_straighten'):
        os.mkdir(path+'_straighten')

    for idx in np.arange(len(flist)):
        time = idx+1
        channel = flist[idx].split('_')[-1]
#        print(flist[time],channel)

        spline = interp_spline(spline3d[idx])[:,:-1]
        print( flist[idx] )

        if len(spline[:,0]) > 1 and not os.path.isfile(path + '_straighten\\straight%.3d_%s'% (time, channel)):

            imgs = loadstack(flist[idx])
            straightenimgs = np.zeros(( len(imgs), 2 * width + 1, len(spline[:,0]) - 2 ))

            # first compute orthogonal vectors and coordinates
            coord = np.zeros(( len(spline[:,0]) - 2, 2 * width + 1, 2 ))

            for i, pos in enumerate(spline[1:-1,:]):
                # for each position on the spline, compute tangent and orthogonal vector (normalized)
                tang = spline[i] - spline[i+2]#np.array([spline[i,0] - spline[i+2,0], spline[i,1] - spline[i+2,1]])
                vort = ortho(tang) / norm(tang)

                # build the line with the coordinates for the interpolation
                coord[i, :, :] = np.array([np.arange(-width, width + 1), np.arange(-width, width + 1)]).T * vort + spline[i + 1,:]
                
            # creates the 3d points on which to evaluate the interpolation function
            coords = np.array( [ 
                    [ np.vstack((np.array([j for i in c]),c[:,1],c[:,0])).T for c in coord ] 
                    for j in np.arange(len(imgs)) ] )

            # define the grid and the interpolation function
            xo = np.arange(len(imgs[0][0]))
            yo = np.arange(len(imgs[0][0]))
            zo = np.arange(len(imgs))
            func = ip.RegularGridInterpolator( (zo, yo, xo), imgs )
            
            # evaluate and rotate the images horizontally
            straightenimgs = func(coords).astype('uint16')
            straightenimgs = np.transpose(straightenimgs,(0,2,1))
            
            # save the images
            imsave(path + '_straighten\\straight%.3d_%s'% (time, channel), straightenimgs.astype('uint16'))


def resize_straighten_worms(path):
    
    '''
    This function takes all the straighten images in the specified path, calculates the larger image and replace all the smaller
    one with a black band in the right side. The idea is to make all the images with the same size.
    '''

    flist = glob.glob(path + '\\straight*.tif')

    size = []
    for filename in flist:
        img = loadstack(filename)[0]
        size.append([img.shape[0], img.shape[1]])

    size = np.max(np.array(size).T, 1)

    for filename in flist:

        imgs = loadstack(filename)
        newimgs = np.zeros(( len(imgs), size[0], size[1] ))

        for iimg, img in enumerate(imgs):
            newimgs[iimg][:img.shape[0], :img.shape[1]] = img

        newimgs = np.array(newimgs, dtype='uint16')
        imsave(filename, newimgs)


def maximum_projection(imgs, lower_slice=0, upper_slice=None, leftlim=0, rightlim=None):
    
    '''
    This function performs a maximum projection of a Z-stack between two defined slices
    '''

    if not rightlim:
        rightlim = imgs.shape[2] + 1

    if not upper_slice:
        upper_slice = imgs.shape[0] + 1

    return np.max(imgs[lower_slice: upper_slice + 1, :, leftlim: rightlim], 0)


def build_worm_movie(path, filetemplate='\\straight*.tif'):
    
    '''
    Function to create a movie of a (straighten) worm. 
    NB: If you want a movie of a full microchamber pass also the filetemplate argument.
    '''

    filename = path + filetemplate
    flist = glob.glob(filename)

    cells = pickle.load(open(path + '\\cells.pickle', 'rb'))

    movie = [[] for i in flist]
    bc = [0, 10000]

    regexp = re.compile('([0-9]{3})_[0-9]{3}')

    for idx, f in enumerate(flist):

        imgs = loadstack(f)

        if len(cells[idx]) != 0:
            slices = [c[1] for c in cells[idx]]

            l = np.clip(np.min(slices), 0, len(imgs))
            u = np.clip(np.max(slices), 0, len(imgs))
        else:
            l = 0
            u = len(imgs) + 1

        movie[idx] = maximum_projection(imgs, l, u)
        movie[idx] = (movie[idx] - np.min(movie[idx])) / (np.max(movie[idx]) - np.min(movie[idx])) * (bc[1] - bc[0])

    movie = np.array(movie, dtype='uint16')

    if not os.path.isfile(path + '\\movie.tif'):
        imsave(path + '\\movie.tif', movie)

    return movie


def build_cell_movie(path, cellname):

    '''
    from a worm and a seam cell name, i.e. '1', this function builds the movie of the cell lineage and divisions.
    NB: use this function only for (V1-V5) seam cells! This functions doesn't work with T and H cells because it needs to
    calculate the distance to the closest left and right cell. Anyway it's only for making a fancy movie, so it doesn't really
    matter which cell you choose. And the V cells are the nicest to make a movie!
    '''

    filename = path + '\\straight*_488nm.tif'
    flist = glob.glob(filename)

    cells = pickle.load(open(path + '\\cells.pickle', 'rb'))

    # compute distance to the closest left and right cell
    centralpos = [[] for c in cells]
    rightspace = []
    leftspace = []

    for idx, cs in enumerate(cells):
        pos = []
        [pos.append(c[2]) for c in cs if c[0][0] == cellname]
        centralpos[idx] = np.mean(pos, 0)
        '''
        TO BE FIXED!!!
        '''

        rightlim = cs[idx_closer_cell(centralpos[idx], clear_cell_list(cs, cellname, centralpos[idx][0], +1))][2]
        leftlim = cs[idx_closer_cell(centralpos[idx], clear_cell_list(cs, cellname, centralpos[idx][0], -1))][2]

        rightspace.append(rightlim[0] - centralpos[idx][0])
        leftspace.append(centralpos[idx][0] - leftlim[0])

    maxrightspace = int(np.max(rightspace))
    maxleftspace = int(np.max(leftspace))

    # center the images
    cellmovie = [[] for i in cells]

    for idx, cs in enumerate(cells):
        slices = []
        [slices.append(c[1]) for c in cs if c[0][0] == cellname]

        imgs = loadstack(flist[idx])
        l = np.clip(np.min(slices), 0, len(imgs))
        u = np.clip(np.max(slices), 0, len(imgs))

        cellmovie[idx] = maximum_projection(imgs, l, u, centralpos[idx][0] - maxleftspace,
                                            centralpos[idx][0] + maxrightspace)

    # adjust b&c
    bc = [0, 10000]

    for idx, celltp in enumerate(cellmovie):
        cellmovie[idx] = (celltp - np.min(celltp)) / (np.max(celltp) - np.min(celltp)) * (bc[1] - bc[0])

    # save the movie
    cellmovie = np.array(cellmovie, dtype='uint16')

    if not os.path.isfile(path + '\\movie_' + cellname + '.tif'):
        imsave(path + '\\movie_' + cellname + '.tif', cellmovie)

    return cellmovie


def create_cell( refpos, bd, side = None, cname = None ):

    if not side:
        side = insert_side( refpos, bd )
    if not cname:
        cname = np.nan
        
    newcell = pd.DataFrame( {
        'rowtype': 'cell',
        'cname': cname,
        'cXpos': refpos[0],
        'cYpos': refpos[1],
        'cZpos': refpos[2],
        'cside': side,
        'tidx': bd.tidx.values[0],
        'times': bd.times.values[0]}, index=[0] )
        
    return newcell





def insert_cell_division( cellstp, ctemplate, refpos, ftp, bd ):
    
    # find the closest cell in the current timepoint
    idx, cell = closer_cell( refpos, cellstp, ftp )
    
    # update the name of the existing cell in ctemplate
    namemask = ctemplate['cname'] == cell.cname
    sidemask = ctemplate['cside'] == cell.cside
    ctemplate['cname'].ix[ namemask & sidemask ] += 'a'
    
    # create the new cell with the new name and the same side and position 
    # the new cell is the posterior so the position is to the right
    newpos = np.array( [ cell.cXpos+10, cell.cYpos, cell.cZpos ] )
    newname = cell.cname + 'p'
    newcell = create_cell( newpos, bd, cell.cside, newname)
    
    # add the new cell to ctemplate and return it
    ctemplate = pd.concat([ctemplate,newcell]).sort(['cside','cXpos']).reset_index(drop=True)

    return ctemplate



def remove_cell( cellstp, ctemplate, refpos, ftp ):

    # find the closest cell in the current timepoint
    idx, cell = closer_cell( refpos, cellstp, ftp )

    # remove the correspondent cell in the template
    namemask = ctemplate['cname'] == cell.cname
    sidemask = ctemplate['cside'] == cell.cside
    jdx = ctemplate[ namemask & sidemask ].index[0]    
    ctemplate = ctemplate.drop([jdx])
    ctemplate = ctemplate.sort(['cside','cXpos']).reset_index(drop=True)

    return ctemplate
    


def automatic_cell_names( cellstp, ctemplate ):
    
    cellstp = cellstp.sort(['cside','cXpos'])
    ctemplate = ctemplate.sort(['cside','cXpos'])
    
    newcells = pd.DataFrame( {} )    
    
    side= ['L','R']
    for s in side:
        cside = cellstp[cellstp['cside']==s].reset_index(drop=True)
        ctside = ctemplate[ctemplate['cside']==s].reset_index(drop=True)
        
        if len(cside.cname.values) <= len(ctside.cname.values):
            for j, c in cside.iterrows():
                cside.ix[j,'cname'] = ctside.ix[j,'cname']
                
        else:
            for j, c in ctside.iterrows():
                cside.ix[j,'cname'] = ctside.ix[j,'cname']
        
        newcells = pd.concat([newcells,cside])
        
    newcells = newcells.sort(['cside','cXpos']).reset_index(drop=True)
    newcells.index = cellstp.index

    return newcells


def insert_side( refpos, bd ):
    
    '''
    Automatically update the left-right side of the cells based on the relative position in the stack
    '''

    spline = bd.spline.values[0]

    if np.abs(refpos[-1] - spline[refpos[0],2]) < 3:

        return 'L'

    else:

        return 'R'


def invert_cells_side(cells, key, pos, ftp, idx = None):
    
    '''
    Swap all sides of the cells, or just one cell if 'ctrl+i' is pressed
    '''

    sides = ['L', 'R']


    if key == 'u':
        
        if not idx:
            idx, cell = closer_cell(pos, cells, ftp)
#        print('\n\n old:',cells.ix[idx])
        cells.ix[idx,'cside'] = sides[sides.index(cells.ix[idx,'cside']) - 1]
#        print('\n\n new:',cells.ix[idx])

    else:

        for idx, cell in cells.iterrows():
            cells.ix[idx,'cside'] = sides[sides.index(cells.ix[idx,'cside']) - 1]

    return cells


def cell_trajectory( path, cellid, cells, relative = False, complete = True ):
    
    '''
    COMMENTA STA CAZZO DI FUNZIONE!!!!!!!!
    '''
    
    pxlsize = 6.5/40
    traj = np.nan*np.zeros( ( len(cellid), len(cells) ) )
	
    if not os.path.isfile( path + '\\maxpos.pickle' ):
        maxpos = find_max_pos(path)
    else:
        maxpos = pickle.load( open( path + '\\maxpos.pickle', 'rb' ) )
    
    for idx2, cn in enumerate( cellid ):
#        print(cn)

        for idx1, cs in enumerate( cells ):
            if cs != []:
#                print(cs, [ c[2][0] for c in cs if c[-1] == cn[1] ])

                minpos = np.min( [ c[2][0] for c in cs if c[-1] == cn[1] ] )
                try:
                    if complete:
                        j = [ c[0] in cn[0] and c[-1] == cn[1] for c in cs].index(True)
                    else:
                        j = [ c[0] == cn[0] and c[-1] == cn[1] for c in cs].index(True)
                    if relative:
                        traj[idx2,idx1] = ( cs[j][2][0] - minpos ) / (maxpos[idx1]-minpos)
                    else:
                        traj[idx2,idx1] = ( cs[j][2][0] - minpos ) * pxlsize
                except:
                    continue
            
        if not complete and len(cn[0])>2:
            firsttp = np.where( np.isfinite( traj[idx2,:]) )[0][0]
#            print(cn,firsttp,cn[0][:-1],cn[1],firsttp-1)
            prevtp = firsttp-1
            while cells[prevtp] == []:
                prevtp -= 1
            j = [ c[0] == cn[0][:-1] and c[-1] == cn[1] for c in cells[prevtp] ].index(True)
            minpos = np.min( [ c[2][0] for c in cells[prevtp] if c[-1] == cn[1] ] )
            if relative:
                traj[idx2,prevtp] = ( cells[prevtp][j][2][0] - minpos ) / (maxpos[prevtp]-minpos)
            else:
                traj[idx2,prevtp] = ( cells[prevtp][j][2][0] - minpos ) * pxlsize

    return traj


def find_max_pos( path ):
    flist = glob.glob(path + '\\straight*.tif')
    flist.sort()
    cells = pickle.load( open( path + '\\cells.pickle', 'rb' ) )
	
    pos = np.zeros(len(flist))
    for idx, i in enumerate( zip(flist,cells) ):
#        if i[1] != []:
#            pos[idx] = np.max([j[2][0] for j in i[1]])
#        print(i[0])
        img = loadstack(i[0])
        #p = np.max( [ c[2][0] for c in i[1] ] )
        try:           
            pos[idx] = np.where(img==0)[2][0]#np.min(img[:,:,p:]))[2][0]
        except:
            pos[idx] = img.shape[2]
	
    return np.array(pos).astype('uint16')


def bytes_to_unicode(ob):
    
    '''
    Convert .p files pickled in python 2 in readable python 3 files
    
    EXAMPLE:
    with open(filename,'rb') as fo:
        data = pickle.load(fo,encoding='bytes')
        ndata = bytes_to_unicode(data)
        pickle.dump(ndata, open(filename[:-7] + '3.pickle', "wb"))
    '''
    
    t = type(ob)
    if t in (list, tuple):
        l = [str(i, 'utf-8') if type(i) is bytes else i for i in ob]
        l = [bytes_to_unicode(i) if type(i) in (list, tuple, dict) else i for i in l]
        ro = tuple(l) if t is tuple else l
    elif t is dict:
        byte_keys = [i for i in ob if type(i) is bytes]
        for bk in byte_keys:
            v = ob[bk]
            del(ob[bk])
            ob[str(bk,'utf-8')] = v
        for k in ob:
            if type(ob[k]) is bytes:
                ob[k] = str(ob[k], 'utf-8')
            elif type(ob[k]) in (list, tuple, dict):
                ob[k] = bytes_to_unicode(ob[k])
        ro = ob
    else:
        ro = ob
        print("unprocessed object: {0} {1}".format(t, ob))
    return ro

