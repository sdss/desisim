'''
Code for quickly generating an output zcatalog given fiber assignment tiles,
a truth catalog, and optionally a previous zcatalog.

This uses templates obtained from redmonster.  Given an input z, the code finds the template with the smallest 
redshift that is greater than the input redhsift.  It then takes the difference betweeb the template true redshift and 
the template output redshift, and gives that interval to the input redshift to produce the output redshift. 

TODO:
- Include magnitudes or [OII] flux as part of parameterizing results
'''

from __future__ import absolute_import, division, print_function

import os.path
from collections import Counter
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column

import sys
from operator import itemgetter, attrgetter, methodcaller

import astropy.constants
c = astropy.constants.c.to('km/s').value

from desitarget.targets import desi_mask



def very_new_get_observed_redshifts(truetype, truez):
    """
    Returns observed z, zerr, zwarn arrays given true object types and redshifts
    
    Args:
        truetype : array of ["ELG","LRG","QSO","STD","SKY","BGS","MWS_STAR"]
        truez: array of true redshifts
        
    Returns tuple of (zout, zerr, zwarn) for each category: elg, etc.
    """
    hdulist=fits.open('/project/projectdirs/desi/spectro/quickcat/z_lookuptab_oak1.fits')
    hdud=hdulist[1].data
    target_types=["ELG","LRG","QSO","STD","SKY","BGS","MWS_STAR"]
    tlist=[]
    zout =  np.zeros(len(truez), dtype=np.float32)
    zerr = np.zeros(len(truez), dtype=np.float32)
    zwarn = np.zeros(len(truez), dtype=np.int32)
    for i in range(len(target_types)):
        x=(hdud['TROBJTYPE']==target_types[i])
        tlist.append(hdud[x])
        print(" %d  templates of type %s "%(len(tlist[i]),target_types[i]))
    

    for t_type in range(len(target_types)):
 
        newzb=tlist[t_type]
        #true z values for galaxies run through redmonster for this t_type
        truezarray=newzb['ZBTRUEZ']
        number_templates=len(truezarray)
        #bin size in z
        #set mean templates per bin
        mean=10.
        bin_delta=mean/number_templates
        zmin=np.amin(truezarray)
        zmax=np.amax(truezarray)
        number_bins=int((zmax-zmin)/bin_delta)+1
        print("%s : number of bins %d  bin_delta %f "%(target_types[t_type],number_bins,bin_delta))

        #   put lookup table lines into bins for this t_type
        print(" min %f max %f "%(zmin,zmax))
        A=[]
        x=zmin
        B=[]
        template=0
        for i in range(number_bins):
            #first template not included previously goes into the next collection
            #it is possible that a template might be the only one in a collection 
            #several times
            B.append(template)
           #B is collection of templates belonging in this bin
            while(template<number_templates-1 and truezarray[template+1]<x+bin_delta):
                template=template+1
                B.append(template)
            #add this collection of templates to A, 
            #which is the list of collections of templates    
            A.append(B)
            x=x+bin_delta
            B=[]
        print(" no of bins %d"%len(A))    
        ii, = np.where(truetype == target_types[t_type])
        
        for i in ii:
            zin=truez[i]

            if ((zin>zmax) or (zin<zmin)):
                #print('bad z  %f'%zin)
                zout[i]=0.
                zerr[i]=0.
                zwarn[i]=4
 
            else:    
                #find the right bin for input z
                bin_find=int(  (zin-zmin  )/bin_delta)  
                if(bin_find>number_bins):
                    print('z %f bin_find %d'%(z,bin_find))
                    zout[i]=0.
                    zerr[i]=0.
                    zwarn[i]=4
                else:
                    #use collection of templates in this bin
                    choices=A[bin_find]
                    right_choice=choices[0]
                    if len(choices)!=1:
                        for j in range(len(choices)):
                            #take first template whose true z is less than input z
                            #if all template-z's are greater than input take choices[0]
                            if zin>truezarray[choices[j]]:
                                right_choice=choices[j]
                            else:
                                break


                    zout[i] = truez[i]+(newzb['ZBFITZ'][right_choice]-newzb['ZBTRUEZ'][right_choice])
                    zerr[i] =  newzb['ZBZERR'][right_choice]
                    zwarn[i] = newzb['ZBZWARN'][right_choice]

        
    return zout, zerr, zwarn    


def quickcat(tilefiles, targets, truth, zcat=None, perfect=False,newversion=True):
    """
    Generates quick output zcatalog

    Args:
        tilefiles : list of fiberassign tile files that were observed
        targets : astropy Table of targets
        truth : astropy Table of input truth with columns TARGETID, TRUEZ, and TRUETYPE
        zcat (Optional): input zcatalog Table from previous observations
        perfect (Optional): if True, treat spectro pipeline as perfect with input=output,
            otherwise add noise and zwarn!=0 flags
        
    Returns:
        zcatalog astropy Table based upon input truth, plus ZERR, ZWARN,
        NUMOBS, and TYPE columns

    TODO: Add BGS, MWS support
    """



    #- convert to Table for easier manipulation
    truth = Table(truth)

    #- Count how many times each target was observed for this set of tiles
    ### print('Reading {} tiles'.format(len(obstiles)))
    nobs = Counter()
    print(" number of tilefiles %d" %len(tilefiles))
    for infile in tilefiles:
        #these files are produced by fiberassign
        fibassign = fits.getdata(infile, 'FIBER_ASSIGNMENTS')
        ii = (fibassign['TARGETID'] != -1)  #- targets with assignments
        nobs.update(fibassign['TARGETID'][ii])

    #- Count how many times each target was observed in previous zcatalog
    #- NOTE: assumes that no tiles have been repeated
    if zcat is not None:
        ### print('Counting targets from previous zobs')
        for targetid, n in zip(zcat['TARGETID'], zcat['NUMOBS']):
            nobs[targetid] += n

    #- Trim truth down to just ones that have already been observed
    ### print('Trimming truth to just observed targets')

    obs_targetids = np.array(list(nobs))
    print("len obs_targetids %d"%len(obs_targetids))

    iiobs = np.in1d(truth['TARGETID'], obs_targetids)
    truth = truth[iiobs]
    targets = targets[iiobs]

    #- Construct initial new z catalog
    newzcat = Table()
    newzcat['TARGETID'] = truth['TARGETID']

    if 'BRICKNAME' in truth.dtype.names:
        newzcat['BRICKNAME'] = truth['BRICKNAME']
    else:
        newzcat['BRICKNAME'] = np.zeros(len(truth), dtype=(str, 8))

    #- Copy TRUEZ -> Z so that we can add errors without altering original
    newzcat['Z'] = truth['TRUEZ'].copy()

    newzcat['SPECTYPE'] = truth['TRUETYPE'].copy()
    # rnc add RA and DEC
    newzcat['RA'] = truth['RA'].copy()
    newzcat['DEC'] = truth['DEC'].copy() 
    #- Add numobs column
    ### print('Adding NUMOBS column')
    nz = len(newzcat)
    print(" nz  %d" %nz)
    newzcat.add_column(Column(name='NUMOBS', length=nz, dtype=np.int32))
    for i in range(nz):
        newzcat['NUMOBS'][i] = nobs[newzcat['TARGETID'][i]]

    #- Add ZERR and ZWARN
    ### print('Adding ZERR and ZWARN')

    if not perfect:
        print("zcat using templates from redmonster")
        #- GALAXY -> ELG or LRG
        objtype = newzcat['SPECTYPE'].copy()
        isLRG = (objtype == 'GALAXY') & ((targets['DESI_TARGET'] & desi_mask.LRG) != 0)
        isELG = (objtype == 'GALAXY') & ((targets['DESI_TARGET'] & desi_mask.ELG) != 0)
        objtype[isLRG] = 'LRG'
        objtype[isELG] = 'ELG'

        z, zerr, zwarn = very_new_get_observed_redshifts(objtype, newzcat['Z'])

        newzcat['Z'] = z  #- update with noisy redshift
    else:
        zerr = np.zeros(nz, dtype=np.float32)
        zwarn = np.zeros(nz, dtype=np.int32)

    newzcat['ZERR'] = zerr
    newzcat['ZWARN'] = zwarn
   
    #- Metadata for header
    newzcat.meta['EXTNAME'] = 'ZCATALOG'
    print(" galaxies in zcat %d" %len(newzcat))
    return newzcat
