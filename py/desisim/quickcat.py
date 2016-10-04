'''
Code for quickly generating an output zcatalog given fiber assignment tiles,
a truth catalog, and optionally a previous zcatalog.

The current redshift errors and ZWARN completeness are based upon Redmonster
performance on the zdc1 training samples, documented by Govinda Dhungana at
https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=1657

TODO:
- Include magnitudes or [OII] flux as part of parameterizing results
'''

from __future__ import absolute_import, division, print_function

import os.path
from collections import Counter
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from desispec.log import get_logger
import sys
from operator import itemgetter, attrgetter, methodcaller

import astropy.constants
c = astropy.constants.c.to('km/s').value

from desitarget.targets import desi_mask

#- redshift errors and zwarn fractions from DESI-1657
#- sigmav = c sigmaz / (1+z)
_sigma_v = {
    'ELG': 19.,
    'LRG': 40.,
    'QSO': 423.,
    'STAR': 18.,
    'SKY': 9999,      #- meaningless
    'UNKNOWN': 9999,  #- meaningless
}

_zwarn_fraction = {
    'ELG': 0.14,       # 1 - 4303/5000
    'LRG': 0.015,      # 1 - 4921/5000
    'QSO': 0.18,       # 1 - 4094/5000 
    'STAR': 0.238,     # 1 - 3811/5000
    'SKY': 1.0,
    'UNKNOWN': 1.0,
}

_z_range = {
    'ELG': (0.6,1.7),
    'LRG': (0.4,1.2),
    'QSO': (0.0,3.5),
    'STAR': (-0.005,0.005),
    'SKY': (-0.005,0.005),
    'UNKNOWN': (0.,3.5),
}
def old_get_observed_redshifts(truetype, truez):
    """
        Returns observed z, zerr, zwarn arrays given true object types and redshifts
        
        Args:
        truetype : array of ELG, LRG, QSO, STAR, SKY, or UNKNOWN
        truez: array of true redshifts
        
        Returns tuple of (zout, zerr, zwarn)
        
        TODO: Add BGS, MWS support
        """
    print('old version')
    zout = truez.copy()
    zerr = np.zeros(len(truez), dtype=np.float32)
    zwarn = np.zeros(len(truez), dtype=np.int32)
    for objtype in _sigma_v.keys():
        ii = (truetype == objtype)
        n = np.count_nonzero(ii)
        zerr[ii] = _sigma_v[objtype] * (1+truez[ii]) / c
        zout[ii] += np.random.normal(scale=zerr[ii])
        #- randomly select some objects to set zwarn
        num_zwarn = int(_zwarn_fraction[objtype] * n)
        if num_zwarn > 0:
            jj = np.random.choice(np.where(ii)[0], size=num_zwarn, replace=False)
            zwarn[jj] = 4

    return zout, zerr, zwarn


def new_get_observed_redshifts(truetype, truez):
    """
    Returns observed z, zerr, zwarn arrays given true object types and redshifts
    
    Args:
        truetype : array of ELG, LRG, QSO, STAR, SKY, or UNKNOWN
        truez: array of true redshifts
        
    Returns tuple of (zout, zerr, zwarn)

    TODO: Add BGS, MWS support     
    """
    print('new version')
    log=get_logger()
    zout = truez.copy()
    zerr = np.zeros(len(truez), dtype=np.float32)
    zwarn = np.zeros(len(truez), dtype=np.int32)

    #- reads lookup table from $DESI_ROOT/spectro/quickcat/zbest-zdc1-redmonster-MIX-5000.fits

    desi_root = os.getenv('DESI_ROOT', '/project/projectdirs/desi')
    quickcat_z_lookup = desi_root+'/'+'spectro/quickcat/'
    file = os.path.join(quickcat_z_lookup,'zbest-zdc1-redmonster-MIX-5000.fits')
    try :
        zb_hdulist=fits.open(file)
    except :
        log.error("can not open file %s:"%file)
        zb_hdulist.close()
        sys.exit(12)

    zb_name = zb_hdulist[1].name
    zb=zb_hdulist[zb_name].data

    #- matches redmonster types with DESI target classes
    for objtype in _z_range.keys():
        if ((objtype == 'ELG') | (objtype == 'LRG')):
            rmtype = 'ssp_em_galaxy'
        elif ((objtype == 'STAR') | (objtype == 'SKY') | (objtype == 'QSO_BAD')):
            rmtype = 'spEigenStar'
        else:
            rmtype = objtype

        ii, = np.where(truetype == objtype)
        jj, = np.where((zb['ZBTYPE'] == rmtype))

    #- matches truez with lookup table by target class by looking for the closest redshift from truez
        for i in ii:
            if (len(jj) != 0): #
                newzb = zb[jj]
                diffz = np.abs(newzb['ZBTRUEZ']-truez[i])
                iz, = np.where(diffz == np.min(diffz))
                if (len(iz)>1): iz = np.random.choice(iz,1)
                assert (len(iz)==1),"There should be only one match" #- if more than one match, choose one at random
                zmatch = newzb['ZBTRUEZ'][iz]
                #- fills output catalog
                zout[i] = truez[i]+(newzb['ZBFITZ'][iz]-newzb['ZBTRUEZ'][iz])
                zerr[i] =  newzb['ZBZERR'][iz]
                zwarn[i] = newzb['ZBZWARN'][iz]

        
    return zout, zerr, zwarn    
def very_new_get_observed_redshifts(truetype, truez):
    """
    Returns observed z, zerr, zwarn arrays given true object types and redshifts
    
    Args:
        truetype : array of ELG, LRG, QSO, STAR, SKY, or UNKNOWN
        truez: array of true redshifts
        
    Returns tuple of (zout, zerr, zwarn)

    TODO: Add BGS, MWS support     
    """
    print('very new version')
    log=get_logger()
    zout = truez.copy()
    zerr = np.zeros(len(truez), dtype=np.float32)
    zwarn = np.zeros(len(truez), dtype=np.int32)

    #- reads lookup table from $DESI_ROOT/spectro/quickcat/sorted_true_elgs.fits

    desi_root = os.getenv('DESI_ROOT', '/project/projectdirs/desi')
    quickcat_z_lookup = desi_root+'/'+'spectro/quickcat/'
    file = os.path.join(quickcat_z_lookup,'sorted_elgs.fits')
    try :
        zb_hdulist=fits.open(file)
    except :
        log.error("can not open file %s:"%file)
        zb_hdulist.close()
        sys.exit(12)

    newzb=zb_hdulist[1].data
    zarray=newzb['ZBTRUEZ']
    number_templates=len(zarray)
    ii=(newzb['ZBTYPE']=='ssp_em_galaxy')
    truezarray=newzb[ii]['ZBTRUEZ']

    bin_delta=0.5e-3
    zmin=np.amin(truezarray)
    zmax=np.amax(truezarray)
    number_bins=int((zmax-zmin)/bin_delta)+1

#   put lookup table lines into bins
    print(" min %f max %f "%(zmin,zmax))
    A=[]
    x=zmin
    B=[]
    template=0
    for i in range(number_bins):
        B.append(template)
    
        while(template<number_templates and zarray[template+1]<x+bin_delta):
            template=template+1
            B.append(template)
 
        A.append(B)
        x=x+bin_delta
        B=[]
    print(' len A %d'%len(A))
    print( A[len(A)-1])
    ii, = np.where(truetype == 'ELG')
    for i in ii:
        zin=truez[i]
        if(zin>zmax or zin<zmin):
            #print('bad z  %f'%zin)
            
            zout[i]=0.
            zerr[i]=0.
            zwarn[i]=4
        else:    
            bin_find=int(  (zin-zmin  )/bin_delta)  
            if(bin_find>number_bins):
                print('z %f bin_find %d'%(z,bin_find))
                zout[i]=0.
                zerr[i]=0.
                zwarn[i]=4
            else:
                choices=A[bin_find]
                right_choice=choices[0]
                if len(choices)!=1:
                    for j in range(len(choices)):
                        if zin>zarray[choices[j]]:
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
        
    Options:
        zcat : input zcatalog Table from previous observations
        perfect : if True, treat spectro pipeline as perfect with input=output,
            otherwise add noise and zwarn!=0 flags
        
    Returns:
        zcatalog astropy Table based upon input truth, plus ZERR, ZWARN,
        NUMOBS, and TYPE columns   
        
    TODO: Add BGS, MWS support     
    """



    #- convert to Table for easier manipulation
    truth = Table(truth)
    print("length of truth %d"%len(truth))
    #- Count how many times each target was observed for this set of tiles
    ### print('Reading {} tiles'.format(len(obstiles)))
    nobs = Counter()
    print(" number of tilefiles %d" %len(tilefiles))
    for infile in tilefiles:
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
    print("length of truth after iiobs %d" %len(truth))
    targets = targets[iiobs]
    print("length of targets after iiobs %d" %len(targets))

    #- Construct initial new z catalog
    newzcat = Table()
    newzcat['TARGETID'] = truth['TARGETID']
    print ("length of newzcat after TARGETID %d"%len(newzcat))
    if 'BRICKNAME' in truth.dtype.names:
        newzcat['BRICKNAME'] = truth['BRICKNAME']
    else:
        newzcat['BRICKNAME'] = np.zeros(len(truth), dtype='S8')

    #- Copy TRUEZ -> Z so that we can add errors without altering original
    newzcat['Z'] = truth['TRUEZ'].copy()
    newzcat['TYPE'] = truth['TRUETYPE'].copy()
    # rnc add RA and DEC
    newzcat['RA'] = truth['RA'].copy()
    newzcat['DEC'] = truth['DEC'].copy() 
    # rnc add CATEGORY
    #newzcat['CATEGORY']=truth['CATEGORY'].copy()
    #- Add numobs column
    ### print('Adding NUMOBS column')
    nz = len(newzcat)
    print(" nz  %d" %nz)
    newzcat.add_column(Column(name='NUMOBS', length=nz, dtype=np.int32))
    for i in range(nz):
        newzcat['NUMOBS'][i] = nobs[newzcat['TARGETID'][i]]

    #- Add ZERR and ZWARN
    ### print('Adding ZERR and ZWARN')
    print ("perfect z? %s "%perfect)
    if not perfect:
        #- GALAXY -> ELG or LRG
        objtype = newzcat['TYPE'].copy()
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

