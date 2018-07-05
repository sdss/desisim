"""
lvmsim.scripts.quickgen
========================

Quickgen quickly simulates pipeline outputs if given input files.

  - must provide simspec and fibermap files via newexp script
  - Number of spectra to be simulated can be given as an argument for quickgen,
    but the number of spectra in the simspec file is taken by default
  - For this option, airmass and exposure time are keywords given to newexp
  - The keywords provided in the examples are all required, additional keywords
    are provided below
  - Collect a set of templates to simulate as a new exposure::

        newexp --nspec 500 --night 20150915 --expid 0 --flavor dark

  - newexp keyword arguments::

        --flavor : arc/flat/dark/gray/bright/bgs/mws/elg/lrg/qso, type=str, default='dark'
        --tileid : tile id, type=int
        --expid : exposure id, type=int
        --exptime : exposure time in seconds, default for arc = 5s, flat = 10s, dark/elg/lrg/qso=1000s, bright/bgs/mws=300s, type=int
        --night : YEARMMDD, type=str
        --nspec : Number of spectra to simulate, type=int, default=5000
        --airmass : type=float, default=1.0
        --seed : random number seed, type=int
        --testslit : test slit simulation with fewer fibers, action="store_true"
        --arc-lines : alternate arc lines filename, type=str, default=None
        --flat-spectrum : alternate flat spectrum filename, type=str

  - Actually do the simulation::

        simdir=$LVM_SPECTRO_SIM/$PIXPROD/20150915
        quickgen --simspec $simdir/simspec-00000000.fits --fibermap $simdir/fibermap-00000000.fits

  - quickgen output (can also provide frame file only as keyword)::

        1. frame-{camera}-{expid}.fits : raw extracted photons with no calibration at all
        2. sky-{camera}-{expid}.fits : the sky model in photons
        3. fluxcalib-{camera}-{expid}.fits] : the flux calibration vectors
        4. cframe-{camera}-{expid}.fits : flux calibrated spectra

        These files are written to $simdir/{expid}

  - nspec, config, seed, moon-phase, moon-angle, moon-zenith
  - simspec, fibermap, nstart, spectrograph, frameonly
    zrange-qso, zrange-elg, zrange-lrg, zrange-bgs, sne-rfluxratiorange,
    add-SNeIa
"""
from __future__ import absolute_import, division, print_function

import argparse
import astropy.units as u
from astropy.io import fits
from astropy.table import Table, Column, vstack
from time import asctime

import os
import os.path
import numpy as np
import scipy.special
import scipy.interpolate
import scipy.constants as const
import sys

import lvmspec
import lvmspec.io
import lvmsim
import lvmsim.io
import lvmsim.templates
import lvmutil.io
from lvmspec.resolution import Resolution
from lvmspec.io import write_flux_calibration, write_fiberflat, read_fibermap, specprod_root, fitsheader, empty_fibermap
from lvmspec.interpolation import resample_flux
from lvmspec.frame import Frame
from lvmspec.fiberflat import FiberFlat
from lvmspec.sky import SkyModel
from lvmspec.fluxcalibration import FluxCalib
from lvmutil.log import get_logger, DEBUG, INFO
from lvmsim.obs import get_night
from lvmsim.targets import sample_objtype
from lvmsim.specsim import get_simulator
from lvmmodel.io import load_lvmparams
from lvmsim.simexp import get_source_types


def _add_truth(hdus, header, meta, trueflux, sflux, wave, channel):
    """Utility function for adding truth to an output FITS file."""
    hdus.append(
        fits.ImageHDU(trueflux[channel], name='_TRUEFLUX', header=header))
    if channel == 'b':
        swave = wave.astype(np.float32)
        hdus.append(fits.ImageHDU(swave, name='_SOURCEWAVE', header=header))
        hdus.append(fits.ImageHDU(sflux, name='_SOURCEFLUX', header=header))
        metatable = lvmutil.io.encode_table(meta, encoding='ascii')
        metahdu = fits.convenience.table_to_hdu(meta)
        metahdu.header['EXTNAME'] = '_TRUTH'
        hdus.append(metahdu)


def parse(options=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--simspec', type=str, required=True, help="input simspec file")
    parser.add_argument('--fibermap', type=str, help='input fibermap file')
    parser.add_argument('-n', '--nspec', type=int, default=100, help='number of spectra to be simulated, starting from first')
    parser.add_argument('--nstart', type=int, default=0, help='starting spectra # for simulation 0-4999')
    parser.add_argument('--spectrograph', type=int, default=None, help='Spectrograph no. 0-9')
    parser.add_argument('--config', type=str, default='lvm', help='specsim configuration')
    parser.add_argument('-b', '--n_fibers', type=int, default=650, help='total number of fibers')
    parser.add_argument('-t', '--telescope', type=str, default='1m', help='telescope', choices=['1m', '160mm'])
    parser.add_argument('-l', '--location', type=str, default='APO', help='site location', choices=['APO', 'LCO'])
    parser.add_argument('-s', '--seed', type=int, default=0, help="random seed")
    # Only produce uncalibrated output
    parser.add_argument('--frameonly', action="store_true", help="only output frame files")

    # Moon options if bright or gray time
    parser.add_argument('--moon-phase', type=float, help='moon phase (0=full, 1=new)', default=None, metavar='')
    parser.add_argument('--moon-angle', type=float, help='separation angle to the moon (0-180 deg)', default=None, metavar='')
    parser.add_argument('--moon-zenith', type=float, help='zenith angle of the moon (0-90 deg)', default=None, metavar='')

    parser.add_argument('--objtype', type=str, help='ELG, LRG, QSO, BGS, MWS, WD, DARK_MIX, or BRIGHT_MIX', default='DARK_MIX', metavar='')
    parser.add_argument('-a', '--airmass', type=float, help='airmass', default=None, metavar='')
    parser.add_argument('-e', '--exptime', type=float, help='exposure time (s)', default=None, metavar='')
    parser.add_argument('-o', '--outdir', type=str, help='output directory', default='.', metavar='')
    parser.add_argument('-v', '--verbose', action='store_true', help='toggle on verbose output')
    parser.add_argument('--outdir-truth', type=str, help='optional alternative output directory for truth files', metavar='')

    # Object type specific options
    parser.add_argument('--zrange-qso', type=float, default=(0.5, 4.0), nargs=2, metavar='',
                        help='minimum and maximum redshift range for QSO')
    parser.add_argument('--zrange-elg', type=float, default=(0.6, 1.6), nargs=2, metavar='',
                        help='minimum and maximum redshift range for ELG')
    parser.add_argument('--zrange-lrg', type=float, default=(0.5, 1.1), nargs=2, metavar='',
                        help='minimum and maximum redshift range for LRG')
    parser.add_argument('--zrange-bgs', type=float, default=(0.01, 0.4), nargs=2, metavar='',
                        help='minimum and maximum redshift range for BGS')
    parser.add_argument('--rmagrange-bgs', type=float, default=(15.0, 19.5), nargs=2, metavar='',
                        help='Minimum and maximum BGS r-band (AB) magnitude range')
    parser.add_argument('--sne-rfluxratiorange', type=float, default=(0.1, 1.0), nargs=2, metavar='',
                        help='r-band flux ratio of the SNeIa spectrum relative to the galaxy')
    parser.add_argument('--add-SNeIa', action='store_true', help='include SNeIa spectra')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)

    log = get_logger()
    if args.simspec:
        args.objtype = None
        if args.fibermap is None:
            dirname = os.path.dirname(os.path.abspath(args.simspec))
            filename = os.path.basename(args.simspec).replace('simspec', 'fibermap')
            args.fibermap = os.path.join(dirname, filename)
            log.warning('deriving fibermap {} from simspec input filename'.format(args.fibermap))

    return args


def setup_envs():
    ''' Set up the environment and directory structure '''

    LVM_SPECTRO_REDUX_DIR = "./quickGen"

    if 'LVM_SPECTRO_REDUX' not in os.environ:
        log.info('LVM_SPECTRO_REDUX environment is not set.')
    else:
        LVM_SPECTRO_REDUX_DIR = os.environ['LVM_SPECTRO_REDUX']

    if os.path.exists(LVM_SPECTRO_REDUX_DIR):
        if not os.path.isdir(LVM_SPECTRO_REDUX_DIR):
            raise RuntimeError("Path %s Not a directory" % LVM_SPECTRO_REDUX_DIR)
    else:
        try:
            os.makedirs(LVM_SPECTRO_REDUX_DIR)
        except:
            raise

    SPECPROD_DIR = 'specprod'
    if 'SPECPROD' not in os.environ:
        log.info('SPECPROD environment is not set.')
    else:
        SPECPROD_DIR = os.environ['SPECPROD']
    prod_Dir = specprod_root()

    if os.path.exists(prod_Dir):
        if not os.path.isdir(prod_Dir):
            raise RuntimeError("Path %s Not a directory" % prod_Dir)
    else:
        try:
            os.makedirs(prod_Dir)
        except:
            raise


def get_fibermap(fmapfile, log=None, nspec=None):
    ''' Get the fibermap table

    The fibermap table contains (simulated) information about the
    position of each target in the focal plane

    The fibermap contains the fiber positioner configuration information for
    each exposure: what fiber is placed where, what target that is, etc.

    http://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_DATA/NIGHT/fibermap-EXPID.html

    Parameters:
        fmapfile (str):
            An input fibermap filename
        log (logger):
            A Python logger
        nspec (int):
            The number of spectra to simulate

    Returns:
        fibermap (Table):
            An Astropy table of the simulated fibers on the focal plane.
        objtype (array):
            An array of object source types, one for each fiber
        night (int):
            The night of observation
        expid (int):
            The exposure id

    '''

    objtype = None
    if fmapfile:
        # read the fibermap file to get the object information
        if log:
            log.info("Reading fibermap file {}".format(fmapfile))
        fibermap = read_fibermap(fmapfile)
        objtype = get_source_types(fibermap)
        stdindx = np.where(objtype == 'STD')  # match STD with STAR
        mwsindx = np.where(objtype == 'MWS_STAR')  # match MWS_STAR with STAR
        bgsindx = np.where(objtype == 'BGS')  # match BGS with LRG
        objtype[stdindx] = 'STAR'
        objtype[mwsindx] = 'STAR'
        objtype[bgsindx] = 'LRG'
        night = fibermap.meta['NIGHT']
        expid = fibermap.meta['EXPID']
    else:
        # create a new empty fake fibermap object
        if log:
            log.info('Creating empty fibermap')
        assert nspec is not None, 'nspec is required to create an empty fibermap'
        fibermap = empty_fibermap(nspec)
        targetids = random_state.randint(2**62, size=nspec)
        fibermap['TARGETID'] = targetids
        night = get_night()
        expid = 0

    return fibermap, objtype, night, expid


def get_simspec(simspecfile, log=None, nspec=None):
    ''' Get the simspec object

    The simspec table holds the "truth" spectra and the intrinsic properties
    of each object (redshift, noiseless photometry, [OII] flux, etc.).
    (Input spectra to simulate with pixsim.)

    http://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_SIM/PIXPROD/NIGHT/simspec-EXPID.html

    Parameters:
        simspecfile (str):
            The filename of the input simspec file

    '''

    minwave = 3533.
    maxwave = 9913.1
    stepwave = 0.2
    scale = 1.e17
    exptime = 5.  # typical BOSS exposure time in s

    if simspecfile:
        if log:
            log.info('Reading input file {}'.format(args.simspec))
        # create SimSpec object
        simspec = lvmsim.io.read_simspec(args.simspec)
        # number of spectra to simulate from SimSpec
        sim_nspec = simspec.nspec

        # get the spectra and wavelengths arrays for different flavors
        if simspec.flavor == 'arc':
            # - TODO: do we need quickgen to support arcs?  For full pipeline
            # - arcs are used to measure PSF but aren't extracted except for
            # - debugging.
            # - TODO: if we do need arcs, this needs to be redone.
            # - conversion from phot to flux doesn't include throughput,
            # - and arc lines are rebinned to nearest 0.2 A.

            # Create full wavelength and flux arrays for arc exposure
            wave_b = np.array(simspec.wave['b'])
            wave_r = np.array(simspec.wave['r'])
            wave_z = np.array(simspec.wave['z'])
            phot_b = np.array(simspec.phot['b'][0])
            phot_r = np.array(simspec.phot['r'][0])
            phot_z = np.array(simspec.phot['z'][0])
            sim_wave = np.concatenate((wave_b, wave_r, wave_z))
            sim_phot = np.concatenate((phot_b, phot_r, phot_z))
            wavelengths = np.arange(minwave, maxwave, stepwave)
            phot = np.zeros(len(wavelengths))
            for i in range(len(sim_wave)):
                wavelength = sim_wave[i]
                flux_index = np.argmin(abs(wavelength - wavelengths))
                phot[flux_index] = sim_phot[i]
            # Convert photons to flux: following specter conversion method
            dw = np.gradient(wavelengths)
            fibarea = const.pi * (1.07e-2 / 2) ** 2  # cross-sectional fiber area in cm^2
            hc = scale * const.h * const.c  # convert to erg A
            spectra = (hc * exptime * fibarea * dw * phot) / wavelengths
        else:
            wavelengths = simspec.wave['brz']
            spectra = simspec.flux

        # check there's enough spectra to simulate from what we ask for
        if sim_nspec < nspec:
            log.info("Only {} spectra in input file".format(sim_nspec))
            nspec = sim_nspec
    else:
        # Initialize the output truth table.
        spectra = []
        wavelengths = qsim.source.wavelength_out.to(u.Angstrom).value
        npix = len(wavelengths)
        truth = dict()
        meta = Table()
        truth['OBJTYPE'] = np.zeros(args.nspec, dtype=(str, 10))
        truth['FLUX'] = np.zeros((args.nspec, npix))
        truth['WAVE'] = wavelengths
        jj = list()

        for thisobj in set(true_objtype):
            ii = np.where(true_objtype == thisobj)[0]
            nobj = len(ii)
            truth['OBJTYPE'][ii] = thisobj
            if log:
                log.info('Generating {} template'.format(thisobj))

            # Generate the templates
            if thisobj == 'ELG':
                elg = lvmsim.templates.ELG(wave=wavelengths, add_SNeIa=args.add_SNeIa)
                flux, tmpwave, meta1 = elg.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_elg,
                                                          sne_rfluxratiorange=args.sne_rfluxratiorange)
            elif thisobj == 'LRG':
                lrg = lvmsim.templates.LRG(wave=wavelengths, add_SNeIa=args.add_SNeIa)
                flux, tmpwave, meta1 = lrg.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_lrg,
                                                          sne_rfluxratiorange=args.sne_rfluxratiorange)
            elif thisobj == 'QSO':
                qso = lvmsim.templates.QSO(wave=wavelengths)
                flux, tmpwave, meta1 = qso.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_qso)
            elif thisobj == 'BGS':
                bgs = lvmsim.templates.BGS(wave=wavelengths, add_SNeIa=args.add_SNeIa)
                flux, tmpwave, meta1 = bgs.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_bgs,
                                                          rmagrange=args.rmagrange_bgs,
                                                          sne_rfluxratiorange=args.sne_rfluxratiorange)
            elif thisobj == 'STD':
                fstd = lvmsim.templates.FSTD(wave=wavelengths)
                flux, tmpwave, meta1 = fstd.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'QSO_BAD':  # use STAR template no color cuts
                star = lvmsim.templates.STAR(wave=wavelengths)
                flux, tmpwave, meta1 = star.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'MWS_STAR' or thisobj == 'MWS':
                mwsstar = lvmsim.templates.MWS_STAR(wave=wavelengths)
                flux, tmpwave, meta1 = mwsstar.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'WD':
                wd = lvmsim.templates.WD(wave=wavelengths)
                flux, tmpwave, meta1 = wd.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'SKY':
                flux = np.zeros((nobj, npix))
                meta1 = Table(dict(REDSHIFT=np.zeros(nobj, dtype=np.float32)))
            elif thisobj == 'TEST':
                flux = np.zeros((args.nspec, npix))
                indx = np.where(wave > 5800.0 - 1E-6)[0][0]
                ref_integrated_flux = 1E-10
                ref_cst_flux_density = 1E-17
                single_line = (np.arange(args.nspec) % 2 == 0).astype(np.float32)
                continuum = (np.arange(args.nspec) % 2 == 1).astype(np.float32)

                for spec in range(nspec):
                    flux[spec, indx] = single_line[spec] * ref_integrated_flux / np.gradient(wavelengths)[indx]  # single line
                    flux[spec] += continuum[spec] * ref_cst_flux_density  # flat continuum

                meta1 = Table(dict(REDSHIFT=np.zeros(args.nspec, dtype=np.float32),
                                   LINE=wave[indx] * np.ones(args.nspec, dtype=np.float32),
                                   LINEFLUX=single_line * ref_integrated_flux,
                                   CONSTFLUXDENSITY=continuum * ref_cst_flux_density))
            else:
                raise RuntimeError('Unknown object type')

            # Pack it in.
            truth['FLUX'][ii] = flux
            meta = vstack([meta, meta1])
            jj.append(ii.tolist())

            # Sanity check on units; templates currently return ergs, not 1e-17 ergs...
            # assert (thisobj == 'SKY') or (np.max(truth['FLUX']) < 1e-6)

        # Sort the metadata table.
        jj = sum(jj, [])
        meta_new = Table()
        for k in range(nspec):
            index = int(np.where(np.array(jj) == k)[0])
            meta_new = vstack([meta_new, meta[index]])
        meta = meta_new

        # Add TARGETID and the true OBJTYPE to the metadata table.
        meta.add_column(Column(true_objtype, dtype=(str, 10), name='TRUE_OBJTYPE'))
        meta.add_column(Column(targetids, name='TARGETID'))

        # Rename REDSHIFT -> TRUEZ anticipating later table joins with zbest.Z
        meta.rename_column('REDSHIFT', 'TRUEZ')

    return spectra, wavelengths, nspec


def main(args):

    # Set up the logger
    if args.verbose:
        log = get_logger(DEBUG)
    else:
        log = get_logger()

    # Make sure all necessary environment variables are set
    setup_envs()

    # Initialize random number generator to use.
    np.random.seed(args.seed)
    random_state = np.random.RandomState(args.seed)

    # Derive spectrograph number from nstart if needed
    if args.spectrograph is None:
        args.spectrograph = args.nstart / args.n_fibers

    # Read fibermapfile to get object type, night and expid
    fibermap, objtype, night, expid = get_fibermap(args.fibermap, log=log, nspec=args.nspec)

    # Initialize the spectral simulator
    log.info("Initializing SpecSim with config {}".format(args.config))
    lvmparams = load_lvmparams(config=args.config, telescope=args.telescope)
    qsim = get_simulator(args.config, num_fibers=1, params=lvmparams)

    if args.simspec:
        # Read the input file
        log.info('Reading input file {}'.format(args.simspec))
        simspec = lvmsim.io.read_simspec(args.simspec)
        nspec = simspec.nspec
        if simspec.flavor == 'arc':
            # - TODO: do we need quickgen to support arcs?  For full pipeline
            # - arcs are used to measure PSF but aren't extracted except for
            # - debugging.
            # - TODO: if we do need arcs, this needs to be redone.
            # - conversion from phot to flux doesn't include throughput,
            # - and arc lines are rebinned to nearest 0.2 A.

            # Create full wavelength and flux arrays for arc exposure
            wave_b = np.array(simspec.wave['b'])
            wave_r = np.array(simspec.wave['r'])
            wave_z = np.array(simspec.wave['z'])
            phot_b = np.array(simspec.phot['b'][0])
            phot_r = np.array(simspec.phot['r'][0])
            phot_z = np.array(simspec.phot['z'][0])
            sim_wave = np.concatenate((wave_b, wave_r, wave_z))
            sim_phot = np.concatenate((phot_b, phot_r, phot_z))
            wavelengths = np.arange(3533., 9913.1, 0.2)
            phot = np.zeros(len(wavelengths))
            for i in range(len(sim_wave)):
                wavelength = sim_wave[i]
                flux_index = np.argmin(abs(wavelength - wavelengths))
                phot[flux_index] = sim_phot[i]
            # Convert photons to flux: following specter conversion method
            dw = np.gradient(wavelengths)
            exptime = 5.  # typical BOSS exposure time in s
            fibarea = const.pi * (1.07e-2 / 2) ** 2  # cross-sectional fiber area in cm^2
            hc = 1.e17 * const.h * const.c  # convert to erg A
            spectra = (hc * exptime * fibarea * dw * phot) / wavelengths
        else:
            wavelengths = simspec.wave['brz']
            spectra = simspec.flux
        if nspec < args.nspec:
            log.info("Only {} spectra in input file".format(nspec))
            args.nspec = nspec

    else:
        # Initialize the output truth table.
        spectra = []
        wavelengths = qsim.source.wavelength_out.to(u.Angstrom).value
        npix = len(wavelengths)
        truth = dict()
        meta = Table()
        truth['OBJTYPE'] = np.zeros(args.nspec, dtype=(str, 10))
        truth['FLUX'] = np.zeros((args.nspec, npix))
        truth['WAVE'] = wavelengths
        jj = list()

        for thisobj in set(true_objtype):
            ii = np.where(true_objtype == thisobj)[0]
            nobj = len(ii)
            truth['OBJTYPE'][ii] = thisobj
            log.info('Generating {} template'.format(thisobj))

            # Generate the templates
            if thisobj == 'ELG':
                elg = lvmsim.templates.ELG(wave=wavelengths, add_SNeIa=args.add_SNeIa)
                flux, tmpwave, meta1 = elg.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_elg,
                                                          sne_rfluxratiorange=args.sne_rfluxratiorange)
            elif thisobj == 'LRG':
                lrg = lvmsim.templates.LRG(wave=wavelengths, add_SNeIa=args.add_SNeIa)
                flux, tmpwave, meta1 = lrg.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_lrg,
                                                          sne_rfluxratiorange=args.sne_rfluxratiorange)
            elif thisobj == 'QSO':
                qso = lvmsim.templates.QSO(wave=wavelengths)
                flux, tmpwave, meta1 = qso.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_qso)
            elif thisobj == 'BGS':
                bgs = lvmsim.templates.BGS(wave=wavelengths, add_SNeIa=args.add_SNeIa)
                flux, tmpwave, meta1 = bgs.make_templates(nmodel=nobj, seed=args.seed, zrange=args.zrange_bgs,
                                                          rmagrange=args.rmagrange_bgs,
                                                          sne_rfluxratiorange=args.sne_rfluxratiorange)
            elif thisobj == 'STD':
                fstd = lvmsim.templates.FSTD(wave=wavelengths)
                flux, tmpwave, meta1 = fstd.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'QSO_BAD':  # use STAR template no color cuts
                star = lvmsim.templates.STAR(wave=wavelengths)
                flux, tmpwave, meta1 = star.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'MWS_STAR' or thisobj == 'MWS':
                mwsstar = lvmsim.templates.MWS_STAR(wave=wavelengths)
                flux, tmpwave, meta1 = mwsstar.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'WD':
                wd = lvmsim.templates.WD(wave=wavelengths)
                flux, tmpwave, meta1 = wd.make_templates(nmodel=nobj, seed=args.seed)
            elif thisobj == 'SKY':
                flux = np.zeros((nobj, npix))
                meta1 = Table(dict(REDSHIFT=np.zeros(nobj, dtype=np.float32)))
            elif thisobj == 'TEST':
                flux = np.zeros((args.nspec, npix))
                indx = np.where(wave > 5800.0 - 1E-6)[0][0]
                ref_integrated_flux = 1E-10
                ref_cst_flux_density = 1E-17
                single_line = (np.arange(args.nspec) % 2 == 0).astype(np.float32)
                continuum = (np.arange(args.nspec) % 2 == 1).astype(np.float32)

                for spec in range(args.nspec):
                    flux[spec, indx] = single_line[spec] * ref_integrated_flux / np.gradient(wavelengths)[indx]  # single line
                    flux[spec] += continuum[spec] * ref_cst_flux_density  # flat continuum

                meta1 = Table(dict(REDSHIFT=np.zeros(args.nspec, dtype=np.float32),
                                   LINE=wave[indx] * np.ones(args.nspec, dtype=np.float32),
                                   LINEFLUX=single_line * ref_integrated_flux,
                                   CONSTFLUXDENSITY=continuum * ref_cst_flux_density))
            else:
                log.fatal('Unknown object type {}'.format(thisobj))
                sys.exit(1)

            # Pack it in.
            truth['FLUX'][ii] = flux
            meta = vstack([meta, meta1])
            jj.append(ii.tolist())

            # Sanity check on units; templates currently return ergs, not 1e-17 ergs...
            # assert (thisobj == 'SKY') or (np.max(truth['FLUX']) < 1e-6)

        # Sort the metadata table.
        jj = sum(jj, [])
        meta_new = Table()
        for k in range(args.nspec):
            index = int(np.where(np.array(jj) == k)[0])
            meta_new = vstack([meta_new, meta[index]])
        meta = meta_new

        # Add TARGETID and the true OBJTYPE to the metadata table.
        meta.add_column(Column(true_objtype, dtype=(str, 10), name='TRUE_OBJTYPE'))
        meta.add_column(Column(targetids, name='TARGETID'))

        # Rename REDSHIFT -> TRUEZ anticipating later table joins with zbest.Z
        meta.rename_column('REDSHIFT', 'TRUEZ')

    # ---------- end simspec

    # explicitly set location on focal plane if needed to support airmass
    # variations when using specsim v0.5
    if qsim.source.focal_xy is None:
        qsim.source.focal_xy = (u.Quantity(0, 'mm'), u.Quantity(100, 'mm'))

    # Set simulation parameters from the simspec header or lvmparams
    bright_objects = ['bgs', 'mws', 'bright', 'BGS', 'MWS', 'BRIGHT_MIX']
    gray_objects = ['gray', 'grey']
    if args.simspec is None:
        object_type = objtype
        flavor = None
    elif simspec.flavor == 'science':
        object_type = None
        flavor = simspec.header['PROGRAM']
    else:
        object_type = None
        flavor = simspec.flavor
        log.warning('Maybe using an outdated simspec file with flavor={}'.format(flavor))

    # Set airmass
    if args.airmass is not None:
        qsim.atmosphere.airmass = args.airmass
    elif args.simspec and 'AIRMASS' in simspec.header:
        qsim.atmosphere.airmass = simspec.header['AIRMASS']
    else:
        qsim.atmosphere.airmass = 1.25   # Science Req. Doc L3.3.2

    # Set site location
    if args.location is not None:
        qsim.observation.observatory = args.location
    else:
        qsim.observation.observatory = 'APO'

    # Set exptime
    if args.exptime is not None:
        qsim.observation.exposure_time = args.exptime * u.s
    elif args.simspec and 'EXPTIME' in simspec.header:
        qsim.observation.exposure_time = simspec.header['EXPTIME'] * u.s
    elif objtype in bright_objects:
        qsim.observation.exposure_time = lvmparams['exptime_bright'] * u.s
    else:
        qsim.observation.exposure_time = lvmparams['exptime_dark'] * u.s

    # Set Moon Phase
    if args.moon_phase is not None:
        qsim.atmosphere.moon.moon_phase = args.moon_phase
    elif args.simspec and 'MOONFRAC' in simspec.header:
        qsim.atmosphere.moon.moon_phase = simspec.header['MOONFRAC']
    elif flavor in bright_objects or object_type in bright_objects:
        qsim.atmosphere.moon.moon_phase = 0.7
    elif flavor in gray_objects:
        qsim.atmosphere.moon.moon_phase = 0.1
    else:
        qsim.atmosphere.moon.moon_phase = 0.5

    # Set Moon Zenith
    if args.moon_zenith is not None:
        qsim.atmosphere.moon.moon_zenith = args.moon_zenith * u.deg
    elif args.simspec and 'MOONALT' in simspec.header:
        qsim.atmosphere.moon.moon_zenith = simspec.header['MOONALT'] * u.deg
    elif flavor in bright_objects or object_type in bright_objects:
        qsim.atmosphere.moon.moon_zenith = 30 * u.deg
    elif flavor in gray_objects:
        qsim.atmosphere.moon.moon_zenith = 80 * u.deg
    else:
        qsim.atmosphere.moon.moon_zenith = 100 * u.deg

    # Set Moon - Object Angle
    if args.moon_angle is not None:
        qsim.atmosphere.moon.separation_angle = args.moon_angle * u.deg
    elif args.simspec and 'MOONSEP' in simspec.header:
        qsim.atmosphere.moon.separation_angle = simspec.header['MOONSEP'] * u.deg
    elif flavor in bright_objects or object_type in bright_objects:
        qsim.atmosphere.moon.separation_angle = 50 * u.deg
    elif flavor in gray_objects:
        qsim.atmosphere.moon.separation_angle = 60 * u.deg
    else:
        qsim.atmosphere.moon.separation_angle = 60 * u.deg

    # Initialize per-camera output arrays that will be saved
    waves, trueflux, noisyflux, obsivar, resolution, sflux = {}, {}, {}, {}, {}, {}

    maxbin = 0
    nmax = args.nspec
    for camera in qsim.instrument.cameras:
        # Lookup this camera's resolution matrix and convert to the sparse format used in lvmspec.
        R = Resolution(camera.get_output_resolution_matrix())
        resolution[camera.name] = np.tile(R.to_fits_array(), [args.nspec, 1, 1])
        waves[camera.name] = (camera.output_wavelength.to(u.Angstrom).value.astype(np.float32))
        nwave = len(waves[camera.name])
        maxbin = max(maxbin, len(waves[camera.name]))
        nobj = np.zeros((nmax, 3, maxbin))  # object photons
        nsky = np.zeros((nmax, 3, maxbin))  # sky photons
        nivar = np.zeros((nmax, 3, maxbin))  # inverse variance (object+sky)
        cframe_observedflux = np.zeros((nmax, 3, maxbin))   # calibrated object flux
        cframe_ivar = np.zeros((nmax, 3, maxbin))  # inverse variance of calibrated object flux
        cframe_rand_noise = np.zeros((nmax, 3, maxbin))  # random Gaussian noise to calibrated flux
        sky_ivar = np.zeros((nmax, 3, maxbin))  # inverse variance of sky
        sky_rand_noise = np.zeros((nmax, 3, maxbin))  # random Gaussian noise to sky only
        frame_rand_noise = np.zeros((nmax, 3, maxbin))  # random Gaussian noise to nobj+nsky
        trueflux[camera.name] = np.empty((args.nspec, nwave))  # calibrated flux
        noisyflux[camera.name] = np.empty((args.nspec, nwave))  # observed flux with noise
        obsivar[camera.name] = np.empty((args.nspec, nwave))  # inverse variance of flux
        if args.simspec:
            dw = np.gradient(simspec.wave[camera.name])
        else:
            sflux = np.empty((args.nspec, npix))

    # - Check if input simspec is for a continuum flat lamp instead of science
    # - This does not convolve to per-fiber resolution
    if args.simspec:
        if simspec.flavor == 'flat':
            log.info("Simulating flat lamp exposure")
            for i, camera in enumerate(qsim.instrument.cameras):
                channel = camera.name
                assert camera.output_wavelength.unit == u.Angstrom
                num_pixels = len(waves[channel])
                dw = np.gradient(simspec.wave[channel])
                meanspec = resample_flux(
                    waves[channel], simspec.wave[channel],
                    np.average(simspec.phot[channel] / dw, axis=0))
                fiberflat = random_state.normal(loc=1.0,
                                                scale=1.0 / np.sqrt(meanspec),
                                                size=(nspec, num_pixels))
                ivar = np.tile(meanspec, [nspec, 1])
                mask = np.zeros((simspec.nspec, num_pixels), dtype=np.uint32)

                for kk in range((args.nspec + args.nstart - 1) // args.n_fibers + 1):
                    camera = channel + str(kk)
                    outfile = lvmspec.io.findfile('fiberflat', night, expid, camera)
                    start = max(args.n_fibers * kk, args.nstart)
                    end = min(args.n_fibers * (kk + 1), nmax)

                    if (args.spectrograph <= kk):
                        log.info("Writing files for channel:{}, spectrograph:{}, spectra:{} to {}".format(channel, kk, start, end))

                    ff = FiberFlat(
                        waves[channel], fiberflat[start:end, :],
                        ivar[start:end, :], mask[start:end, :], meanspec,
                        header=dict(CAMERA=camera))
                    write_fiberflat(outfile, ff)
                    filePath = lvmspec.io.findfile("fiberflat", night, expid, camera)
                    log.info("Wrote file {}".format(filePath))

            sys.exit(0)

    # Repeat the simulation for all spectra
    scale = 1e-17
    fluxunits = scale * u.erg / (u.s * u.cm ** 2 * u.Angstrom)
    for j in range(args.nspec):

        thisobjtype = objtype[j]
        sys.stdout.flush()
        if flavor == 'arc':
            qsim.source.update_in(
                'Quickgen source {0}'.format, 'perfect',
                wavelengths * u.Angstrom, spectra * fluxunits)
        else:
            qsim.source.update_in(
                'Quickgen source {0}'.format(j), thisobjtype.lower(),
                wavelengths * u.Angstrom, spectra[j, :] * fluxunits)
        qsim.source.update_out()

        qsim.simulate()
        qsim.generate_random_noise(random_state)

        for i, output in enumerate(qsim.camera_output):
            assert output['observed_flux'].unit == 1e17 * fluxunits
            # Extract the simulation results needed to create our uncalibrated
            # frame output file.
            num_pixels = len(output)
            nobj[j, i, :num_pixels] = output['num_source_electrons'][:, 0]
            nsky[j, i, :num_pixels] = output['num_sky_electrons'][:, 0]
            nivar[j, i, :num_pixels] = 1.0 / output['variance_electrons'][:, 0]

            # Get results for our flux-calibrated output file.
            cframe_observedflux[j, i, :num_pixels] = 1e17 * output['observed_flux'][:, 0]
            cframe_ivar[j, i, :num_pixels] = 1e-34 * output['flux_inverse_variance'][:, 0]

            # Fill brick arrays from the results.
            camera = output.meta['name']
            trueflux[camera][j][:] = 1e17 * output['observed_flux'][:, 0]
            noisyflux[camera][j][:] = 1e17 * (output['observed_flux'][:, 0] +
                                              output['flux_calibration'][:, 0] *
                                              output['random_noise_electrons'][:, 0])
            obsivar[camera][j][:] = 1e-34 * output['flux_inverse_variance'][:, 0]

            # Use the same noise realization in the cframe and frame, without any
            # additional noise from sky subtraction for now.
            frame_rand_noise[j, i, :num_pixels] = output['random_noise_electrons'][:, 0]
            cframe_rand_noise[j, i, :num_pixels] = 1e17 * (
                output['flux_calibration'][:, 0] * output['random_noise_electrons'][:, 0])

            # The sky output file represents a model fit to ~40 sky fibers.
            # We reduce the variance by a factor of 25 to account for this and
            # give the sky an independent (Gaussian) noise realization.
            sky_ivar[j, i, :num_pixels] = 25.0 / (
                output['variance_electrons'][:, 0] - output['num_source_electrons'][:, 0])
            sky_rand_noise[j, i, :num_pixels] = random_state.normal(
                scale=1.0 / np.sqrt(sky_ivar[j, i, :num_pixels]), size=num_pixels)

    armName = {"b": 0, "r": 1, "z": 2}
    for channel in 'brz':

        # Before writing, convert from counts/bin to counts/A (as in Pixsim output)
        # Quicksim Default:
        # FLUX - input spectrum resampled to this binning; no noise added [1e-17 erg/s/cm2/s/Ang]
        # COUNTS_OBJ - object counts in 0.5 Ang bin
        # COUNTS_SKY - sky counts in 0.5 Ang bin

        num_pixels = len(waves[channel])
        dwave = np.gradient(waves[channel])
        nobj[:, armName[channel], :num_pixels] /= dwave
        frame_rand_noise[:, armName[channel], :num_pixels] /= dwave
        nivar[:, armName[channel], :num_pixels] *= dwave**2
        nsky[:, armName[channel], :num_pixels] /= dwave
        sky_rand_noise[:, armName[channel], :num_pixels] /= dwave
        sky_ivar[:, armName[channel], :num_pixels] /= dwave**2

        # Now write the outputs in DESI standard file system. None of the output file can have more than args.n_fibers spectra

        # Looping over spectrograph
        for ii in range((args.nspec + args.nstart - 1) // args.n_fibers + 1):

            start = max(args.n_fibers * ii, args.nstart)  # first spectrum for a given spectrograph
            end = min(args.n_fibers * (ii + 1), nmax)  # last spectrum for the spectrograph

            if (args.spectrograph <= ii):
                camera = "{}{}".format(channel, ii)
                log.info("Writing files for channel:{}, spectrograph:{}, spectra:{} to {}".format(channel, ii, start, end))
                num_pixels = len(waves[channel])

                # Write frame file
                framefileName = lvmspec.io.findfile("frame", night, expid, camera)

                frame_flux = nobj[start:end, armName[channel], :num_pixels] + \
                    nsky[start:end, armName[channel], :num_pixels] + \
                    frame_rand_noise[start:end, armName[channel], :num_pixels]
                frame_ivar = nivar[start:end, armName[channel], :num_pixels]

                # required for slicing the resolution metric, resolusion matrix has (nspec, ndiag, wave)
                # for example if nstart =400, nspec=150: two spectrographs:
                # 400-499=> 0 spectrograph, 500-549 => 1
                sh1 = frame_flux.shape[0]

                if (args.nstart == start):
                    resol = resolution[channel][:sh1, :, :]
                else:
                    resol = resolution[channel][-sh1:, :, :]

                # must create lvmspec.Frame object
                frame = Frame(waves[channel], frame_flux, frame_ivar,
                              resolution_data=resol, spectrograph=ii,
                              fibermap=fibermap[start:end],
                              meta=dict(CAMERA=camera, FLAVOR=simspec.flavor))
                lvmspec.io.write_frame(framefileName, frame)

                framefilePath = lvmspec.io.findfile("frame", night, expid, camera)
                log.info("Wrote file {}".format(framefilePath))

                if args.frameonly or simspec.flavor == 'arc':
                    continue

                # Write cframe file
                cframeFileName = lvmspec.io.findfile("cframe", night, expid, camera)
                cframeFlux = cframe_observedflux[start:end, armName[channel], :num_pixels] + \
                    cframe_rand_noise[start:end, armName[channel], :num_pixels]
                cframeIvar = cframe_ivar[start:end, armName[channel], :num_pixels]

                # must create lvmspec.Frame object
                cframe = Frame(waves[channel], cframeFlux, cframeIvar,
                               resolution_data=resol, spectrograph=ii,
                               fibermap=fibermap[start:end],
                               meta=dict(CAMERA=camera, FLAVOR=simspec.flavor))
                lvmspec.io.frame.write_frame(cframeFileName, cframe)

                cframefilePath = lvmspec.io.findfile("cframe", night, expid, camera)
                log.info("Wrote file {}".format(cframefilePath))

                # Write sky file
                skyfileName = lvmspec.io.findfile("sky", night, expid, camera)
                skyflux = nsky[start:end, armName[channel], :num_pixels] + \
                    sky_rand_noise[start:end, armName[channel], :num_pixels]
                skyivar = sky_ivar[start:end, armName[channel], :num_pixels]
                skymask = np.zeros(skyflux.shape, dtype=np.uint32)

                # must create lvmspec.Sky object
                skymodel = SkyModel(waves[channel], skyflux, skyivar, skymask,
                                    header=dict(CAMERA=camera))
                lvmspec.io.sky.write_sky(skyfileName, skymodel)

                skyfilePath = lvmspec.io.findfile("sky", night, expid, camera)
                log.info("Wrote file {}".format(skyfilePath))

                # Write calib file
                calibVectorFile = lvmspec.io.findfile("calib", night, expid, camera)
                flux = cframe_observedflux[start:end, armName[channel], :num_pixels]
                phot = nobj[start:end, armName[channel], :num_pixels]
                calibration = np.zeros_like(phot)
                jj = (flux > 0)
                calibration[jj] = phot[jj] / flux[jj]

                # - TODO: what should calibivar be?
                # - For now, model it as the noise of combining ~10 spectra
                calibivar = 10 / cframe_ivar[start:end, armName[channel], :num_pixels]
                # mask=(1/calibivar>0).astype(int)??
                mask = np.zeros(calibration.shape, dtype=np.uint32)

                # write flux calibration
                fluxcalib = FluxCalib(waves[channel], calibration, calibivar, mask)
                write_flux_calibration(calibVectorFile, fluxcalib)

                calibfilePath = lvmspec.io.findfile("calib", night, expid, camera)
                log.info("Wrote file {}".format(calibfilePath))
