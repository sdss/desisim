'''
lvmsim.specsim
===============

DESI wrapper functions for external specsim classes.
'''

from __future__ import absolute_import, division, print_function
import numpy as np
from astropy.units import Quantity
from specsim.config import Configuration
import lvmutil.log
import specsim.simulator
import six
import difflib
import re


# - Cached simulators, keyed by config string
_simulators = dict()

# - Cached defaults after loading a new simulator, to be used to reset a
# - simulator back to a reference state before returning it as a cached copy
_simdefaults = dict()


log = lvmutil.log.get_logger()


PARAM_MAP = {'area.M1_diameter': 'instrument.primary_mirror_diameter',
             'fibers.diameter_um': 'instrument.fiber_diameter',
             'area.obscuration_diameter': 'instrument.obscuration_diameter',
             'area.M2_support_width': 'instrument.support_width',
             'ccd.*.readnoise': 'instrument.cameras.*.read_noise',
             'ccd.*.darkcurrent': 'instrument.cameras.*.dark_current',
             'ccd.*.gain': 'instrument.cameras.*.gain',
             'exptime': 'observation.exposure_time'}


def get_simulator(config='lvm', num_fibers=1, params=None):
    '''
    returns new or cached specsim.simulator.Simulator object

    Also adds placeholder for BGS fiberloss if that isn't already in the config
    '''
    if isinstance(config, Configuration):
        w = config.wavelength
        wavehash = (np.min(w), np.max(w), len(w))
        key = (config.name, wavehash, num_fibers)
    else:
        key = (config, num_fibers)

    parammsg = 'with telescope {0}'.format(params['telescope']) if params else ''
    msg = '{0} Simulator for {1} {2}'.format(key, config, parammsg)

    if key in _simulators:
        log.debug('Returning cached {0}'.format(msg))
        qsim = _simulators[key]
        defaults = _simdefaults[key]
        qsim.source.focal_xy = defaults['focal_xy']
        qsim.atmosphere.airmass = defaults['airmass']
        qsim.observation.exposure_time = defaults['exposure_time']
        qsim.atmosphere.moon.moon_phase = defaults['moon_phase']
        qsim.atmosphere.moon.separation_angle = defaults['moon_angle']
        qsim.atmosphere.moon.moon_zenith = defaults['moon_zenith']

    else:
        log.debug('Creating new {0}'.format(msg))

        # - New config; create Simulator object
        qsim = specsim.simulator.Simulator(config, num_fibers, params=params)

        # - Cache defaults to reset back to original state later
        defaults = dict()
        defaults['focal_xy'] = qsim.source.focal_xy
        defaults['airmass'] = qsim.atmosphere.airmass
        defaults['exposure_time'] = qsim.observation.exposure_time
        defaults['moon_phase'] = qsim.atmosphere.moon.moon_phase
        defaults['moon_angle'] = qsim.atmosphere.moon.separation_angle
        defaults['moon_zenith'] = qsim.atmosphere.moon.moon_zenith

        _simulators[key] = qsim
        _simdefaults[key] = defaults

        # update the parameters
        if params:
            update_params(qsim, params)

    return qsim


def _parse_param(name):
    ''' Parse the PARAM_MAP value into group/parameter names

    Takes a PARAM_MAP key/value and split it into groups and parameters
    E.g., instrument.primary_mirror_diameter becomes
    (['instrument'], 'primary_mirror_diameter') or
    'instrument.cameras.*.read_noise' becomes (['instrument', 'cameras', '*'], 'read_noise')

    Parameters:
        name (str):
            The name of the parameter in the mapping

    Returns:
        A tuple containing a list of groups names, and a parameter name

    '''

    # assume last dot is the parameter, split into groups + parameter
    groups, par = name.rsplit('.', 1)

    # check for nested groups
    dotcount = groups.count('.')
    grp = groups.split('.')

    return grp, par


def _set_parameter(base, param_name, param_value):
    ''' Sets a new parameter value '''

    base_obj = _retrieve_object(base, param_name)
    if base_obj:
        newq = Quantity(param_value, base_obj.unit)
        base.__setattr__(param_name, newq)
    return base


def _retrieve_object(base, name):
    ''' Retrieve an object from the simulator '''

    # loop over nested groups
    if isinstance(name, list):
        if '*' in name:
            name.remove('*')
        for n in name:
            base = _retrieve_object(base, n)
        return base

    # retrieve the group object
    base_has_name = hasattr(base, name)
    if base_has_name:
        return base.__getattribute__(name)
    else:
        return None


def _update_obj(sim, key, origkey, value, wild=None):
    ''' Update the simulator object values

    Parameters:
        sim (inst):
            the Simulator object
        key (str):
            The simulator mapping key (PARAM_MAP value)
        origkey (str):
            The original full parameter name (e.g. area.M1_diameter)
        value (int|str|float):
            The new parameter value
        wild (str):
            The original wildcard PARAM_MAP key (e.g. ccd.*.readnoise)
    '''

    # parse the key into groups and parameters
    groups, par = _parse_param(key)
    # get the group out of the simulator
    grpobj = _retrieve_object(sim, groups)

    # if it's a wildcard parameter, identify the name of the wildcard, e.g. ccd cameras
    if wild:
        name = [i[-1] for i in difflib.ndiff(wild, origkey) if '+' in i][0]
        grpobj = [g for g in grpobj if g.name == name][0]

    if grpobj:
        # set the parameter
        grpobj = _set_parameter(grpobj, par, value)


def update_params(sim, params, base=None):
    ''' Updates the simulator parameters

    Updates the simulator parameteres with those from the LVMMODEL yaml file
    using the PARAM_MAP dictionary as a mapping.  Parameters to be
    updated that have a corresponding value into the built-in specsim config
    must have a mapping available.

    Loops over all parameters in the yaml file, checks if a mapping is available and attempts
    to update the build in config with the new value.

    Parameters:
        sim (inst):
            The Simualator object
        params (dict):
            A dictionary of parameters
        base (str):
            the base part of a parameter name
    '''

    # loop over the parameters
    for key, values in params.items():

        # if nested, go deeper
        if isinstance(values, dict):
            # join the base group with the param (e.g. area.M1_diameter)
            newbase = base + key + '.' if base else key + '.'
            update_params(sim, values, base=newbase)
        else:
            # build the full parameter name
            param = base + key if base else key
            # check if parameter in PARAM_MAP
            if param in PARAM_MAP:
                # get the specsim key
                simpar_key = PARAM_MAP[param]
                # update a simulator object
                _update_obj(sim, simpar_key, param, values)
            else:
                # check for wildcard params
                wildparams = [p for p in PARAM_MAP if '*' in p]
                if wildparams:
                    wildmatch = [w for w in wildparams if re.match(w.replace('*', '[a-z]'), param)]
                    if wildmatch:
                        simpar_key = PARAM_MAP[wildmatch[0]]
                        # update a simulator object
                        _update_obj(sim, simpar_key, param, values, wild=wildmatch[0])



