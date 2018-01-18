from __future__ import absolute_import, division, print_function

import unittest

def test_suite():
    """Returns unittest.TestSuite of lvmsim tests for use by setup.py"""

    #- DEBUG Travis test failures
    # return unittest.defaultTestLoader.loadTestsFromNames([
    #     # 'lvmsim.test.test_batch',      #- OK
    #     # 'lvmsim.test.test_io',         #- OK
    #     # 'lvmsim.test.test_obs',        #- OK
    #     'lvmsim.test.test_pixsim',
    #     # 'lvmsim.test.test_quickcat',   #- OK
    #     # 'lvmsim.test.test_targets',    #- OK
    #     # 'lvmsim.test.test_templates',  #- OK
    #     # 'lvmsim.test.test_top_level',  #- OK
    #     ])
    #- DEBUG Travis test failures

    from os.path import dirname
    lvmsim_dir = dirname(dirname(__file__))
    print(lvmsim_dir)
    return unittest.defaultTestLoader.discover(lvmsim_dir,
        top_level_dir=dirname(lvmsim_dir))
