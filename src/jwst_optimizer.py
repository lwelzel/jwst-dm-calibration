#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Practical for course 'Space Instrumentation',
TU Delft, The Netherlands
2022
By L. Welzel, Z. Burr
"""

import numpy as np
import matlab.engine
from os import getcwd
eng = matlab.engine.start_matlab()
eng.addpath(eng.genpath(getcwd()), nargout=0)

def objective():
    eng.JWST_Simulation(nargout=0)
    return
    

def optimize_jwst_mirror_segments():
    print(eng.sqrt(4.))
    return


if __name__ == "__main__":
    optimize_jwst_mirror_segments()