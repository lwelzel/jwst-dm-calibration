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
eng = matlab.engine.start_matlab()
eng.JWST_Simulation(nargout=0)

def objective():
    return
    

def optimize_jwst_mirror_segments():
    return


if __name__ == "__main__":
    optimize_jwst_mirror_segments()