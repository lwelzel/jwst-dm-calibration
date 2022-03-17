# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 12:19:50 2022

@author: Zach
"""
import numpy as np

with open("JWST_aligmnent_results.txt", "r") as f:
    lines = f.readlines()
    
outputs = np.zeros((450,3))
for i in range(450):
    line = lines[2+i*41]
    line = line[1:-3].split()
    for j in range(len(line)):
        outputs[i, j] = float(line[j])

