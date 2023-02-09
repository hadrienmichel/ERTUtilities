# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:22:45 2021

@author: tomde
"""

import numpy as np

from Utilities import ArrayToXML, SampleReciprocals, SortArray, CreateReciprocals

points = np.loadtxt('C:/Users/tomde/OneDrive/Bureau/github/ERTProtocolsCreator/GD7.txt',dtype=int)

reciprocals = CreateReciprocals(points)
sampled_rec= SampleReciprocals(points, 0.1)

ArrayToXML(points)
ArrayToXML(reciprocals)
ArrayToXML(sampled_rec)