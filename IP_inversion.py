# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:57:24 2022

@author: tomde
"""


import pybert as pb
import pygimli as pg


tdip=pb.TDIPdata(filename='C:/Users/tomde/Desktop/Bistade_monitoring/bistade_moni_d1/Second_phase/IP/TL/L2/data_bert_TL/BERT_T3_L2.ohm', verbose=True)
data = tdip.data
mgr = pb.ERTManager(data = tdip.data)


#tdip.showRhoa()


mesh = mgr.createMesh(data, paraDepth = 15, quality = 33.6, paraMaxCellSize = 1, paraDX = 0.6)


# # Invert resistivity

tdip.invertRhoa(mesh=mesh,lam = 10, robustData= False, verbose = True)
tdip.showResistivity(cMap='jet')



# # # Invert chargeability

tdip.invertMa(nr = 0, lam= 10, ma = tdip.data['ip']*0.001)
tdip.showChargeability(cMap='jet',cMin = 0, cMax = 25)


### Export VTK 
tdip.pd.addData("Res", tdip.res)
tdip.pd.addData("Coverage", tdip.coverage)
tdip.pd.addData("IP", tdip.m*1000)
tdip.pd.exportVTK('C:/Users/tomde/Desktop/Bistade_monitoring/bistade_moni_d1/Second_phase/IP/TL/L2/data_bert_TL/IP_inv_T3_L2')