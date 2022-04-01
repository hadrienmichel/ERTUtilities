# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 17:07:25 2022

@author: tomde
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 15:07:14 2022

@author: tomde
"""

#Imports required modules

import pygimli as pg
import numpy as np
from os import listdir
from math import ceil
from shutil import rmtree
from os.path import isfile, join, exists
from matplotlib import pyplot as plt

plt.close('all')

## Directories initialization
directory = './HBHData'
save_dir = './HBHData'

# %% Baseline function for initilisation of data, inversion, plot & saving results
def baseline_inversion (baseline, paraDepth):  
    
    # %%% data initialisation & instance creation (ERTmanager & InversionFramework)
    
    data=baseline
    data['k'] = pg.physics.ert.createGeometricFactors(data, numerical=True)
    data['rhoa']= data['r']*data['k']
    print(data)
    mgr=pg.physics.ert.ERTManager(data)
    

    inv=mgr.fw.inv
    # data['err'] = mgr.estimateError(data, absoluteError=0.001, relativeError=2.0)
    mgr.checkData(data)
    # %%% Inversions of the baseline dataset  

    # First inversion using invert() to get an homogeneous model using med(rhoa) & to control the mesh parameters
    print(' ')
    print('Rough inversion of the baseline dataset')
    print(' ')
    #mesh = mgr.createMesh(data, quality=33.6, paraMaxCellSize=2,paraDepth=12,paraDX=0.6)
    mod_default=mgr.invert(data, lam=10,verbose=True,paraDepth=paraDepth/2,paraDX=1, paraMaxCellSize=20, quality=33.6, maxIter=1) # Force a first inversion to initialize the inversion
    #mod_default=mgr.invert(data, lam=10,verbose=True,mesh=mesh)
    # parameters initialization 
    #iv.setOptimizeLambda(True)
    inv.setRobustData(True)
    #inv.setBlockyModel(True)
    inv.setVerbose(True)
    #inv.stopAtChi1(True)
    inv.setLambda(10)
    inv.setDeltaPhiAbortPercent(1)
    
    # inversion using run ()
    print(' ')
    print('Inversion of the baseline dataset')
    print(' ')
    mod = inv.run()
    # mgr.showFit()
    
    # %%% plot of the computed vs mesured rhoa 

    x = np.linspace(min(data['rhoa']),max(data['rhoa']),100)
    y= x    

    figComp = plt.figure()
    ax2 = figComp.add_subplot(122)
    ax2.plot(data['rhoa'],inv.response(),'bo', label='Full Data')
    ax2.plot(x, y, 'k', label='Unity line')
    ax2.set_xlabel('Measured Apparent resistivities [Ohm.m]')
    ax2.set_ylabel('Computed Apparent resistivities [Ohm.m]')
    ax2.set_title('Measured VS Computed Rhoa')
    ax1 = figComp.add_subplot(121)
    relativeError = ((abs(data['rhoa']-inv.response()))/data['rhoa'] * 100)
    ax1.hist(relativeError, color='b')
    ax1.set_xlabel('Relative error [%]')
    ax1.set_ylabel('Number of data points')
    # %% data trim   
    if inv.relrms() > 10:
        indexRemove= relativeError > 10 # If error above 10% remove the point.
        data.remove(indexRemove)
        data.save(save_dir + '/dataTrimmed.ohm')
        indexKeep = [not(i) for i in indexRemove]
        relativeErrorBis = ((abs(data['rhoa']-np.asarray(inv.response())[indexKeep]))/data['rhoa'] * 100)
        ax1.hist(relativeErrorBis, color='r')
        ax2.plot(data['rhoa'], inv.response()[indexKeep], 'ro', label='Trimmed Dataset')
        plt.show()

        print(' ')
        print('The number of remaing data is %d' %len(data['rhoa']))
        print(' ')
        print(' ')
        print(indexRemove)
        print(' ')
        mod_default=mgr.invert(data, lam=10,verbose=True,paraDepth=paraDepth,paraDX=1, paraMaxCellSize=5, quality=33.6, maxIter=1)
        #mesh_2 = mgr.createMesh(data, quality=33.6, paraMaxCellSize=2,paraDepth=12,paraDX=0.6)
        #fig=pg.show(mesh)
        #mod_default=mgr.invert(data, lam=10,verbose=True,mesh=mesh_2)
        # Using L2 norm because outliers were removed
        inv.setRobustData(True)
        inv.setVerbose(True)
        inv.setLambda(10)
        mod = inv.run()

        mgr.showResult(mod, cMin=5, cMax=500, cMap="Spectral_r")
        # plt.plot(data['rhoa'],inv.response(),'go',label='Trimmed dataset')
        # plt.xlabel('Measured Apparent resistivities [Ohm.m]')
        # plt.ylabel('Computed Apparent resistivities [Ohm.m]')
        # plt.title('Measured VS Computed Rhoa')
        # plt.legend()
        # plt.show()
        
    else:
        plt.show()
        mgr.showResult(mod, cMin=5, cMax=500, cMap="Spectral_r")
        print(' ')
        print('Everything is gonna be alright')
        print(' ')
        indexRemove=False
        # Using L2 norm because they are no outliers
        inv.setRobustData(False)
        mod = inv.run()
    plt.savefig(join(save_dir,'Images/BG_compVS_meas'))
    # %%% saving results of the baseline inversion
    if exists(join(save_dir,'Base')):
        rmtree(join(save_dir,'Base'))
        mgr.saveResult(folder= join(save_dir,'Base'))
    else:
        mgr.saveResult(folder= join(save_dir,'Base'))
    
    
    ## DOI Computation 
    # mod_default_DOI=mgr.invert(data, lam=10,verbose=True,paraDepth=50,paraDX=0.6, paraMaxCellSize=2, quality=33.6)
    mean_rhoa = np.mean(data['rhoa'])
    
    ref1 = mod_default*0+ mean_rhoa*10
    # Ref model 1
    inv.setModel(ref1) 
    rm = mgr.inv.fop.regionManager()
    rm.setConstraintType(0) # 0 = reference model, 1 = first-order smoothing (default), 10 = both
    rm.fillConstraints(mgr.inv.fop.constraints())
    inv.setReferenceModel(ref1)
    mod_ref1 = inv.run()
    # Ref model 2
    ref2 = mod_default*0+ mean_rhoa/10
    inv.setModel(ref2)
    inv.setReferenceModel(ref2)
    mod_ref2 = inv.run()
    
    DOI = abs((mod_ref1 - mod_ref2)/(ref1-ref2))
    Norm_DOI = DOI/max(DOI)
    print(max(DOI))
    print(max(Norm_DOI))
    m = pg.Mesh(mgr.paraDomain)  
    m['Res']=mgr.paraModel(mod)
    m['Normalized DOI']= mgr.paraModel(Norm_DOI)
    m['DOI']= mgr.paraModel(DOI)
    m.exportVTK(join(save_dir, 'Inversion_t0'))  
    
    fig=pg.show(mgr.paraDomain, mod,cMin=5,cMax=25,label='Res [Ohm.m]')
    fig=pg.show(mgr.paraDomain,Norm_DOI,cMin=0,cMax=1,label='DOI')
    return mod, mgr ,inv, data, indexRemove

    
    
# # %% Timelapse function   
# def TLfunction(timestep_count,mod,mgr,inv,index_remove,data,timestep):

# # %%% data initialisation & instance creation (ERTmanager & InversionFramework)
#   # filtered_min=-1.5
#   # filtered_max=1
#   data_tl= pg.load(timestep)
#   if len(data_tl['a'])!=len (data['a']):
#       data_tl.remove(index_remove)
#   else:
#        pass
  
#   print(' ')
#   print('The number of remaing data is %d' %len(data_tl['a']))
#   print(' ')


#   data_tl['k']= pg.physics.ert.createGeometricFactors(data_tl, numerical=True)
#   data_tl['rhoa']= data_tl['r']*data_tl['k']
#   print(data_tl)
#   mgr_tl=pg.physics.ert.ERTManager(data_tl)
#   data_tl['err'] = mgr.estimateError(data_tl, absoluteError=0.001, relativeError=2.0)
#   mgr_tl.checkData(data_tl)
#   inv_tl=mgr_tl.fw.inv
# # %%% Visualisation of baseline & step rhoa 
#   # fig, ax = plt.subplots(3)
#   # fig.suptitle('Data visualization', fontsize=12)
#   # # Display the data as pseudosection
#   # ## Displays apparent resistivities (Ohm.m)
#   # mgr.showData(data,ax=ax[0],vals=data['rhoa'],label='Apparent resistivities')
#   # mgr_tl.showData(data,ax=ax[1],vals=data_tl['rhoa'],label='Apparent resistivities')
#   # ## Displays apparent resistivities differences (Ohm.m)
#   # mgr.showData(data,ax=ax[2],vals=data['rhoa']/data_tl['rhoa'],label='Apparent resistivities ratio')
#   ## Displays TL error computed with TL error model
#   # mgr_tl.showData(data,vals=data_tl['err'],label='Error')

# # %%% Inversions of the timestep dataset

# # %%%% First inversion using invert() to get an homogeneous model using med(rhoa) & to control the mesh parameters then inversion using run()
#   print('Rough inversion of the timestep %d' %timestep_count)
#   print(' ')
#   #mod_default_tl=mgr_tl.invert(data_tl,lam=10,verbose=True,paraDepth=12,paraDX=0.6, paraMaxCellSize=2, quality=33.6)
#   mod_default_tl=mgr_tl.invert(data, lam=10,verbose=True,paraDepth=12,paraDX=0.6, paraMaxCellSize=2, quality=33.6)
#   # pg.show(mgr.paraDomain,mod_default_tl)
 
# ## Inversion parameters control
#   #inv_tl.setOptimizeLambda(True)
#   inv_tl.setRobustData(False)
#   inv_tl.setVerbose(True)
#   #inv_tl.stopAtChi1(True)
#   inv_tl.setDeltaPhiAbortPercent(1)
#   inv_tl.setLambda(10)
# ## Independant inversion of the timestep dataset using run()
#   print(' ')
#   print('Basic inversion of the timestep %d' %timestep_count)
#   print(' ')  
#   mod_tl=inv_tl.run()
#   # mgr_tl.showFit()

#   # Shows the inversion results

#   fig, ax = plt.subplots(3,2)
#   fig.set_figheight(11)
#   fig.set_figwidth(9)
#   fig.suptitle('Basic inversion for timestep %d' %timestep_count, fontsize=12)
#   pg.show(mgr.paraDomain, mod,ax=ax[0,0],cMin=5,cMax=25 ,label='Res [Ohm.m]')
#   ax[0,0].text(30,34, 'Rrms = %f' %inv.relrms())
#   ax[0,0].text(30,32, 'Chi2= %f' %inv.chi2())
#   pg.show(mgr.paraDomain, mod_tl,ax=ax[0,1],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,1].text(30,34, 'Rrms = %f' %inv_tl.relrms())
#   ax[0,1].text(30,32, 'Chi2= %f' %inv_tl.chi2())
#   ## Shows the inversion results as a % change and absolute value difference
#   pg.show(mgr.paraDomain, 100*(mod-mod_tl)/mod,ax=ax[1,0], label='Ratio change [%]')

#   Abs_ch = mod - mod_tl
#   for i in range(len(Abs_ch)):
#          if (Abs_ch [i] >-1) and (Abs_ch[i]<1):
#              Abs_ch[i]=0
#          else:
#              pass
#   pg.show(mgr.paraDomain, Abs_ch,ax=ax[1,1], cMin=-5,cMax=5,cMap='seismic', label='Absolute difference [Ohm.m]') 
#   ax[2,0].hist(100*(mod-mod_tl)/mod,bins=ceil(np.sqrt(len(mod))),facecolor='b')
#   ax[2,0].set_xlabel('Ratio change [%]')
#   ax[2,0].set_ylabel('Occurence')
#   ax[2,1].hist(mod-mod_tl,bins=ceil(np.sqrt(len(mod))),facecolor='g')
#   ax[2,1].set_xlabel('Absolute difference [Ohm.m]')
#   ax[2,1].set_ylabel('Occurence')
#   fig.savefig(join(save_dir,'Images/basicinv_TL_%d' %timestep_count))
#   # %%%% Inversion using the baseline inverted model as starting model 

#   # We can change the starting model of the timelapse to the baseline model
#   inv_tl.setModel(mod) 
 
#   print(' ')
#   print('Inversion with the baseline as starting model of the timestep %d' %timestep_count)
#   print(' ')

#   mod_tl_advanced_start=inv_tl.run()
#   # Shows the inversion results
#   fig, ax = plt.subplots(3,2)
#   fig.set_figheight(11)
#   fig.set_figwidth(9)
#   fig.suptitle('Inversion with baseline as starting model for timestep %d' %timestep_count, fontsize=12)
#   pg.show(mgr.paraDomain, mod,ax=ax[0,0],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,0].text(30,34, 'Rrms = %f' %inv.relrms())
#   ax[0,0].text(30,32, 'Chi2= %f' %inv.chi2())
#   pg.show(mgr.paraDomain, mod_tl_advanced_start,ax=ax[0,1],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,1].text(30,34, 'Rrms = %f' %inv_tl.relrms())
#   ax[0,1].text(30,32, 'Chi2= %f' %inv_tl.chi2())
#   ## Shows the inversion results as a % change and absolute value difference
#   pg.show(mgr.paraDomain, 100*(mod-mod_tl_advanced_start)/mod,ax=ax[1,0],label='Ratio change [%]')
#   Abs_ch_adv_st = mod-mod_tl_advanced_start
#   for i in range(len(Abs_ch_adv_st)):
#          if (Abs_ch_adv_st [i] > -1) and (Abs_ch_adv_st[i]<1):
#              Abs_ch_adv_st[i]=0
#          else:
#              pass
#   pg.show(mgr.paraDomain,Abs_ch_adv_st,ax=ax[1,1], cMin=-5,cMax=5, cMap='seismic', label='Absolute difference [Ohm.m]')
#   ax[2,0].hist(100*(mod-mod_tl_advanced_start)/mod,bins=ceil(np.sqrt(len(mod))),facecolor='b')
#   ax[2,0].set_xlabel('Ratio change [%]')
#   ax[2,0].set_ylabel('Occurence')
#   ax[2,1].hist(mod-mod_tl_advanced_start,bins=ceil(np.sqrt(len(mod))),facecolor='g')
#   ax[2,1].set_xlabel('Absolute difference [Ohm.m]')
#   ax[2,1].set_ylabel('Occurence')
#   fig.savefig(join(save_dir,'Images/startinv_TL_%d' %timestep_count))
#   # %%%% Inversion using the baseline inverted model as starting & reference model
#   ## We can also set the reference model as the baseline model
#   rm_tl = mgr_tl.inv.fop.regionManager()
#   rm_tl.setConstraintType(0) # 0 = reference model, 1 = first-order smoothing (default), 10 = both
#   rm_tl.fillConstraints(mgr_tl.inv.fop.constraints())
#   inv_tl.setReferenceModel(mod)
#   inv_tl.setVerbose(True)
#   # inv_tl.stopAtChi1(True)
#   # inv.setRobustData(True)
#   #inv_tl.setOptimizeLambda(False)
#   # inv_tl.setLambda(1)
#   # # inv_tl.setOptimizeLambda(True)
#   print(' ')
#   print('Inversion with the baseline as starting & reference model of the timestep %d' %timestep_count)
#   print(' ')
#   mod_tl_advanced=inv_tl.run()

#   fig, ax = plt.subplots(3,2)
#   fig.set_figheight(11)
#   fig.set_figwidth(9)
#   fig.suptitle('Inversion with baseline as starting & reference model for timestep %d' %timestep_count, fontsize=12)
#   pg.show(mgr.paraDomain, mod,ax=ax[0,0],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,0].text(30,34, 'Rrms = %f' %inv.relrms())
#   ax[0,0].text(30,32, 'Chi2= %f' %inv.chi2())
#   pg.show(mgr.paraDomain, mod_tl_advanced,ax=ax[0,1],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,1].text(30,34, 'Rrms = %f' %inv_tl.relrms())
#   ax[0,1].text(30,32, 'Chi2= %f' %inv_tl.chi2())
 
#   pg.show(mgr.paraDomain, 100*(mod-mod_tl_advanced)/mod,ax=ax[1,0],label='Ratio change [%]')
#   Abs_ch_adv = mod-mod_tl_advanced
#   for i in range(len(Abs_ch_adv)):
#          if (Abs_ch_adv [i] >-1) and (Abs_ch_adv[i]<1):
#              Abs_ch_adv[i]=0
#          else:
#              pass
#   pg.show(mgr.paraDomain,Abs_ch_adv,ax=ax[1,1], cMin=-5,cMax=5, cMap='seismic', label='Absolute difference [Ohm.m]')
#   ax[2,0].hist(100*(mod-mod_tl_advanced)/mod,bins=ceil(np.sqrt(len(mod))),facecolor='b')
#   ax[2,0].set_xlabel('Ratio change [%]')
#   ax[2,0].set_ylabel('Occurence')
#   ax[2,1].hist(mod-mod_tl_advanced,bins=ceil(np.sqrt(len(mod))),facecolor='g')
#   ax[2,1].set_xlabel('Absolute difference [Ohm.m]')
#   ax[2,1].set_ylabel('Occurence')
#   fig.savefig(join(save_dir,'Images/ref_inv_TL_%d' %timestep_count))

#   # %%%% Ratio Inversion using the baseline inverted model as starting & reference model

#   inv_tl.setData(inv.response()/(data['rhoa']) * data_tl['rhoa'])
#   inv_tl.setRelativeError(0.2)
#   # inv_tl.setLambda(10)
#   print(' ')
#   print('Ratio inversion with the baseline as starting & reference model of the timestep %d' %timestep_count)
#   print(' ')
#   mod_tl_advanced_ratio = inv_tl.run()


#   fig, ax = plt.subplots(3,2)
#   fig.set_figheight(11)
#   fig.set_figwidth(9)
#   fig.suptitle('Ratio inversion with baseline as starting & reference model for timestep %d' %timestep_count, fontsize=12)
#   pg.show(mgr.paraDomain, mod,ax=ax[0,0],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,0].text(30,34, 'Rrms = %f' %inv.relrms())
#   ax[0,0].text(30,32, 'Chi2= %f' %inv.chi2())
#   pg.show(mgr.paraDomain, mod_tl_advanced_ratio,ax=ax[0,1],cMin=5,cMax=25,label='Res [Ohm.m]')
#   ax[0,1].text(30,34, 'Rrms = %f' %inv_tl.relrms())
#   ax[0,1].text(30,32, 'Chi2= %f' %inv_tl.chi2())
#   pg.show(mgr.paraDomain,100*(mod-mod_tl_advanced_ratio)/mod,ax=ax[1,0], label='Ratio change [%]')  
#   Abs_ch_adv_rat = mod-mod_tl_advanced_ratio
#   for i in range(len(Abs_ch_adv_rat)):
#          if (Abs_ch_adv_rat [i] > -0.01) and (Abs_ch_adv_rat[i]<0.01):
#              Abs_ch_adv_rat[i]=0
#          else:
#              pass
#   pg.show(mgr.paraDomain,Abs_ch_adv_rat,ax=ax[1,1], cMin=-5,cMax=5, cMap='seismic', label='Absolute difference [Ohm.m]')
#   ax[2,0].hist(100*(mod-mod_tl_advanced_ratio)/mod,bins=ceil(np.sqrt(len(mod))),facecolor='b')
#   ax[2,0].set_xlabel('Ratio change [%]')
#   ax[2,0].set_ylabel('Occurence')
#   ax[2,1].hist(mod-mod_tl_advanced_ratio,bins=ceil(np.sqrt(len(mod))),facecolor='g')
#   ax[2,1].set_xlabel('Absolute difference [Ohm.m]')
#   ax[2,1].set_ylabel('Occurence')
#   fig.savefig(join(save_dir,'Images/ratioinv_TL_%d' %timestep_count))
  
# # %%% saving results of the baseline inversion
#   if exists(join(save_dir,'TL_%d' %timestep_count)):
#       rmtree(join(save_dir,'TL_%d' %timestep_count))
#       mgr_tl.saveResult(folder= join(save_dir,'TL_%d' %timestep_count))
#   else:
#       mgr_tl.saveResult(folder= join(save_dir,'TL_%d' %timestep_count))
#   m_tl = pg.Mesh(mgr_tl.paraDomain)   
#   m_tl['m1']= mgr_tl.paraModel(Abs_ch)
#   m_tl['m2']= mgr_tl.paraModel( Abs_ch_adv_st)
#   m_tl['m3']= mgr_tl.paraModel(Abs_ch_adv)
#   m_tl['m4']= mgr_tl.paraModel(Abs_ch_adv_rat)
#   m_tl.exportVTK(join(save_dir, 'Inversion_t%d' %timestep_count))
# # %% Complet run of baseline_inversion & TL_function

## Initialization of the baseline dataset 
baseline = pg.load('./HBHData/BERT_ERT_P1_err.dat')
# mesh= pg.Mesh('C:/Users/tomde/Desktop/Bistade_monitoring/bistade_moni_d1/2D_pygimli/P2/data_bert_TL/mesh_BERT/mesh/mesh.bms')
mod,mgr,inv,data,index_remove= baseline_inversion (baseline, 100)

#Getting the timestep files in the chosen directory for the TL_function()
# if __name__ == "__main__":

#     filenames = [f for f in listdir(directory) if isfile(join(directory, f))]
#     timestep_count=1
#     for filename in filenames:
#         TLfunction(timestep_count,mod,mgr,inv,index_remove,data,timestep=join(directory,filename))
#         timestep_count+=1