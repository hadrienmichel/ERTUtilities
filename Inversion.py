from pygimli.physics import ert# import pygimli as pg

mgr = ert.ERTManager(data='./HBHData/BERT_ERT_P1_err.dat', verbose=True)
# dataERT = mgr.data
# dataERT['k'] = ert.createGeometricFactors(dataERT, numerical=True)
# dataERT['rhoa'] = dataERT['k']*dataERT['r']
# dataERT.remove(dataERT['rhoa']<0)
# mgr.setData(dataERT)
data = mgr.data
data['k'] = ert.createGeometricFactors(data)
data['rhoa'] = data['r'] * data['k']
data.remove(data['rhoa'] < 1)
data.remove(data['rhoa'] > 5000)
print(data)
mod = mgr.invert()
mgr.showResultAndFit()
mgr.showModel(mod)
