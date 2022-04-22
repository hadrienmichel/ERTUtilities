import pybert as pb
from pybert.tdip.tdipdata import TDIP
tdip = TDIP("AgatheData\Project2_G7_ABMN_2.txt")
print(tdip)
tdip.setGates(dt=0.01, delay=0.01)
tdip.showDecay(ab=[1, 4])
tdip.fitDecays()
# %% some unfiltered data
tdip.showData("rhoa")
tdip.showData("m0")
tdip.showData("tau")
tdip.showData("fit")
# %% do some filtering
tdip.filter(fitmax=0.5, taumin=0.2, taumax=5, rmax=500, m0max=30)
tdip.showData("m0")
tdip.showData("tau")
