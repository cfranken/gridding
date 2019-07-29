import numpy as np
import matplotlib
import glob
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import h5py as h5
files = glob.glob("/net/fluo/data2/projects/MPI_FLUXCOM/GPP.RS.FP-ALL.MLM-ALL.METEO-NONE.4320_2160.8daily.*.nc")
allGPP = np.zeros((46,2160,4320,len(files)))
print(len(files))
counter = 0
for file in files:
	d = Dataset(file,"r")
	allGPP[:,:,:,counter] = d.variables["GPP"][:]
	counter = counter+1
	d.close()
	print(file)
mGPP = np.mean(allGPP, axis=3)
sGPP = np.std(allGPP, axis=3)
GPP = np.zeros((2160,4320,46))
GPP_std = np.zeros((2160,4320,46))
for i in range(46):
	GPP[:,:,i]=mGPP[i,::-1,:]
	GPP_std[:,:,i]=sGPP[i,::-1,:]
f1 = h5.File('/net/fluo/data2/projects/MPI_FLUXCOM/meanGPP.h5')
f2 = h5.File('/net/fluo/data2/projects/MPI_FLUXCOM/stdGPP.h5')
f1["GPP"]=GPP
f2["stdGPP"]=GPP_std
f1.close()
f2.close()
