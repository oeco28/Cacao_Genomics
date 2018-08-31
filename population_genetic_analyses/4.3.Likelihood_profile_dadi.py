from numpy import array
import dadi
import my_demo_modelsC
import pylab
import numpy

dd = dadi.Misc.make_data_dict("cacao.dadi")
data = dadi.Spectrum.from_data_dict(dd, ['curaray', 'criollo'],[10,8])
ns = data.sample_sizes
pts_l = [40,50,60]

data.mask[1,1] = True
data.mask[1,2] = True
data.mask[2,1] = True
data.mask[2,2] = True
data.mask[0,2] = True
data.mask[2,0] = True
data.mask[1,0] = True
data.mask[0,1] = True

data_f = data.fold()

func = my_demo_modelsC.split_grow
func_ex = dadi.Numerics.make_extrap_log_func(func)

params = array([ 0.902, 6.00, 0.0224, 1.04, 0.179, 0.3073, 0.0314, 2.092 ])

my_nuCur0 = numpy.arange(0.1,10.1,0.5)
my_TCur = numpy.arange(0.5,25.5, 1)
my_nuCurF = numpy.arange(0.1,1.3,0.1)
my_nuCriF = numpy.arange(0.1,2,0.1)
#my_s = numpy.arange(0.1,1,0.1)
my_T2 = numpy.arange(0.1,2,0.1)
#my_m12 = numpy.arange(0.01,1,0.01)
#my_m21 = numpy.arange(0.2,4,0.2)

ll_prof = numpy.zeros(shape=(len(my_nuCur0),len(my_TCur),len(my_nuCurF),len(my_nuCriF),len(my_T2)))
Thetas = numpy.zeros(shape=(len(my_nuCur0),len(my_TCur),len(my_nuCurF),len(my_nuCriF),len(my_T2)))
my_params = numpy.zeros(1625184)

for i in range(len(my_nuCur0)):
    for j in range(len(my_TCur)):
        for l in range(len(my_nuCurF)):
            for m in range(len(my_nuCriF)):
                for n in range(len(my_T2)):
                    pts_l = [40,50,60]
                    ns = data.sample_sizes
                    func_split = my_demo_modelsC.split_grow
                    params = array([my_nuCur0[i], my_TCur[j], my_nuCurF[l], my_nuCriF[m], 0.179, my_T2[n], 0.0314, 2.092])
                    func_ex = dadi.Numerics.make_extrap_log_func(func_split)
                    my_model_g = func_ex(params, ns, pts_l)
                    ll_my_model_g = dadi.Inference.ll_multinom(my_model_g, data_f)
                    #theta_growth2 = dadi.Inference.optimal_sfs_scaling(my_model_g, data_f)
                    ll_prof[i,j,l,m,n] = ll_my_model_g
                    #Thetas[i,j,l,m,n] = theta_growth2
from numpy import *
from pylab import load
from pylab import save
numpy.savetxt('my_ll_prof_base.txt', ll_prof)
