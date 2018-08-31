from numpy import array
import dadi
import my_demo_modelsC
import pylab

dd = dadi.Misc.make_data_dict("CriolloCuraray.dadi")
data = dadi.Spectrum.from_data_dict(dd, ['curaray', 'criollo'],[10,8])
ns = data.sample_sizes
data_f = data.fold()
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

params = array([  1.30262462e+00,   5.07789590e+00,   2.44631180e-02, 1.99996373e+00,   1.09330231e-01,   3.49366683e-01, 1.94998738e-03,   1.99998470e+00])
upper_bound = [5 , 6, 2, 2, 0.2, 5, 2, 2]
lower_bound = [1e-3, 1e-3, 1e-3, 1e-3, 0.001, 0, 0, 0]

model = func_ex(params, ns, pts_l)
model_f = model.fold()
ll_model = dadi.Inference.ll_multinom(model_f, data_f)

p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)

popt = dadi.Inference.optimize_log_fmin(p0, data_f, func_ex, pts_l,
                                lower_bound=lower_bound,
                                upper_bound=upper_bound,
                                verbose=len(params))

model = func_ex(popt, ns, pts_l)
model.mask[1,1] = True
model.mask[1,2] = True
model.mask[2,1] = True
model.mask[2,2] = True
model.mask[0,2] = True
model.mask[2,0] = True
model.mask[1,0] = True
model.mask[0,1] = True
model_f = model.fold()

ll_optF = dadi.Inference.ll_multinom(model_f, data_f)

thetaF = dadi.Inference.optimal_sfs_scaling(model_f, data_f)

expected_folded = model_f*thetaF
residAnscombe = dadi.Inference.Anscombe_Poisson_residual(expected_folded,data_f)
residLinear = dadi.Inference.linear_Poisson_residual(expected_folded,data_f)
