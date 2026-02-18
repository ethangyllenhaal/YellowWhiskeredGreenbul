import moments
import numpy as np

def model_func(params, ns):
	nu_1, nu_2, t1, nu11, nu12, m1_12, m1_21, nu12_1, nu12_2, t2, nu21, nu22, nu23, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32 = params
	_Nanc_size = 1.0  # This value can be used in splits with fractions
	sts = moments.LinearSystem_1D.steady_state_1D(np.sum(ns))
	fs = moments.Spectrum(sts)
	fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
	nu2_func = lambda t: nu_2 * (nu12 / nu_2) ** (t / t1)
	migs = np.array([[0, m1_12], [m1_21, 0]])
	fs.integrate(tf=t1, Npop=lambda t: [nu11, nu2_func(t)], m=migs, dt_fac=0.01)
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	nu1_func = lambda t: nu11 + (nu21 - nu11) * (t / t2)
	nu2_func = lambda t: nu12_1 * (nu22 / nu12_1) ** (t / t2)
	nu3_func = lambda t: nu12_2 * (nu23 / nu12_2) ** (t / t2)
	migs = np.array([[0, m2_12, m2_13], [m2_21, 0, m2_23], [m2_31, m2_32, 0]])
	fs.integrate(tf=t2, Npop=lambda t: [nu1_func(t), nu2_func(t), nu3_func(t)], m=migs, dt_fac=0.01)
	return fs

data = moments.Spectrum.from_file('/lustre/scratch/egyllenh/redo_greenbul/gadma/easysfs_noZ_thin1k/dadi/East-Central-West.sfs')
data = data.project([34, 38, 8])
data = np.transpose(data, [2, 1, 0])
data.pop_ids = ['East', 'Central', 'West']
ns = data.sample_sizes

p0 = [0.9908792926332146, 0.10431730074038632, 0.8396418992266634, 0.18548199705205884, 3.3601260507000252, 4.062715476271032, 0.0, 10.0, 0.1, 0.3190315514934076, 10.0, 10.0, 6.003154491162586, 0.30822397128678836, 0.050428132107367685, 0.12077831516790752, 0.09477780271727647, 0.001, 0.6998904486196827]
lower_bound = [0.1, 0.1, 1e-15, 0.1, 0.1, 0.0, 0.0, 0.1, 0.1, 1e-15, 0.1, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
upper_bound = [10.0, 10.0, 5.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 5.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
model = model_func(p0, ns)
ll_model = moments.Inference.ll_multinom(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

Nanc = 564112.5223388347
mu = 4.6e-09
L = 632399
theta0 = 4 * mu * L
Nanc = int(theta / theta0)
print('Size of ancestral population: {0}'.format(Nanc))


plot_ns = [4 for _ in ns]  # small sizes for fast drawing
gen_mod = moments.ModelPlot.generate_model(model_func,
                                           p0, plot_ns)
moments.ModelPlot.plot_model(gen_mod,
                             save_file='model_from_GADMA.svg',
                             fig_title='Demographic model from GADMA',
                             draw_scale=True,
                             pop_labels=['West', 'Central', 'East'],
                             nref=564112,
                             gen_time=3.27,
                             gen_time_units='years',
                             reverse_timeline=True)
