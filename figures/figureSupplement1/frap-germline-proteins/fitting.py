import numpy as np
from scipy.optimize import curve_fit


def exp_model(t, I0, Iinf, tau):
    return Iinf - (Iinf - I0) * np.exp(-t / tau)


# t: time points, y: normalized recovery
p0 = [y[0], y[-1], np.median(t)]
pars, cov = curve_fit(exp_model, t, y, p0=p0, bounds=([0, 0, 0], [1.5, 2, np.inf]))

I0, Iinf, tau = pars
tau_sd = np.sqrt(np.diag(cov))[2]
