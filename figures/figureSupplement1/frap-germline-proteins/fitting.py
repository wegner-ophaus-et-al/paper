import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


def one_phase_exp_model(t, I0, Iinf, tau):
    return Iinf - (Iinf - I0) * np.exp(-t / tau)


def one_phase_exp_model_fixed_zero(t, Iinf, tau):
    return Iinf - Iinf * np.exp(-t / tau)


def two_phase_exp_model(t, I0, Iinf, tau1, tau2, A):
    return Iinf - (Iinf - I0) * (A * np.exp(-t / tau1) + (1 - A) * np.exp(-t / tau2))


def recovery_fit(
    time, intensity, model="one_phase", low_limit_index=None, high_limit_index=None
):
    if model == "one_phase":
        popt, _ = curve_fit(
            one_phase_exp_model,
            time[low_limit_index:high_limit_index],
            intensity[
                low_limit_index:high_limit_index
            ],  # , bounds=(0, time.max() * 10)
        )
        y_pred = one_phase_exp_model(time[low_limit_index:high_limit_index], *popt)
        rsquared = r2_score(intensity[low_limit_index:high_limit_index], y_pred)
        params = {"tau": popt[2], "I0": popt[0], "Iinf": popt[1], "r_squared": rsquared}

        fit_values = {
            "time": time[low_limit_index:high_limit_index],
            "fitted_intensity": y_pred,
        }

        return params, fit_values

    elif model == "one_phase_fixed_zero":
        popt, _ = curve_fit(
            one_phase_exp_model_fixed_zero,
            time[low_limit_index:high_limit_index],
            intensity[low_limit_index:high_limit_index],
        )
        y_pred = one_phase_exp_model_fixed_zero(
            time[low_limit_index:high_limit_index], *popt
        )
        rsquared = r2_score(intensity[low_limit_index:high_limit_index], y_pred)
        params = {"tau": popt[1], "Iinf": popt[0], "r_squared": rsquared}

        fit_values = {
            "time": time[low_limit_index:high_limit_index],
            "fitted_intensity": y_pred,
        }

        return params, fit_values

    elif model == "two_phase":
        popt, _ = curve_fit(
            two_phase_exp_model,
            time[low_limit_index:high_limit_index],
            intensity[
                low_limit_index:high_limit_index
            ],  # , bounds=(0, time.max() * 10)
        )
        y_pred = two_phase_exp_model(time[low_limit_index:high_limit_index], *popt)

        params = {
            "tau1": popt[2],
            "tau2": popt[3],
            "A": popt[4],
            "I0": popt[0],
            "Iinf": popt[1],
            "r_squared": r2_score(intensity[low_limit_index:high_limit_index], y_pred),
        }

        fit_values = {
            "time": time[low_limit_index:high_limit_index],
            "fitted_intensity": y_pred,
        }

        return params, fit_values
    else:
        raise ValueError("Model must be 'one_phase' or 'two_phase'")
