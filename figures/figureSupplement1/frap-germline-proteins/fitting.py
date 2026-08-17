import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


def one_phase_exp_model(t, I0, Iinf, tau):
    return Iinf - (Iinf - I0) * np.exp(-t / tau)


def one_phase_exp_model_fixed_zero(t, Iinf, tau):
    return Iinf - Iinf * np.exp(-t / tau)


def two_phase_exp_model(t, I0, Iinf, tau1, tau2, A):
    return Iinf - (Iinf - I0) * (A * np.exp(-t / tau1) + (1 - A) * np.exp(-t / tau2))


def _early_time_sigma(t, offset=0.5, power=1.0):
    """
    Build sigma values for curve_fit to prioritize early recovery points.
    Smaller sigma => larger optimization weight in curve_fit.
    """
    if offset <= 0:
        raise ValueError("offset must be > 0")
    if power <= 0:
        raise ValueError("power must be > 0")

    t = np.asarray(t, dtype=float)
    t_relative = t - np.nanmin(t)
    weights = 1.0 / np.power(t_relative + offset, power)
    sigma = 1.0 / weights
    return sigma


def recovery_fit(
    time,
    intensity,
    model="one_phase",
    low_limit_index=None,
    high_limit_index=None,
    weighted=False,
    weight_offset=0.5,
    weight_power=1.0,
):
    t = time[low_limit_index:high_limit_index]
    y = intensity[low_limit_index:high_limit_index]
    sigma = (
        _early_time_sigma(t, offset=weight_offset, power=weight_power)
        if weighted
        else None
    )

    if model == "one_phase":
        popt, _ = curve_fit(
            one_phase_exp_model,
            t,
            y,  # , bounds=(0, time.max() * 10)
            sigma=sigma,
            absolute_sigma=False,
        )
        y_pred = one_phase_exp_model(t, *popt)
        rsquared = r2_score(y, y_pred)
        params = {"tau": popt[2], "I0": popt[0], "Iinf": popt[1], "r_squared": rsquared}

        fit_values = {
            "time": t,
            "fitted_intensity": y_pred,
        }

        return params, fit_values

    elif model == "one_phase_fixed_zero":
        popt, _ = curve_fit(
            one_phase_exp_model_fixed_zero,
            t,
            y,
            sigma=sigma,
            absolute_sigma=False,
        )
        y_pred = one_phase_exp_model_fixed_zero(t, *popt)
        rsquared = r2_score(y, y_pred)
        params = {"tau": popt[1], "Iinf": popt[0], "r_squared": rsquared}

        fit_values = {
            "time": t,
            "fitted_intensity": y_pred,
        }

        return params, fit_values

    elif model == "two_phase":
        popt, _ = curve_fit(
            two_phase_exp_model,
            t,
            y,  # , bounds=(0, time.max() * 10)
            sigma=sigma,
            absolute_sigma=False,
        )
        y_pred = two_phase_exp_model(t, *popt)

        params = {
            "tau1": popt[2],
            "tau2": popt[3],
            "A": popt[4],
            "I0": popt[0],
            "Iinf": popt[1],
            "r_squared": r2_score(y, y_pred),
        }

        fit_values = {
            "time": t,
            "fitted_intensity": y_pred,
        }

        return params, fit_values
    else:
        raise ValueError(
            "Model must be 'one_phase', 'one_phase_fixed_zero', or 'two_phase'"
        )
