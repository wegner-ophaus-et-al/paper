import numpy as np
from scipy import stats


def parametric(data_set1, data_set2):
    set1_normal = stats.normaltest(data_set1).pvalue > 0.05
    set2_normal = stats.normaltest(data_set2).pvalue > 0.05

    return min(set1_normal, set2_normal)


def statistical_analysis(data_set1: list, data_set2: list):
    if parametric(data_set1, data_set2):
        t_statistic, p_value = stats.ttest_ind(data_set1, data_set2)
        test_type = "t-test"
    else:
        t_statistic, p_value = stats.mannwhitneyu(data_set1, data_set2)
        test_type = "Mann-Whitney U test"
    return test_type, t_statistic, p_value


def get_stars(p_value):
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


def print_nested(d, indent=0):
    for key, value in d.items():
        prefix = "  " * indent
        if isinstance(value, dict):
            print(f"{prefix}{key}:")
            print_nested(value, indent + 1)
        elif (
            isinstance(value, list)
            and value
            and all(isinstance(i, dict) for i in value)
        ):
            print(f"{prefix}{key}: list ({len(value)})")
            print_nested(value[0], indent + 1)
        elif isinstance(value, np.ndarray):
            print(f"{prefix}{key}: ndarray {value.shape}")
        else:
            print(f"{prefix}{key}: {type(value).__name__}")
