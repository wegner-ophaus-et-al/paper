from scipy import stats


def parametric(data_set1, data_set2):
    set1_normal = stats.normaltest(data_set1).pvalue > 0.05
    set2_normal = stats.normaltest(data_set2).pvalue > 0.05

    return min(set1_normal, set2_normal)


def statistical_analysis(data_set1, data_set2):
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
