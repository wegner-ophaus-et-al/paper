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


def get_representative_ids(df, feature, rep_min=40, rep_max=60):
    """
    For a given feature returns the IDs of samples that are representative all conditions, meaning they are around the median of the distribution for that feature. This can be used to select samples for visualizations that are representative of the overall trends in the data.
    params:
    df: DataFrame containing results of an expeiment
    feature: Feature, meaning column, of the dataframe in the
    rep_min: lower percentile threshold for selecting representative samples
    rep_max: upper percentile threshold for selecting representative samples
    """
    condition_dict = {}
    for condition in df["condition"].unique():
        representative_ids = []
        for rna_type in df["construct_id"].unique():
            sub_df = df.query("construct_id == @rna_type and condition == @condition")
            min_val = np.percentile(sub_df[feature].values, rep_min)
            max_val = np.percentile(sub_df[feature].values, rep_max)
            if len(representative_ids) == 0:
                print(sub_df.head(10))
                print(sub_df.columns)
                representative_ids = sub_df.query(f"@min_val <= {feature} <= @max_val")[
                    "uid"
                ].values
            else:
                representative_ids = np.intersect1d(
                    representative_ids,
                    sub_df.query(f"@min_val <= {feature} <= @max_val")["UID"].values,
                )

            condition_dict[condition] = representative_ids

    return condition_dict
