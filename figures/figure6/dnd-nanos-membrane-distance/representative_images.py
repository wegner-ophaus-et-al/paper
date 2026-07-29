import numpy as np


def export_representative_uids(
    df, features, output_path, low=40, high=60, uid_col="uid", condition_col="condition"
):
    """Write UIDs of samples whose feature value falls within [low, high]
    percentile of their condition group, for each feature, to a text file."""
    lines = []
    for feature in features:
        lines.append(f"Feature: {feature}")
        for condition, group in df.groupby(condition_col):
            values = group[feature].dropna()
            if values.empty:
                lines.append(f"  {condition}: none")
                continue

            p_low, p_high = np.percentile(values, [low, high])
            mask = values.between(p_low, p_high)
            uids = group.loc[mask.index[mask], uid_col].tolist()

            lines.append(f"  {condition}: {', '.join(uids) if uids else 'none'}")
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))
