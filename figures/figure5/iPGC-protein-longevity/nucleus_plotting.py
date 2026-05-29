from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv(Path(__file__).parent / "data" / "nucleus.csv")

# Normalize the dnd and nanos fluorophore levels to the globin
df["dnd_normalized"] = df["dnd"] / df["globin"]
df["nanos_normalized"] = df["nanos"] / df["globin"]

# Normalize to condition where the
df["dnd_normalized_to_full_mix"] = (
    df["dnd_normalized"] / df[df["condition"] == "full_mix"]["dnd_normalized"].mean()
)
df["nanos_normalized_to_full_mix"] = (
    df["nanos_normalized"]
    / df[df["condition"] == "full_mix"]["nanos_normalized"].mean()
)

fig, ax = plt.subplots(2, 2, figsize=(5, 10))
sns.lineplot(data=df, x="stage", y="dnd_normalized", ax=ax[0, 0], hue="condition")
sns.lineplot(data=df, x="stage", y="nanos_normalized", ax=ax[0, 1], hue="condition")
sns.lineplot(
    data=df, x="stage", y="dnd_normalized_to_full_mix", ax=ax[1, 0], hue="condition"
)
sns.lineplot(
    data=df, x="stage", y="nanos_normalized_to_full_mix", ax=ax[1, 1], hue="condition"
)
fig.savefig(Path(__file__).parent / "output" / "dnd_normalized.png")
