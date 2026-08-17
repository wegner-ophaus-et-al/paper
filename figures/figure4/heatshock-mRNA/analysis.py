import marimo

__generated_with = "0.23.10"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    from pathlib import Path
    import tifffile as tiff
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from utils import get_confocal_pixel_size, compute_masks, statistical_analysis, get_stat_stars

    return (
        Path,
        compute_masks,
        get_confocal_pixel_size,
        get_stat_stars,
        np,
        pd,
        plt,
        sns,
        statistical_analysis,
        tiff,
    )


@app.cell
def _(Path):
    root = Path("/Volumes/Kur/paper_data_sorted/heatshock_rna/data")
    repeat_dirs = [pth for pth in root.iterdir() if pth.is_dir() and not pth.name.startswith('.')]
    list_of_sample_paths = []

    for repeat_dir in repeat_dirs:
        list_of_sample_paths.extend([pth for pth in repeat_dir.iterdir() if pth.is_dir() and not pth.name.startswith('.')])
    print(f'Found {len(list_of_sample_paths)} samples across {len(list_of_sample_paths)} repeats.')
    return list_of_sample_paths, root


@app.cell
def _():
    heat_shock_palette = {"control": "#888888", "heatshock": "#b04a46"}
    return (heat_shock_palette,)


@app.cell
def _(compute_masks, get_confocal_pixel_size, list_of_sample_paths, np, tiff):
    measurement_dicts = []
    for sample_path in list_of_sample_paths:
        sample_dict = {}
        splits = sample_path.name.split('__')
        if len(splits) == 2:
            uid, name = splits
            date = '2026-01-14'  # print(f'Processing sample: {sample_path.name}')
        elif len(splits) == 3:
            uid, name, date = splits
        else:
            print(f'Unexpected sample name format: {sample_path.name}')
            continue
        proto_condition, sample_number = name.split('-', 1)
        _condition = None
        if 'soma' in proto_condition:
            continue
        elif 'heatshock' in proto_condition:
            _condition = 'heatshock'
        elif 'control' in proto_condition:
            _condition = 'control'
        else:
            print(f'Unexpected proto_condition format: {proto_condition}')
            continue
        sample_dict.update({'sample_path': sample_path, 'uid': uid, 'condition': _condition, 'repeat_date': date, 'sample_number': sample_number, 'is_heatshock': True if 'heat' in _condition.lower() else False})
        print(f'Processing sample: {sample_path.name}')
        pixel_size = get_confocal_pixel_size(sample_path)
        img_rna_actin_path = sample_path / 'raw' / 'bact_rna.tiff'
        img_granule_path = sample_path / 'raw' / 'granule.tiff'
        img_rna_nanos_path = sample_path / 'raw' / 'nos_rna.tiff'
        img_nucleus_path = sample_path / 'raw' / 'nucleus.tiff'
        seg_granule_path = sample_path / 'segmentations' / 'granule.tiff'
        seg_nucleus_path = sample_path / 'segmentations' / 'nucleus.tiff'
        seg_cell_path = sample_path / 'segmentations' / 'cell.tiff'
        img_rna_actin = tiff.imread(img_rna_actin_path)
        img_granule = tiff.imread(img_granule_path)
        img_rna_nanos = tiff.imread(img_rna_nanos_path)
        img_seg_granule = tiff.imread(seg_granule_path)
        img_seg_nucleus = tiff.imread(seg_nucleus_path)
        img_seg_cell = tiff.imread(seg_cell_path)
        mask_cell, mask_nucleus, mask_granule, mask_cytoplasm = compute_masks(img_seg_cell, img_seg_nucleus, img_seg_granule, safety_margin=2)
        areas_to_measure = {'cytoplasm': mask_cytoplasm, 'nucleus': mask_nucleus, 'granules': mask_granule, 'cell': mask_cell}
        images_to_measure = {'bactin_rna': img_rna_actin, 'nanos_rna': img_rna_nanos, 'granules': img_granule}
        for area_name, area_mask in areas_to_measure.items():
            for image_name, image in images_to_measure.items():
                mean_intensity = np.mean(image[area_mask])
                total_intensity = np.sum(image[area_mask])
                area_size = np.sum(area_mask)
                measurement_dict = sample_dict.copy()
                measurement_dict.update({'area_name': area_name, 'marker_name': image_name, 'mean_intensity': mean_intensity, 'total_intensity': total_intensity, 'area_size': area_size, 'pixel_size': pixel_size})
                measurement_dicts.append(measurement_dict)  # img_nucleus = tiff.imread(img_nucleus_path)
    return (measurement_dicts,)


@app.cell
def _(measurement_dicts, pd):
    df = pd.DataFrame(measurement_dicts)
    return (df,)


@app.cell
def _(df):
    df.head(5)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Show differnt pixel sizes across repeats to check for consistency
    """)
    return


@app.cell
def _(df, plt, sns):
    line_thickness = 1.5

    # sns.set_style("ticks"), rc={
    #     "axes.linewidth": line_thickness,
    #     "axes.labelsize": 15,
    #     "legend.fontsize": 12,
    #     "legend.title_fontsize": 13,
    #     "xtick.labelsize": 15,
    #     "ytick.labelsize": 15,
    #     "xtick.major.width": line_thickness,
    #     "ytick.major.width": line_thickness,
    # })

    sns.set_style("ticks")
    # plt.rcParams["axes.linewidth"] = line_thickness
    # plt.rcParams["axes.labelsize"] = 15
    # plt.rcParams["legend.fontsize"] = 12
    # plt.rcParams["legend.title_fontsize"] = 13
    # plt.rcParams["xtick.labelsize"] = 15
    # plt.rcParams["ytick.labelsize"] = 15
    # plt.rcParams["xtick.major.width"] = line_thickness
    # plt.rcParams["ytick.major.width"] = line_thickness
    # plt.rcParams["axes.labelweight"] = "bold"


    df_pixel_size = df.groupby(["uid"])[['pixel_size', "repeat_date"]].first().reset_index()
    fig_px_size, ax_px_size = plt.subplots(figsize=(4, 6))
    sns.violinplot(data=df_pixel_size, y="pixel_size", ax=ax_px_size, inner="quart", color="#124585")
    sns.swarmplot(data=df_pixel_size, y="pixel_size", hue="repeat_date", dodge=False, size=3, palette="grey")
    sns.despine()
    return


@app.cell
def _(df):
    df_mean_int = df.pivot_table(
        index=['uid', 'condition', 'sample_number', 'is_heatshock', 'marker_name'],
        columns='area_name',
        values='mean_intensity'
    ).reset_index()
    df_mean_int['cytoplasm_to_granule_ratio'] = (
        df_mean_int['cytoplasm'] / df_mean_int['granules']
    )
    df_mean_int["cytoplasm_to_total_ratio"] = (
        df_mean_int['cytoplasm'] / (df_mean_int['cytoplasm'] + df_mean_int['granules'])
    )

    df_mean_int = df_mean_int.query("marker_name != 'granules'")
    return (df_mean_int,)


@app.cell
def _(df_mean_int):
    df_mean_int.head(5)
    return


@app.cell
def _(
    df_mean_int,
    get_stat_stars,
    heat_shock_palette,
    np,
    plt,
    root,
    sns,
    statistical_analysis,
):
    fig_mean_cyto_gran, ax_mean_cyto_gran = plt.subplots(figsize=(3, 4))
    sns.violinplot(data=df_mean_int, x='marker_name', y='cytoplasm_to_granule_ratio', hue='condition', ax=ax_mean_cyto_gran, palette=heat_shock_palette, inner='quart', split=True, density_norm='width')
    sns.swarmplot(data=df_mean_int, x='marker_name', y='cytoplasm_to_granule_ratio', hue='condition', ax=ax_mean_cyto_gran, dodge=True, size=2, palette='dark:k')
    sns.despine()
    ax_mean_cyto_gran.set_title('Cytoplasm to Granule Mean Intensity Ratio')
    ax_mean_cyto_gran.set_ylabel('$\\frac{\\overline{I}_{Cytoplasm}}{\\overline{I}_{Granule}}$', fontsize=15)
    ax_mean_cyto_gran.set_xlabel('Marker')
    ax_mean_cyto_gran.legend(title='Condition')
    plt.tight_layout()
    plt.savefig(root.parent / 'cytoplasm_to_granule_ratio_violinplot.pdf', bbox_inches='tight')
    stat_results = {}
    for _marker in df_mean_int['marker_name'].unique():
        _heatshock_values = df_mean_int[(df_mean_int['marker_name'] == _marker) & df_mean_int['is_heatshock']]['cytoplasm_to_granule_ratio']
        _control_values = df_mean_int[(df_mean_int['marker_name'] == _marker) & ~df_mean_int['is_heatshock']]['cytoplasm_to_granule_ratio']
        stat_rna_results = {'Mean_ratio': round(_heatshock_values.mean() / _control_values.mean(), 4) if _control_values.mean() != 0 else np.inf, 'test_type': statistical_analysis(_heatshock_values, _control_values)[0], 'statistic': statistical_analysis(_heatshock_values, _control_values)[1], 'p_value': statistical_analysis(_heatshock_values, _control_values)[2], 'Sample_size_heatshock': len(_heatshock_values), 'Sample_size_control': len(_control_values), 'star_significance': get_stat_stars(statistical_analysis(_heatshock_values, _control_values)[2])}
        stat_results[_marker] = stat_rna_results
    print('Statistical Test Results (Cytoplasm to Granule Ratio):')
    for _marker, _results in stat_results.items():
        print(f'{_marker}:')
        for _key, _value in _results.items():
            print(f'    {_key}: {_value}')
    rep_min = 38
    rep_max = 62
    # Statistical tests
    representative_uids = {}
    for _condition in df_mean_int['condition'].unique():
        condition_df = df_mean_int[df_mean_int['condition'] == _condition]
        condition_representatives = None
        for _marker in ['bactin_rna', 'nanos_rna']:
            marker_df = condition_df[condition_df['marker_name'] == _marker]
            marker_min = np.percentile(marker_df['cytoplasm_to_granule_ratio'], rep_min)
            marker_max = np.percentile(marker_df['cytoplasm_to_granule_ratio'], rep_max)
            representative_df = marker_df[(marker_df['cytoplasm_to_granule_ratio'] >= marker_min) & (marker_df['cytoplasm_to_granule_ratio'] <= marker_max)]
            if condition_representatives is None:
                condition_representatives = representative_df['uid'].tolist()
            else:
                condition_representatives = list(set(condition_representatives) & set(representative_df['uid'].tolist()))
        representative_uids[_condition] = condition_representatives
    print('Representative UIDs for each condition:')
    for _condition, uids in representative_uids.items():
    # Get uids for representative images in each condition
        print(f'{_condition}: {uids}')  # Take the intersection of representative samples across markers
    return


@app.cell
def _(df_mean_int, heat_shock_palette, plt, root, sns):
    plt.rcParams.update(
        {
            "font.size": 6,          # default text
            "axes.titlesize": 8,     # subplot titles
            "axes.labelsize": 6,     # x/y axis labels
            "xtick.labelsize": 6,    # x tick labels
            "ytick.labelsize": 6,    # y tick labels
            "legend.fontsize": 6,    # legend
            "figure.titlesize": 8,   # suptitle
        }
    )

    import matplotlib as mpl

    mpl.rcParams["font.family"] = "Arial"

    fig_mean_cyto_gran_newstyle, ax_mean_cyto_gran_newstyle = plt.subplots(1, 1, figsize=(1.63, 2.5))

    sns.stripplot(
        data=df_mean_int,
        x='marker_name',
        y='cytoplasm_to_granule_ratio',
        hue='condition',
        dodge=True,
        palette=heat_shock_palette,
        alpha=0.4,
        size=3,
        ax=ax_mean_cyto_gran_newstyle,
        jitter=0.3,
    )
    sns.pointplot(
        data=df_mean_int,
        x='marker_name',
        y='cytoplasm_to_granule_ratio',
        hue='condition',
        dodge=0.4,
        errorbar="sd",
        estimator="median",
        capsize=0.075,
        linestyle="none",
        markersize=10,
        marker="_",
        err_kws=dict(linewidth=0.4, color="black"),
        markeredgewidth=1,
        palette="dark:black",
        zorder=5,
        ax=ax_mean_cyto_gran_newstyle,
    )

    # Fix legend (pointplot + stripplot will both add hue entries)
    handles, labels = ax_mean_cyto_gran_newstyle.get_legend_handles_labels()
    n_cond = df_mean_int['condition'].nunique()
    ax_mean_cyto_gran_newstyle.legend(handles[:n_cond], labels[:n_cond], title='Condition')

    ax_mean_cyto_gran_newstyle.set_title('Cytoplasm to Granule Mean Intensity Ratio')
    ax_mean_cyto_gran_newstyle.set_ylabel('$\\frac{\\bar{I}_{Cytoplasm}}{\\bar{I}_{Granule}}$', fontsize=6)
    ax_mean_cyto_gran_newstyle.set_xlabel('Marker')

    plt.tight_layout()
    sns.despine()
    fig_mean_cyto_gran_newstyle.savefig(root.parent / 'cytoplasm_to_granule_ratio_pointplot.pdf', bbox_inches='tight')
    return


@app.cell
def _(
    df_mean_int,
    get_stat_stars,
    heat_shock_palette,
    np,
    plt,
    root,
    sns,
    statistical_analysis,
):
    fig_mean_cyto_total, ax_mean_cyto_total = plt.subplots(figsize=(4, 6))
    sns.violinplot(data=df_mean_int, x='marker_name', y='cytoplasm_to_total_ratio', hue='condition', ax=ax_mean_cyto_total, palette=heat_shock_palette, inner='quart', split=True, density_norm='width')
    sns.swarmplot(data=df_mean_int, x='marker_name', y='cytoplasm_to_total_ratio', hue='condition', ax=ax_mean_cyto_total, dodge=True, size=3, palette='dark:k')
    sns.despine()
    ax_mean_cyto_total.set_title('Cytoplasm to Total Mean Intensity Ratio')
    ax_mean_cyto_total.set_ylabel('$\\frac{\\overline{I}_{Cytoplasm}}{\\overline{I}_{Granule} + \\overline{I}_{Cytoplasm}}$', fontsize=15)
    ax_mean_cyto_total.set_xlabel('Marker')
    ax_mean_cyto_total.legend(title='Condition')
    plt.savefig(root.parent / 'cytoplasm_to_total_ratio_violinplot.pdf', bbox_inches='tight')
    stat_results_total = {}
    for _marker in df_mean_int['marker_name'].unique():
        _heatshock_values = df_mean_int[(df_mean_int['marker_name'] == _marker) & df_mean_int['is_heatshock']]['cytoplasm_to_total_ratio']
        _control_values = df_mean_int[(df_mean_int['marker_name'] == _marker) & ~df_mean_int['is_heatshock']]['cytoplasm_to_total_ratio']
        stat_rna_results_total = {'Mean_ratio': round(_heatshock_values.mean() / _control_values.mean(), 4) if _control_values.mean() != 0 else np.inf, 'test_type': statistical_analysis(_heatshock_values, _control_values)[0], 'statistic': statistical_analysis(_heatshock_values, _control_values)[1], 'Sample_size_heatshock': len(_heatshock_values), 'Sample_size_control': len(_control_values), 'p_value': statistical_analysis(_heatshock_values, _control_values)[2], 'star_significance': get_stat_stars(statistical_analysis(_heatshock_values, _control_values)[2])}
        stat_results_total[_marker] = stat_rna_results_total
    print('Statistical Test Results (Cytoplasm to Total Ratio):')
    for _marker, _results in stat_results_total.items():
        print(f'{_marker}:')
        for _key, _value in _results.items():
    # Get statistical test results for cytoplasm to total ratio
            print(f'        {_key}: {_value}')
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
