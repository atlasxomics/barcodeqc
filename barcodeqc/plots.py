import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings

from matplotlib.ticker import FuncFormatter
from pathlib import Path

logger = logging.getLogger(__name__)


def create_density_plot(
    dataframe: pd.DataFrame,
    outPath: str,
    data_column: str,
    group_column: str,
    log10: bool = True,
    x_label: str = '',
    y_label: str = '',
    plotTitle: str = ''
):
    """Generate and save a density plot; used downstream to plot the
    distribution of counts per barcode (_denseOnOff.png').
    """
    plt.figure(figsize=(6, 6))

    # Use the original palette for up to two groups
    palette_full = sns.color_palette("bright", 10)
    palette = [palette_full[0], palette_full[3]]

    # Get unique, non-null groups
    groups = [g for g in dataframe[group_column].dropna().unique()
              if len(dataframe[dataframe[group_column] == g]) > 0]

    if len(groups) == 0:
        logger.debug("No non-empty groups to plot in density plot.")
        return

    # Plot each non-empty group with the correct color
    for idx, group in enumerate(groups):
        group_df = dataframe[dataframe[group_column] == group]
        sns.kdeplot(
            data=group_df,
            x=data_column,
            bw_adjust=1.5,
            log_scale=log10,
            fill=True,
            color=palette[idx],
            label=f"{group}"
        )

    if log10:
        x_label = x_label + " (log scaled)"

    plt.title(plotTitle, fontsize=16, fontweight='bold')
    plt.xlabel(x_label, fontsize=14, fontweight='bold')
    plt.ylabel(y_label, fontsize=14, fontweight='bold')
    plt.gca().yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{y:.2f}'))
    plt.legend(loc='upper right', frameon=False)
    if log10:
        plt.xscale('symlog', linthresh=1)
        plt.xlim(0, 10000)
    else:
        plt.xlim(0, 10000)
    plt.savefig(f'{outPath}', dpi=300)
    plt.close()


def create_heatmap(
    df,
    outPath,
    log10=True,
    axesOff=False,
    countCol='nCellHits',
    vmin=None,
    vmax=None,
    colorMap='RdYlBu_r'
):
    '''Create spatial heatmap arranged by row/column, colored by read counts'''

    # Pivot the dataframe to create a matrix
    matrix = df.pivot(index='row', columns='col', values=countCol)
    cbLbl = "counts"

    # Take the log10 of the values in the matrix
    with warnings.catch_warnings():
        # if there is a zero, log10 throws a warning.
        # ignoring that warning
        warnings.simplefilter("ignore", RuntimeWarning)
        if log10:
            matrix = np.log10(matrix)
            cbLbl = "log10 of counts"
            if vmin is not None:
                vmin = np.log10(vmin)
            if vmax is not None:
                vmax = np.log10(vmax)

    # Create the heatmap
    fig, ax = plt.subplots(figsize=(15, 12))

    im = ax.imshow(matrix, vmin=vmin, vmax=vmax, cmap=colorMap)

    # Add labels
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=6)
    ax.set_yticklabels(matrix.index, fontsize=6)

    if axesOff:
        plt.axis('off')

    # Add colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    # cbar.set_label(cbLbl, rotation=270)  # changed to add labelpad
    cbar.set_label(cbLbl, rotation=270, labelpad=15)

    plt.tight_layout()  # Automatically adjust subplot parameters

    # save the plot
    # plt.savefig(outPath)  # change to reduce whitespace
    plt.savefig(outPath, bbox_inches='tight', pad_inches=0.1)

    logger.debug(f"create_heatmap: Saving to {outPath}")
    plt.close()

    return matrix


def hilo_plot(
    df: pd.DataFrame,
    xval: str,
    yval: str,
    label_col: str,
    output_dir: Path,
    file_name: str,
) -> Path:
    """
    Create a bar plot with mean, 0.5x mean, and 2x mean reference lines.

    Parameters
    ----------
    df : pd.DataFrame
        Input data.
    xval : str
        Column used for bar positions.
    yval : str
        Numeric column to plot.
    label_col : str
        Column used for x-axis tick labels.
    output_dir : Path
        Directory to write plot.
    file_name : str
        Output filename.

    Returns
    -------
    Path
        Path to saved plot.
    """

    # ---- Validate y values ----
    y = pd.to_numeric(df[yval], errors="coerce")
    if y.isna().all():
        raise ValueError(f"Column '{yval}' must contain numeric values")

    mean = y.mean()
    lower = mean / 2
    upper = mean * 2

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(20, 5))

    df.plot(
        kind="bar",
        x=xval,
        y=yval,
        ax=ax,
        legend=False,
    )

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df[label_col], rotation=90)

    ax.axhline(
        mean,
        color="tab:gray",
        linestyle="--",
        linewidth=0.8,
        label=f"mean: {mean:.3f}",
    )
    ax.axhline(
        lower,
        color="tab:red",
        linestyle="--",
        linewidth=0.8,
        label=f"0.5× mean: {lower:.3f}",
    )
    ax.axhline(
        upper,
        color="tab:red",
        linestyle="--",
        linewidth=0.8,
        label=f"2× mean: {upper:.3f}",
    )

    ax.legend(loc="upper right")
    fig.tight_layout()

    # ---- Save ----
    output_dir.mkdir(parents=True, exist_ok=True)
    outpath = output_dir / file_name
    fig.savefig(outpath)
    plt.close(fig)

    return outpath


def pareto_plot(
    df: pd.DataFrame,
    y_col1: str,
    y_col2: str,
    colorby: str,
    channel_label: str,
    wildcard_path: Path,
    max_to_ninety: int,
    output_dir: Path,
    file_name: str,
):

    startI = 0
    endI = max_to_ninety + int(max_to_ninety / 2)

    # Ensure endI doesn't exceed dataframe length
    endI = min(endI, len(df))

    # Filter to only valid (non-NaN, finite) values
    slice_data = df[startI:endI].copy()
    slice_data = slice_data[
        (slice_data[y_col1].notna()) & 
        (slice_data[y_col1] > 0) &
        (np.isfinite(slice_data[y_col1]))
    ]

    if len(slice_data) == 0:
        logger.warning(f"No valid data for pareto plot: {file_name}")
        # Create empty plot
        fig, ax = plt.subplots(figsize=(20, 10))
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center')
        outpath = output_dir / file_name
        fig.savefig(outpath)
        plt.close(fig)
        return outpath
    
    fig, ax = plt.subplots(figsize=(20, 10))

    ax.plot(
        slice_data[y_col1],
        marker="x",
        linewidth=1,
        label="Fraction Total"
    )

    # define second y-axis that shares x-axis with current plot
    ax2 = ax.twinx()

    # add second line to plot
    ax2.plot(
        slice_data[y_col2],
        marker='.',
        linewidth=1,
        color='tab:orange',
        label='Cumulative Fraction'
    )

    # add second y-axis label
    labels0 = list(slice_data.index)

    chanLabels = list(
        "_" + (slice_data[channel_label] + 1).astype(str).str.zfill(2)
    )
    labels3 = [f"{lb}{cb}" for lb, cb in zip(labels0, chanLabels)]

    ax.set_xticks(range(len(labels0)))
    ax.set_xticklabels(labels3, rotation=90)

    s = slice_data[colorby]
    makeGreen = np.where(s)[0].tolist()
    for mg in makeGreen:
        if mg < len(ax.get_xticklabels()):
            ax.get_xticklabels()[mg].set_color("green")

    # Only add axhline if ax2 has valid data limits
    try:
        ax2.axhline(
            y=.9,
            color='tab:gray',
            linestyle='--',
            linewidth=0.5,
            label='90 percent'
        )
    except np.linalg.LinAlgError:
        # Skip if axis transform is singular (e.g., empty data)
        pass
    
    fig.legend(loc="center right")
    plt.title(f"plot: {wildcard_path}")
    fig.subplots_adjust(bottom=0.2)

    pareto_path = output_dir / file_name
    fig.savefig(pareto_path)
    plt.close(fig)

    return pareto_path
