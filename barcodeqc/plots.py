import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.stats import gaussian_kde

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
    groups = [
        g for g in dataframe[group_column].dropna().unique()
        if len(dataframe[dataframe[group_column] == g]) > 0
    ]

    if len(groups) == 0:
        logger.debug("No non-empty groups to plot in density plot.")
        return

    palette = ["#1f77b4", "#2ca02c"]
    fig = go.Figure()

    for idx, group in enumerate(groups):
        group_df = dataframe[dataframe[group_column] == group]
        values = pd.to_numeric(group_df[data_column], errors="coerce").dropna()
        values = values[values > 0]
        if values.empty:
            continue

        if log10:
            log_vals = np.log10(values)
            kde = gaussian_kde(log_vals)
            xs = np.linspace(
                max(log_vals.min(), 0),
                log_vals.max(),
                300,
            )
            ys = kde(xs)
            xs_plot = np.power(10, xs)
        else:
            kde = gaussian_kde(values)
            xs_plot = np.linspace(0, min(values.max(), 10000), 300)
            ys = kde(xs_plot)

        fig.add_trace(
            go.Scatter(
                x=xs_plot,
                y=ys,
                mode="lines",
                line=dict(color=palette[idx % len(palette)], width=2),
                fill="tozeroy",
                name=str(group),
            )
        )

    if log10:
        x_label = x_label + " (log scaled)"

    fig.update_layout(
        title=plotTitle,
        xaxis_title=x_label,
        yaxis_title=y_label,
        template="plotly_white",
        height=360,
        margin=dict(l=40, r=20, t=50, b=40),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.92,
            xanchor="right",
            x=0.96,
            font=dict(size=14),
        ),
        autosize=True,
    )
    fig.update_yaxes(tickformat=".2f")
    if log10:
        fig.update_xaxes(type="log", range=[0, 4])
    else:
        fig.update_xaxes(range=[0, 10000])

    outpath = Path(outPath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(outpath),
        include_plotlyjs=False,
        full_html=False,
        config={"responsive": True, "displaylogo": False},
    )


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

    colorscale = colorMap
    reversescale = False
    if isinstance(colorMap, str) and colorMap.endswith("_r"):
        colorscale = colorMap[:-2]
        reversescale = True

    fig = go.Figure(
        data=go.Heatmap(
            z=matrix.values,
            x=matrix.columns.astype(str),
            y=matrix.index.astype(str),
            colorscale=colorscale,
            reversescale=reversescale,
            zmin=vmin,
            zmax=vmax,
            colorbar=dict(title=cbLbl),
        )
    )
    fig.update_layout(
        template="plotly_white",
        height=500,
        margin=dict(l=40, r=20, t=40, b=40),
        autosize=True,
    )
    if axesOff:
        fig.update_xaxes(visible=False)
        fig.update_yaxes(visible=False)

    outpath = Path(outPath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    pio.write_html(
        fig,
        file=str(outpath),
        include_plotlyjs=False,
        full_html=False,
        config={"responsive": True, "displaylogo": False},
    )

    logger.debug(f"create_heatmap: Saving to {outPath}")

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

    labels = df[label_col].astype(str).tolist()
    x_vals = list(range(len(df)))

    channel_raw = pd.to_numeric(df[xval], errors="coerce")
    channel_vals = channel_raw.map(
        lambda v: str(int(v)) if pd.notna(v) else "NA"
    )
    customdata = np.column_stack([labels, channel_vals])

    fig = go.Figure(
        data=go.Bar(
            x=x_vals,
            y=y,
            marker_color="#225789",
            customdata=customdata,
            hovertemplate=(
                "Barcode: %{customdata[0]}<br>"
                "Channel: %{customdata[1]}<br>"
                "%{y:.4f}<extra></extra>"
            ),
        )
    )
    if "L1" in file_name:
        title_text = "Barcode A: Read Fraction, sorted by channel"
    elif "L2" in file_name:
        title_text = "Barcode B: Read Fraction, sorted by channel"
    else:
        title_text = None

    n_labels = len(labels)
    if n_labels > 0:
        step = max(1, n_labels // 40)
        tickvals = [i for i in range(n_labels) if i % step == 0]
        ticktext = [labels[i] for i in tickvals]
        fig.update_xaxes(
            tickmode="array", tickvals=tickvals, ticktext=ticktext
        )

    fig.add_hline(
        y=mean,
        line_dash="dash",
        line_color="gray",
        annotation_text=f"mean: {mean:.3f}",
        annotation_position="top right",
        annotation_bgcolor="rgba(255,255,255,0.75)",
    )
    fig.add_hline(
        y=lower, line_dash="dash",
        line_color="red",
        annotation_text=f"0.5× mean: {lower:.3f}",
        annotation_position="top right",
        annotation_bgcolor="rgba(255,255,255,0.75)",
    )
    fig.add_hline(
        y=upper,
        line_dash="dash",
        line_color="red",
        annotation_text=f"2× mean: {upper:.3f}",
        annotation_position="top right",
        annotation_bgcolor="rgba(255,255,255,0.75)",
    )

    fig.update_layout(
        template="plotly_white",
        height=320,
        margin=dict(l=60, r=20, t=50, b=60),
        autosize=True,
    )
    if title_text:
        fig.update_layout(title=dict(text=title_text, x=0.5, xanchor="center"))
    fig.update_yaxes(title_text="Fraction of Reads")

    # ---- Save ----
    output_dir.mkdir(parents=True, exist_ok=True)
    outpath = output_dir / file_name
    pio.write_html(
        fig,
        file=str(outpath),
        include_plotlyjs=False,
        full_html=False,
        config={"responsive": True, "displaylogo": False},
    )

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
        fig = go.Figure()
        fig.add_annotation(
            text="No valid data",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
        )
        fig.update_layout(
            template="plotly_white",
            height=420,
            margin=dict(l=40, r=40, t=40, b=40),
        )
        outpath = output_dir / file_name
        pio.write_html(
            fig,
            file=str(outpath),
            include_plotlyjs=False,
            full_html=False,
            config={"responsive": True, "displaylogo": False},
        )
        return outpath

    x_vals = list(range(len(slice_data)))

    labels0 = list(slice_data.index)
    chan_labels = list(
        "_" + (slice_data[channel_label] + 1).astype(str).str.zfill(2)
    )
    labels3 = [f"{lb}{cb}" for lb, cb in zip(labels0, chan_labels)]

    s = slice_data[colorby].astype(bool)
    marker_colors = np.where(s, "green", "#1f2937")
    channel_raw = pd.to_numeric(slice_data[channel_label], errors="coerce")
    channel_vals = channel_raw.map(
        lambda v: str(int(v + 1)) if pd.notna(v) else "NA"
    )
    status_vals = np.where(s, "Expected", "Unexpected")
    customdata = np.column_stack([labels0, channel_vals, status_vals])

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=slice_data[y_col1],
            mode="lines+markers",
            marker=dict(color=marker_colors, size=4),
            line=dict(color="#4c78a8", width=1),
            name="Fraction Total",
            customdata=customdata,
            hovertemplate=(
                "Barcode: %{customdata[0]}<br>"
                "Channel: %{customdata[1]}<br>"
                "Status: %{customdata[2]}<br>"
                "%{y:.2%}<extra></extra>"
            ),
        ),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(
            x=x_vals,
            y=slice_data[y_col2],
            mode="lines+markers",
            marker=dict(color="gray", size=4),
            line=dict(color="gray", width=1),
            name="Cumulative Fraction",
            customdata=customdata,
            hovertemplate=(
                "Barcode: %{customdata[0]}<br>"
                "Channel: %{customdata[1]}<br>"
                "Status: %{customdata[2]}<br>"
                "%{y:.2%}<extra></extra>"
            ),
        ),
        secondary_y=True,
    )

    n_labels = len(labels3)
    if n_labels > 0:
        step = max(1, n_labels // 40)
        tickvals = [i for i in range(n_labels) if i % step == 0]
        ticktext = [labels3[i] for i in tickvals]
        fig.update_xaxes(tickmode="array", tickvals=tickvals, ticktext=ticktext)

    if "L1" in wildcard_path.name:
        title_text = "Barcode A 8mers by fraction of total"
    elif "L2" in wildcard_path.name:
        title_text = "Barcode B 8mers by fraction of total"
    else:
        title_text = f"plot: {wildcard_path}"

    fig.update_layout(
        title=dict(text=title_text, x=0.5, xanchor="center"),
        template="plotly_white",
        height=420,
        margin=dict(l=40, r=40, t=40, b=80),
        legend=dict(
            orientation="v",
            yanchor="top",
            y=0.75,
            xanchor="right",
            x=0.9,
            font=dict(size=12),
        ),
        autosize=True,
    )
    y1_max = float(slice_data[y_col1].max()) if len(slice_data) else 0.0
    y1_max = y1_max * 1.1 if y1_max > 0 else 0.05
    fig.update_yaxes(
        title_text="Fraction Total",
        secondary_y=False,
        range=[0, y1_max],
        tickformat=".2%",
        showgrid=False
    )
    fig.update_yaxes(
        title_text="Cumulative Fraction",
        secondary_y=True,
        range=[0, 1],
        tickformat=".0%",
        showgrid=False
    )

    pareto_path = output_dir / file_name
    pio.write_html(
        fig,
        file=str(pareto_path),
        include_plotlyjs=False,
        full_html=False,
        config={"responsive": True, "displaylogo": False},
    )

    return pareto_path
