import logging
from typing import Any, List

import numpy as np
import pandas as pd
import RNA  # Requires the ViennaRNA python package

import bokeh
import bokeh.palettes
import bokeh.transform
import bokeh.models
from bokeh.plotting import figure

# =============================================================================
# SPATIAL & DATA TRANSFORMATION UTILITIES
# =============================================================================

def mfe_to_ppair_table(mfe: str) -> pd.DataFrame:
    """
    Generate a table of nucleic acid base pairs from a secondary structure 
    in dot-bracket notation.
    
    Args:
        mfe (str): Nucleic acid secondary structure in dot-bracket notation 
                   (e.g., '..((..))..').
                   
    Returns:
        pd.DataFrame: A dataframe with columns 'pos1' and 'pos2' representing 
                      the 0-based indices of paired nucleotides. Unpaired 
                      nucleotides are mapped to themselves (pos1 == pos2).
    """
    res = RNA.ptable(mfe, options=RNA.BRACKETS_DEFAULT)
    L = res.get(0)
    
    # RNA.ptable returns 1-based indices, and 0 for unpaired.
    # We convert to 0-based indices. Unpaired becomes -1.
    out = [[i, res.get(i + 1) - 1] for i in range(L)]
    df = pd.DataFrame(out, columns=["pos1", "pos2"])
    
    # Map unpaired nucleotides to themselves to avoid -1 index issues
    unpaired_mask = df["pos2"] == -1
    df.loc[unpaired_mask, "pos2"] = df.loc[unpaired_mask, "pos1"]
    
    return df


def get_mfe_xy(df: pd.DataFrame, mfe: str) -> pd.DataFrame:
    """
    Calculate x,y positions needed for drawing nucleic acid secondary structures.
    
    Args:
        df (pd.DataFrame): DataFrame containing at least 'pos1' and 'pos2' columns.
        mfe (str): Minimum free energy structure in dot-bracket notation.
                   
    Returns:
        pd.DataFrame: The original dataframe merged with calculated 'x' and 'y' coordinates.
    """
    pp_table = mfe_to_ppair_table(mfe)
    df_merged = pp_table.merge(df, on=["pos1", "pos2"], how="left")
    
    coord = RNA.get_xy_coordinates(mfe)
    df_merged["x"] = [coord.get(i).X for i in range(len(df_merged))]
    df_merged["y"] = [coord.get(i).Y for i in range(len(df_merged))]
    
    return df_merged


def rotate_xy(df: pd.DataFrame, xcol: str, ycol: str, rad: float) -> pd.DataFrame:
    """
    Rotate the coordinates of the nucleic acid structure around the origin.
    
    Args:
        df (pd.DataFrame): DataFrame containing the coordinates.
        xcol (str): Name of the column containing X coordinates.
        ycol (str): Name of the column containing Y coordinates.
        rad (float): Rotation angle in radians.
        
    Returns:
        pd.DataFrame: A new dataframe with rotated X and Y coordinates.
    """
    df_plot = df.copy()
    cos_rad, sin_rad = np.cos(rad), np.sin(rad)
    
    x, y = df_plot[xcol].values, df_plot[ycol].values
    
    # Vectorized 2D rotation
    df_plot[xcol] = x * cos_rad - y * sin_rad
    df_plot[ycol] = x * sin_rad + y * cos_rad
    
    return df_plot


# =============================================================================
# PLOTTING UTILITIES
# =============================================================================

def draw_pair_probability(
    p: figure, 
    df: pd.DataFrame, 
    nuc_size: int = 5, 
    palette: List[str] = bokeh.palettes.Turbo256
) -> bokeh.models.GlyphRenderer:
    """
    Draw a pair probability matrix as a scatter plot.
    
    Args:
        p (figure): The Bokeh figure object to draw on.
        df (pd.DataFrame): DataFrame containing 'pos1', 'pos2', and 'pair_probability'.
        nuc_size (int): Size of the scatter points.
        palette (List[str]): Color palette for the probability mapping.
        
    Returns:
        bokeh.models.GlyphRenderer: The renderer object for the scatter plot.
    """
    pp_color = bokeh.transform.linear_cmap(field_name="pair_probability", palette=palette, low=0, high=1)        
    img = p.circle(
        source=df, x="pos1", y="pos2", size=nuc_size, 
        line_color=None, fill_color=pp_color, line_width=1
    )
    
    tooltips = [("position", "@pos1 to @pos2"), ("probability", "@pair_probability")]
    hvr = bokeh.models.HoverTool(tooltips=tooltips, renderers=[img])
    p.add_tools(hvr)

    color_mapper = bokeh.models.LinearColorMapper(palette=palette, low=0, high=1)
    color_bar = bokeh.models.ColorBar(color_mapper=color_mapper, location=(0, 0), title="pair probability")
    p.add_layout(color_bar, "right")

    p.xaxis.axis_label = "position1 (nt)"
    p.yaxis.axis_label = "position2 (nt)"
    
    return img


def draw_secondary_structure(
    p: figure, 
    df: pd.DataFrame, 
    line_width: int = 2, 
    nuc_color: str = "black", 
    nuc_size: int = 5, 
    nuc_letters: bool = False
) -> figure:
    """
    Draw a nucleic acid secondary structure with backbone, base pairs, and nucleotides.
    
    Args:
        p (figure): The Bokeh figure object to draw on.
        df (pd.DataFrame): DataFrame containing 'x', 'y', 'position', 'pair_position', 
                           and optionally 'nucleotide' and other color mapping columns.
        line_width (int): Width of the backbone and base pair lines.
        nuc_color (str): Column name for nucleotide coloring, or a literal color string.
        nuc_size (int): Size of the nucleotide circles.
        nuc_letters (bool): If True, draws the nucleotide letters (A, U, C, G) on the plot.
        
    Returns:
        figure: The modified Bokeh figure object.
    """
    # 1. Draw backbone helices
    df_bb = df[["x", "y"]].copy()
    p.line(source=df_bb, x="x", y="y", line_color="black", line_width=line_width)
    
    # Add end arrowhead tick to indicate 5' -> 3' direction
    oh = bokeh.models.OpenHead(size=20, fill_color="black")
    x_vals, y_vals = df_bb["x"].values, df_bb["y"].values
    
    if len(x_vals) >= 2:
        x_start, y_start = x_vals[-2], y_vals[-2]
        x_end, y_end = x_vals[-1], y_vals[-1]
        x_arrow = x_end + (x_end - x_start)
        y_arrow = y_end + (y_end - y_start)
        p.add_layout(bokeh.models.Arrow(
            end=oh, line_color="black", line_width=2, 
            x_start=x_end, y_start=y_end, x_end=x_arrow, y_end=y_arrow
        ))
    
    # 2. Draw base pairs
    df_plot = df.copy()
    
    # Optimized lookup for base pair coordinates
    x_lookup = dict(zip(df_plot["position"], df_plot["x"]))
    y_lookup = dict(zip(df_plot["position"], df_plot["y"]))
    
    df_plot["xs"] = [[x, x_lookup.get(pair)] for x, pair in zip(df_plot["x"], df_plot["pair_position"])]
    df_plot["ys"] = [[y, y_lookup.get(pair)] for y, pair in zip(df_plot["y"], df_plot["pair_position"])]
    
    # Filter to draw each base pair only once
    df_bp = df_plot[df_plot["position"] < df_plot["pair_position"]].copy()
    
    palette = bokeh.palettes.Turbo256
    bp_cmap = bokeh.transform.linear_cmap(field_name="pair_probability", palette=palette, low=0, high=1)
    
    bpairs = p.multi_line(
        source=df_bp, xs="xs", ys="ys", 
        line_color=bp_cmap, line_width=line_width
    )

    bp_tooltips = [("pair_position", "@position to @pair_position"), ("probability", "@pair_probability")]
    p.add_tools(bokeh.models.HoverTool(tooltips=bp_tooltips, renderers=[bpairs]))
    
    cbar_cmap = bokeh.models.LinearColorMapper(palette=palette, low=0, high=1)
    cbar = bokeh.models.ColorBar(color_mapper=cbar_cmap, location=(0, 0), title="pair_probability")
    p.add_layout(cbar, "right")
    
    # 3. Draw nucleotides
    # --- INLINED LINEAR CMAP LOGIC ---
    # Maps a dataframe column to a color gradient, or falls back to a literal color string.
    if nuc_color in df_plot.columns:
        y = df_plot[nuc_color]
        low = 0.0 if "probability" in nuc_color.lower() else float(y.min())
        high = 1.0 if "probability" in nuc_color.lower() else float(y.max())
        nuc_color_mapped = bokeh.transform.linear_cmap(
            field_name=nuc_color, 
            palette=bokeh.palettes.Turbo256, 
            low=low, 
            high=high
        )
    else:
        nuc_color_mapped = nuc_color
    # ---------------------------------

    nuc = p.circle(
        source=df_plot, size=nuc_size, 
        line_color="black", fill_color=nuc_color_mapped, line_width=2
    )
    
    nuc_tooltips = [
        ("position", "@position"),
        ("probability", "@pair_probability"),
        ("nucleotide", "@nucleotide")
    ]
    if isinstance(nuc_color, str) and nuc_color in df_plot.columns:
        nuc_tooltips.append((nuc_color, f"@{nuc_color}"))
        
    p.add_tools(bokeh.models.HoverTool(tooltips=nuc_tooltips, renderers=[nuc]))

    # 4. Optionally draw nucleotide letters
    if nuc_letters and "nucleotide" in df_plot.columns:
        p.text(
            source=df_plot, x="x", y="y", 
            x_offset=0, y_offset=10, text_align="center", 
            text_color="black", text="nucleotide"
        )
        
    return p