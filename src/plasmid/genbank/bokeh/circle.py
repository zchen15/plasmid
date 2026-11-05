# science libs
import numpy as np
import pandas as pd

from ..CirclePlots import CircleLayout

# plotting libraries
import bokeh.plotting
import bokeh.models

class Plots(CircleLayout):
    """Visualizes GenBank annotation information in a circular format.

    This class provides methods to render circular maps of plasmids, phages,
    or circular genomes using the Bokeh plotting library. Layout and coordinate
    calculations are provided by :class:`CirclePlots.CircleLayout`.
    """
    
    def autoscale_range(self):
        """Auto scale plot ranges for a circular plot.

        Adjusts the x and y ranges of the Bokeh figure based on the plot's
        aspect ratio and the circle's radius to ensure the plot is properly
        centered and scaled.

        Returns:
            bokeh.plotting.Figure: The updated Bokeh plotting figure object.
        """
        hw_ratio = self.p.height / self.p.width
        R = self.radius0
        Rmin = 1.5 * R
        w = R * 2
        
        h = w * hw_ratio
        if hw_ratio < 1 and h < Rmin:
            h = Rmin
            w = h / hw_ratio
        
        x0, y0 = self.xy0
        self.p.x_range = bokeh.models.Range1d(-w + x0, w + x0)
        self.p.y_range = bokeh.models.Range1d(-h + y0, h + y0)
        return self.p
        
    def add_bounds(self, color="black", alpha=1):
        """Add maximum and minimum boundary lines for the circular plot.

        Draws dashed circular lines representing the inner and outer bounds
        of the annotation track.

        Args:
            color (str, optional): The color of the boundary lines. 
                Defaults to "black".
            alpha (float, optional): The transparency of the boundary lines. 
                Defaults to 1.

        Returns:
            bokeh.plotting.Figure: The updated Bokeh plotting figure object.
        """
        seqlen = len(self.sequence)
        tval = np.linspace(0, seqlen, 200) / seqlen * 2 * np.pi
        
        # draw lower bound line
        xy = self.polar2xy(tval, self.radius2)
        df1 = pd.DataFrame(xy, columns=["xs", "ys"])
        
        # draw upper bound line
        xy = self.polar2xy(tval, self.radius3) 
        df2 = pd.DataFrame(xy, columns=["xs", "ys"])
        
        self.p.line(source=df1, x="xs", y="ys", line_color=color, line_dash="dashed", line_width=2, alpha=alpha)
        self.p.line(source=df2, x="xs", y="ys", line_color=color, line_dash="dashed", line_width=2, alpha=alpha)
        return self.p
        
    def add_features(self):
        """Draw sequence features as circular arcs.

        Renders the sequence backbone and plots features (such as genes) as
        colored arrow patches. Also attaches a HoverTool to display feature
        metadata on mouseover.

        Returns:
            bokeh.plotting.Figure: The updated Bokeh plotting figure object.
        """
        # get base dimensions
        seqlen = len(self.sequence)
        
        # set plotting origin
        x0, y0 = self.xy0
        
        # draw backbone line
        xy = self.get_arc_xy(0, seqlen, self.radius0, self.radius0, seqlen)
        xy[:, 0] = xy[:, 0] + x0
        xy[:, 1] = xy[:, 1] + y0
        self.p.line(x=xy[:, 0], y=xy[:, 1], line_color="black", line_width=2)
        
        # plot arrows as features
        df = self.get_arrow()
        arcs = self.p.patches(source=df, xs="xs", ys="ys", fill_color="color", fill_alpha="alpha",
                              line_color="black", line_width=3)
        
        # enable hovertools for arcs only
        TOOLTIPS = [
            ("name", "@name"),
            ("feature", "@feature"),
            ("color", "@color"),
            ("location", "@location")
        ]
        hvr = bokeh.models.HoverTool(tooltips=TOOLTIPS)
        hvr.renderers = [arcs]
        self.p.add_tools(hvr)
        return self.p
        
    def add_labels(self):
        """Draw text labels and leader lines for sequence features.

        Calculates label positions and draws multi-lines (ticks) connecting
        the features to their corresponding text labels. Text alignment is
        dynamically adjusted based on the label's position relative to the
        center of the plot.

        Returns:
            bokeh.plotting.Figure: The updated Bokeh plotting figure object.
        """
        # draw feature label ticks
        df = self.get_label_positions(self.features)
        self.p.multi_line(source=df, xs="Lx", ys="Ly", line_color="color", line_width=2, alpha="alpha")
        
        # draw feature label text
        c = df["text_alignment"] == "left"
        self.p.text(source=df[c], x="Tx", y="Ty", x_offset=3, y_offset=8, text_align="text_alignment", text="name")
        self.p.text(source=df[~c], x="Tx", y="Ty", x_offset=-5, y_offset=8, text_align="text_alignment", text="name")
        return self.p
        
    def add_radial_ticks(self):
        """Draw nucleotide position ticks and the total sequence length.

        Adds tick marks around the perimeter of the circular plot to indicate
        nucleotide positions. Also displays the total sequence length (e.g.,
        "5000 bp") in the center of the plot.

        Returns:
            bokeh.plotting.Figure: The updated Bokeh plotting figure object.
        """
        # Standardized to self.xy0 to match add_features
        x0, y0 = self.xy0 
        seqlen = len(self.sequence)
        
        # print size of plasmid
        val = [[x0, y0, "center", str(seqlen) + " bp"]]
        df = pd.DataFrame(val, columns=["x", "y", "text_alignment", "text"])
        self.p.text(source=df, x="x", y="y", text_align="text_alignment", text="text")
        
        # draw tick marks
        df = self.get_arc_ticks()
        self.p.multi_line(source=df, xs="Lx", ys="Ly", line_color="black", line_width=2, alpha=0.5)
        
        # draw label text
        c = df["text_alignment"] == "left"
        self.p.text(source=df[c], x="Tx", y="Ty", x_offset=-5, y_offset=8, text_align="text_alignment", text="name")
        self.p.text(source=df[~c], x="Tx", y="Ty", x_offset=5, y_offset=8, text_align="text_alignment", text="name")
        return self.p
        
    def crossplot(self, source, **kwargs):
        """Overlay crossing lines (chords) onto the circular plot.

        Typically used to visualize sequence alignments, synteny, or structural
        variations between different tracks.

        Args:
            source (pandas.DataFrame): A DataFrame containing the line
                coordinate information (must include 'xs' and 'ys' columns).
            **kwargs: Additional keyword arguments passed directly to
                Bokeh's `multi_line` glyph method (e.g., `line_color`,
                `line_width`, `alpha`).

        Returns:
            bokeh.models.GlyphRenderer: The Bokeh glyph renderer object for
                the multi-line plot.
        """        
        # Set default plot properties
        plot_kwargs = {
            "line_color": "black",
            "line_width": 2,
            "alpha": 0.1,
        }
        # Allow kwargs to override defaults
        plot_kwargs.update(kwargs)
        
        glyph = self.p.multi_line(source=source, xs="xs", ys="ys", **plot_kwargs)
        return glyph
        
    def polarplot(self, source, x, y, **kwargs):
        """Plot amplitude data radially around the origin.

        Maps linear amplitude data (such as GC content or coverage) to a
        normalized radius and plots it as a line around the circular track.

        Args:
            source (pandas.DataFrame): A DataFrame containing the data to plot.
            x (str): The column name containing nucleotide position information.
            y (str): The column name containing amplitude information to be
                converted to a normalized radius.
            **kwargs: Additional keyword arguments passed directly to
                Bokeh's `line` glyph method (e.g., `line_color`, `line_dash`).

        Returns:
            bokeh.models.GlyphRenderer: The Bokeh glyph renderer object for
                the line plot.
        """
        # get plotting parameters
        seqlen = len(self.sequence)
        
        # transform the data
        tval = source[x].values / seqlen * 2 * np.pi
        rval = source[y].values
        
        # normalize the radius
        rnorm = (rval - np.min(rval)) / (np.max(rval) - np.min(rval))
        Rrange = self.radius2 - self.radius3
        
        # compute line
        rval = rnorm * Rrange + self.radius3
        xy = self.polar2xy(tval, rval)
        
        # format the dataframe
        df = pd.DataFrame(xy, columns=["xs", "ys"])
        
        # Set default plot properties
        plot_kwargs = {
            "line_color": "black",
            "line_dash": "solid",
            "line_width": 2,
            "alpha": 1,
        }
        # Allow kwargs to override defaults
        plot_kwargs.update(kwargs)
        
        glyph = self.p.line(source=df, x="xs", y="ys", **plot_kwargs)
        return glyph
