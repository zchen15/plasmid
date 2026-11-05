# Module for browsing and editing genbank files

# science libs
import numpy as np

# plotting libraries
import bokeh.plotting
import bokeh.models

from ..LinearPlots import LinearLayout

class Plots(LinearLayout):
    """
    This class holds functions pertaining to visualizing genbank annotation
    information in a linear format.

    Layout and coordinate calculations are provided by
    :class:`LinearPlots.LinearLayout`.
    """
    def autoscale_range(self):
        """
        Auto scale plot ranges for a linear plot
        
        Args:
            p: bokeh plotting object

        Returns:
            bokeh plotting object
        """
        arrow_width = self.arrow_width
        hw_ratio = p.height/p.width
        R = len(self.sequence)

        L = arrow_width*10
        w = R + L
        h = w * hw_ratio

        h_min = -1.2*np.min(self.features["Ty"])
        if h < h_min:
            h = h_min
            w = h/hw_ratio

        [x0,y0] = self.xy0
        L = w - R
        p.x_range = bokeh.models.Range1d(x0 - L/2, x0 + w - L/2)
        p.y_range = bokeh.models.Range1d(y0 - h + arrow_width, y0 + arrow_width)
        return p

    def add_features(self):
        """
        Draw features as linear arrows
        
        Args:
            p: bokeh.plotting figure handle
            legends: whether to add feature legends
            xy0: optional plotting origin

        Return
            bokeh.plotting figure handle
        """
        [x0, y0] = self.xy0

        # draw the backbone line
        L = len(self.sequence)
        x = [x0, x0+L]
        y = [y0, y0]
        p.line(x, y, color="black", line_width=2)
        
        # plot arrows as features
        df = self.get_arrow()
        arcs = p.patches(source=df, xs="xs", ys="ys", fill_color="color", fill_alpha="alpha",
                         line_color="black", line_width=3)
        
        # enable hovertools for arcs only
        TOOLTIPS = [("name", "@name"),
                    ("feature","@feature"),
                    ("color","@color"),
                    ("location", "@location")]
        hvr = bokeh.models.HoverTool(tooltips=TOOLTIPS)
        hvr.renderers = [arcs]
        p.add_tools(hvr)
        return p
        
    def add_labels(self):
        """
        Draw labels to features
        
        Return
            bokeh.plotting figure handle
        """
        # draw feature label ticks
        df = self.get_label_positions(self.features)
        p.multi_line(source=df, xs="Lx", ys="Ly", line_color="color", line_width=2, alpha="alpha")
        
        # draw feature label text
        p.text(source=df, x="Tx", y="Ty", x_offset=0, y_offset=20, text_align="center", text="name")
        return p
        