import numpy as np
import pandas as pd
import warnings
from typing import Any, List

def compute_label_positions(df: pd.DataFrame, arrow_width: float, x0: float, y0: float) -> pd.DataFrame:
    """
    Calculates the x and y coordinates for linear feature labels.
    
    This function computes the placement of text labels along the sequence, 
    handling boundary clamping for out-of-bounds labels, calculating vertical 
    stacking levels to avoid overlap, and generating the final coordinate 
    arrays for rendering.

    Args:
        df (pd.DataFrame): DataFrame containing feature information, including 
            'start', 'end', 'length', 'plotting_level', 'label_x_offset', 
            and 'label_y_offset'.
        arrow_width (float): The width of the feature arrows, used to calculate 
            vertical spacing for label stacking.
        x0 (float): The global x-axis offset for the plot.
        y0 (float): The global y-axis offset for the plot.
    
    Returns:
        pd.DataFrame: The input DataFrame augmented with text alignment and 
        coordinate columns: 'text_t1', 'text_t0', 'label_level', 'text_r1', 
        'text_r0', 'Tx', 'Ty', 'Lx', and 'Ly'.
    """    
    # compute label positions along the sequence
    df["text_t1"] = df["start"] + df["length"]/2 + df["label_x_offset"]
    
    # Handle features that wrap around the origin (start > end)
    is_wrapping = df["start"] > df["end"]
    
    # set label out of bound positions
    c1 = (df["text_t1"] > df["end"]) & (~is_wrapping)
    c2 = (df["text_t1"] < df["start"]) & (~is_wrapping)
    c3 = (df["text_t1"] > df["end"]) & (df["text_t1"] < df["start"]) & is_wrapping
    
    df["text_t0"] = df["text_t1"]
    df.loc[c1, "text_t0"] = df.loc[c1, "end"]
    df.loc[c2, "text_t0"] = df.loc[c2, "start"]
    
    # For wrapping features in the "gap", clamp to the nearest boundary
    df.loc[c3, "text_t0"] = np.where(
        (df.loc[c3, "text_t1"] - df.loc[c3, "end"]) < (df.loc[c3, "start"] - df.loc[c3, "text_t1"]),
        df.loc[c3, "end"],
        df.loc[c3, "start"]
    )
    
    # compute label level
    label_level: List[int] = []
    for i in range(len(df)):
        idx = df["text_t0"].iat[i]
        c1 = np.sum((df["start"] <= idx) & (df["end"] >= idx))
        idx = df["text_t1"].iat[i]
        c2 = np.sum((df["start"] <= idx) & (df["end"] >= idx))
        label_level.append(int(np.max([c1, c2])))
        
    df["label_level"] = label_level
    df["text_r1"] = (df["label_level"]+df["label_y_offset"])*arrow_width + 1.5*arrow_width
    df["text_r0"] = df["plotting_level"]*arrow_width + arrow_width/2
    
    # compute label x,y positions
    df["Tx"] = df["text_t1"] + x0
    df["Ty"] = -df["text_r1"] + y0
    df["Lx"] = [[x0+t0, x0+t1] for t0, t1 in df[["text_t0","text_t1"]].values]
    df["Ly"] = [[y0-r0, y0-r1] for r0, r1 in df[["text_r0","text_r1"]].values]
    return df

class LinearLayout:
    """
    Calculates GenBank annotation coordinates for linear sequence plots.
    
    This class manages the layout logic for visualizing sequence features on 
    a linear axis. It handles feature stacking (to prevent overlaps), calculates 
    arrow geometries for forward/reverse strands, and determines label placements.
    """

    # Class-level type annotations for instance variables
    text_offset: int
    tick_interval: int
    tick_width_ratio: float
    arrow_width_ratio: float
    arrow_length_ratio: float
    arrow_alpha: float
    xy0: List[float]
    sequence: str
    features: pd.DataFrame
    arrow_width: float
    arrow_length: float

    def __init__(self, data: Any) -> None:
        """
        Initializes the linear layout with plot parameters and feature data.

        Args:
            data (Any): A sequence record object containing 'features' 
                (DataFrame or list) and 'sequence' (str) attributes. If the 
                record has no features, a default 'vector' feature is created.
        """
        # Initialize plot parameters as direct class variables
        self.text_offset = 0
        self.tick_interval = 500
        self.tick_width_ratio = 0.02
        self.arrow_width_ratio = 0.02
        self.arrow_length_ratio = 0.02 * 0.5
        self.arrow_alpha = 1.0
        self.xy0 = [0.0, 0.0]
        
        if len(data.features) == 0:
            data = data.annotate(name="vector", seq=data.sequence)
        self.sequence = data.sequence
        self.features = data.features.copy()
        self.update()
        
    def update(self) -> None:
        """
        Applies current plot parameters to the feature data.
        
        This method recalculates absolute arrow dimensions based on the sequence 
        length, resolves multi-part feature locations into individual rows, 
        assigns vertical stacking levels, and initializes label offsets.
        """
        # set arrow length and width
        L = len(self.sequence)
        self.arrow_width = float(L * self.arrow_width_ratio)
        self.arrow_length = float(L * self.arrow_length_ratio)
        self.features["alpha"] = self.arrow_alpha
        
        # apply start and end positions
        self.features["fidx"] = range(len(self.features))
        out: List[List[Any]] = []
        for fidx, loc in self.features[["fidx", "location"]].values:
            for start, stop, strand in loc.positions:
                out.append([fidx, start, stop, strand])
        df = pd.DataFrame(out, columns=["fidx", "start", "end", "strand"])
        
        # correction for circularity if there is
        df["length"] = df["end"] - df["start"]
        c = df["end"] < df["start"]
        df.loc[c, "length"] += len(self.sequence)
        
        self.features = self.features.merge(df, on="fidx", how="left")
        
        # Store string representation in a new column to preserve original Location objects
        self.features["location_str"] = [str(i) for i in self.features["location"].values]
        
        # set plotting levels for each feature
        self.set_level()
        
        # set label offset positions
        self.features["label_x_offset"] = 0
        self.features["label_y_offset"] = 0
        
    def set_level(self) -> None:
        """
        Calculates and assigns vertical stacking levels for each feature.
        
        Iterates through the features to determine the 'plotting_level' (y-axis 
        stack) for each, ensuring that overlapping features are placed on 
        different vertical levels. Uses an O(N^2) pairwise comparison.
        """
        df = self.features.copy()
        
        # assign plotting levels
        cols = ["length", "start", "fidx"]
        df = df[cols].sort_values(by=cols, ascending=[False, True, True]).reset_index(drop=True)
        df["plotting_level"] = 0
        
        for i in range(len(df)):
            q1 = df["start"].iat[i]
            q2 = df["end"].iat[i]
            qL = df["plotting_level"].iat[i]
            
            for j in range(i+1, len(df)):
                t1 = df["start"].iat[j]
                t2 = df["end"].iat[j]
                tL = df["plotting_level"].iat[j]

                c1 = tL == qL
                if q2 < q1 and t2 < t1:
                    c2 = True
                elif q2 < q1 or t2 < t1:
                    c2 = (q2 > t1) | (q1 < t2)
                else:
                    c2 = ~((q2 <= t1) | (t2 <= q1))
                    
                if c1 & c2:
                    df.loc[j, "plotting_level"] += 1
        
        self.features = self.features.merge(df[["fidx", "plotting_level"]], on="fidx", how="left")
        self.features = self.features.sort_values(by=["start", "plotting_level"]).reset_index(drop=True)
        
    def get_label_positions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Wrapper to compute label positions using the module-level function.

        Args:
            df (pd.DataFrame): DataFrame containing feature information.

        Returns:
            pd.DataFrame: DataFrame with computed label coordinates and 
            alignment information.
        """
        df = compute_label_positions(df, self.arrow_width, self.xy0[0], self.xy0[1])
        return df
        
    def get_arrow(self) -> pd.DataFrame:
        """
        Generates the x and y coordinate arrays for drawing linear feature arrows.
        
        Calculates the polygon vertices for each feature arrow, handling forward 
        and reverse strand orientations. Features that wrap around the origin 
        trigger a warning and are adjusted for linear plotting.

        Returns:
            pd.DataFrame: The features DataFrame augmented with 'xs' and 'ys' 
            columns, where each cell contains a list of x and y coordinates 
            defining the arrow polygon.
        """
        df = self.features.copy()
        L = len(self.sequence)
        w = self.arrow_width / 2
        h = self.arrow_length

        # iterate through each feature
        xs: List[np.ndarray] = []
        ys: List[np.ndarray] = []
        for i in range(len(df)):
            start = df["start"].iat[i]
            end = df["end"].iat[i]
            strand = df["strand"].iat[i]
            y_level = df["plotting_level"].iat[i] * w * 2

            # handle features that wrap around the origin
            if end < start:
                end += L
                warnings.warn(f"Feature index {df['fidx'].iat[i]} wraps around the origin. This is not a linear construct.")

            # draw top half of arrow in fwd orientation
            if strand == 1:
                mid = end - h
                if mid < start:
                    mid = start
                x = [start, start, mid, end]

            # draw top half of arrow in rev orientation
            else:
                mid = start + h
                if mid > end:
                    mid = end
                x = [start, mid, end, end]

            # set y coord of the arrow
            y = [0, w, w, 0]

            # add bottom half of arrow
            x = x + x[::-1]
            y = y + [-i for i in y[::-1]]

            # shift plotting axis
            x0, y0 = self.xy0
            x_arr = np.array(x) + x0
            y_arr = np.array(y) - y_level + y0

            xs.append(x_arr)
            ys.append(y_arr)
            
        df["xs"] = xs
        df["ys"] = ys
        return df