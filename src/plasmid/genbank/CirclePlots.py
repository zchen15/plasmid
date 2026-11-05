# Module for browsing and editing genbank files

# science libs
import numpy as np
import pandas as pd
from .LinearPlots import LinearLayout

def polar2xy(pos: list, R: list) -> np.ndarray:
    """Compute Cartesian (x,y) positions from polar coordinates.
    
    Converts a list of angular positions and radii into Cartesian
    coordinates, adjusting the angle to start from the top (12 o'clock).

    Args:
        pos (list or numpy.ndarray): List of radial (angular) positions.
        R (list or numpy.ndarray): List of radii corresponding to the
            positions.

    Returns:
        numpy.ndarray: An array of shape (N, 2) containing the x and y
        coordinates.
    """
    R = np.array(R)
    t = -np.array(pos) + np.pi/2
    return np.transpose([R*np.cos(t), R*np.sin(t)])

def get_arc_xy(start, end, r1, r2, seqlen, N=None):
    """Compute arc coordinates from polar to Cartesian.
    
    Generates a set of x,y coordinates representing an arc from a starting
    position and radius to an ending position and radius.

    Args:
        start (int or float): Starting nucleotide position.
        end (int or float): Ending nucleotide position.
        r1 (float): Radius of the starting coordinate.
        r2 (float): Radius of the ending coordinate.
        seqlen (int): Total length of the sequence.
        N (int, optional): Number of data points to generate for the arc.
            Defaults to the absolute difference between end and start.

    Returns:
        numpy.ndarray: An array of shape (N, 2) containing the x and y
        coordinates.
    """
    if N is None:
        N = abs(int(end - start))

    t = np.linspace(start, end, N)/seqlen*2*np.pi
    R = np.linspace(r1, r2, N)
    return polar2xy(t, R)

def compute_label_positions(df: pd.DataFrame, arrow_width: float, x0: float, y0: float, seqlen: int) -> pd.DataFrame:
    """Get x,y coordinates for feature labels.
    
    Calculates the radial and angular positions for feature labels and
    generates the leader lines connecting the features to their labels.

    Args:
        df (pd.DataFrame): DataFrame containing feature information.
        arrow_width (float): Width of the arrows.
        x0 (float): X coordinate of the plot center.
        y0 (float): Y coordinate of the plot center.
        seqlen (int): Total length of the sequence.

    Returns:
        pd.DataFrame: DataFrame with computed label coordinates and 
        alignment information.
    """
    n = len(df)
    
    if n == 0:
        df["Lx"] = []
        df["Ly"] = []
        df["Tx"] = []
        df["Ty"] = []
        df["text_alignment"] = []
        return df

    # 1. Compute label angular positions (vectorized)
    starts = df["start"].to_numpy()
    ends = df["end"].to_numpy()
    t1 = (df["start"] + df["length"] / 2 + df["label_offset"]).to_numpy()
    
    t0 = t1.copy()
    t0[t1 > ends] = ends[t1 > ends]
    t0[t1 < starts] = starts[t1 < starts]
    t1_mod = t1 % seqlen

    df["text_t0"] = t0
    df["text_t1"] = t1_mod

    # 2. Compute label level: O(N log N) overlap counting
    valid = ends >= starts
    s = np.sort(starts[valid])
    e = np.sort(ends[valid])
    
    c1 = np.searchsorted(s, t0, side="right") - np.searchsorted(e, t0, side="left")
    c2 = np.searchsorted(s, t1_mod, side="right") - np.searchsorted(e, t1_mod, side="left")
    
    plotting_level = df["plotting_level"].to_numpy()
    label_level = np.maximum(np.maximum(c1, c2), plotting_level)
    df["label_level"] = label_level

    # 3. Compute label radial distance
    radius = df["radius"].to_numpy()
    text_r0 = radius + arrow_width / 2
    text_r1 = radius + (label_level - plotting_level) * arrow_width + arrow_width / 2

    # 4. Compute label x,y positions in a single vectorized pass
    pos = np.empty(2 * n)
    pos[0::2] = t0
    pos[1::2] = t1_mod
    
    rad = np.empty(2 * n)
    rad[0::2] = text_r0
    rad[1::2] = text_r1

    # Use the module-level polar2xy function
    theta = pos / seqlen * 2 * np.pi
    z = polar2xy(theta, rad).reshape(n, 2, 2)
    
    z[:, :, 0] += x0
    z[:, :, 1] += y0

    df["Lx"] = list(z[:, :, 0])
    df["Ly"] = list(z[:, :, 1])
    df["Tx"] = z[:, 1, 0]
    df["Ty"] = z[:, 1, 1]

    # 5. Compute text alignment vectorized
    df["text_alignment"] = np.where(z[:, 1, 0] > x0, "left", "right")
    
    return df

class CircleLayout(LinearLayout):
    """Calculates GenBank annotation coordinates for circular plots.

    This class holds functions pertaining to calculating and formatting
    coordinates for visualizing sequence features in a circular layout.
    """

    def __init__(self, data):
        """Initialize the circular layout.

        Args:
            data (Bio.SeqRecord or similar): Record object containing the
                sequence and feature data.
        """
        # Plot dimensions and parameters
        self.height = 600
        self.width = 650
        self.text_offset = 0
        self.tick_interval = 500
        self.tick_width_ratio = 0.02
        self.arrow_width_ratio = 0.02
        self.arrow_length_ratio = 0.02 * 0.5
        self.arrow_alpha = 1
        self.xy0 = [0, 0]
        
        if len(data.features) == 0:
            data = data.annotate(name="vector", seq=data.sequence)
        self.seqlen = len(data.sequence)
        self.features = data.features.copy()
        self.update()
        
    def update(self):
        """Apply current plot parameters to the object data.
        
        Calculates and sets arrow lengths, widths, plotting radii, tick
        intervals, and feature levels based on the sequence length and
        current parameters.
        """
        # set arrow length and width
        self.arrow_width = self.seqlen * self.arrow_width_ratio
        self.arrow_length = self.seqlen * self.arrow_length_ratio
        self.tick_width = self.seqlen * self.tick_width_ratio
        
        # set plotting radius
        tick_width = self.tick_width
        arrow_width = self.arrow_width

        self.radius0 = self.seqlen/(2*np.pi)
        self.radius1 = self.radius0 - arrow_width
        self.radius2 = self.radius1 - tick_width
        self.radius3 = self.radius2 - tick_width
        
        self.features["radius"] = self.radius0
        self.features["alpha"] = self.arrow_alpha
        
        # set tick interval
        self.set_tick_interval()
        
        # apply start and end positions
        self.features["fidx"] = range(len(self.features))
        out = []
        for fidx, loc in self.features[["fidx", "location"]].values:
            i = 0
            while i < len(loc.positions):
                [start1, stop1, strand1] = loc.positions[i]
                if stop1==self.seqlen and i < len(loc.positions)-1:
                    [start2, stop2, strand2] = loc.positions[i+1]
                    if strand2 == strand1 and start2 == 0:
                        stop1 = stop2
                        i+=1
                out.append([fidx, start1, stop1, strand1])
                i+=1
        df = pd.DataFrame(out, columns = ["fidx", "start", "end", "strand"])
        
        # correction for circularity if there is
        df["length"] = df["end"] - df["start"]
        c = df["end"] < df["start"]
        df.loc[c,"length"]+= self.seqlen
        
        self.features = self.features.merge(df, on="fidx", how="left")
        self.features["location"] = [str(i) for i in self.features["location"].values]
        
        # set plotting levels for each feature
        self.set_level()
        self.features["radius"]+= self.features["plotting_level"] * arrow_width
        self.features["label_offset"] = 0
        
    def set_tick_interval(self):
        """Set the tick interval for the circular plot.

        Chooses a "nice" round interval (1, 2, 2.5, or 5 × 10^n) that
        produces roughly 12-24 ticks around the circle, ensuring angular
        spacing stays above the minimum readable threshold (~0.2 rad).

        Returns:
            int: The calculated tick interval.
        """
        # Target ~16 ticks around the circle for readability.
        # Angular spacing = 2π / 16 ≈ 0.39 rad, well above the 0.2 rad floor.
        target_ticks = 16
        raw = self.seqlen / target_ticks
        
        # Very short sequences: just use 1 bp per tick
        if raw < 1:
            self.tick_interval = 1
            return self.tick_interval
        
        # Find the magnitude (power of 10) of the raw interval
        magnitude = 10 ** int(np.floor(np.log10(raw)))
        normalized = raw / magnitude
        
        # Snap up to the nearest "nice" multiplier
        if normalized <= 1:
            nice = 1
        elif normalized <= 2:
            nice = 2
        elif normalized <= 2.5:
            nice = 2.5
        elif normalized <= 5:
            nice = 5
        else:
            nice = 10
        
        self.tick_interval = int(nice * magnitude)
        return self.tick_interval

    def get_arc_ticks(self):
        """Get x,y coordinates for arc ticks and their labels.
        
        Calculates the positions for tick marks around the circle and their
        corresponding text labels.

        Returns:
            pandas.DataFrame: A DataFrame containing tick coordinates, text
            positions, and text alignment information.
        """
        tick_interval = self.tick_interval
        x0, y0 = self.xy0

        # generate ticks at interval
        t = np.arange(0, self.seqlen, tick_interval)
        if len(t) > 1 and (1 - t[-1] / self.seqlen) < (t[1] / self.seqlen):
            t = t[:-1]
            
        df = pd.DataFrame(t, columns=["start"])
        df["name"] = df["start"]
        n = len(df)

        # calculate radial coordinates as local variables only
        r0 = self.radius0 - self.arrow_width / 2
        r1 = r0 - self.tick_width / 2
        r2 = r1 - self.tick_width / 4

        # compute all angles at once
        theta = t / self.seqlen * 2 * np.pi

        # tick marks: 2 points per tick (r0 -> r1) at a fixed angle
        # np.repeat repeats each angle twice, np.tile repeats the radii array n times
        z = polar2xy(np.repeat(theta, 2), np.tile([r0, r1], n))
        z = z.reshape(n, 2, 2)
        
        df["Lx"] = list(z[:, :, 0] + x0)
        df["Ly"] = list(z[:, :, 1] + y0)

        # tick labels: 1 point per tick at r2
        z2 = polar2xy(theta, np.full(n, r2))
        df["Tx"] = z2[:, 0] + x0
        df["Ty"] = z2[:, 1] + y0

        # compute text alignment vectorized
        align = np.full(n, "right", dtype=object)
        align[t > self.seqlen / 2] = "left"
        if n:
            align[0] = "center"
        df["text_alignment"] = align

        return df

    def get_label_positions(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Wrapper to compute label positions using the module-level function.

        Args:
            df (pd.DataFrame): DataFrame containing feature information.

        Returns:
            pd.DataFrame: DataFrame with computed label coordinates and 
            alignment information.
        """
        df = compute_label_positions(df, self.arrow_width, self.xy0[0], self.xy0[1], self.seqlen)
        return df

    def get_scatter_line(self, df):
        """Get x,y coordinates for linear segments on a circular plot.
        
        Converts alignment or read mapping information into Cartesian
        coordinates for plotting as lines or chords across the circular plot.

        Args:
            df (pandas.DataFrame): DataFrame with alignment information,
                including query/target IDs, start/end positions, and
                orientation.
        
        Returns:
            pandas.DataFrame: A DataFrame with "xs", "ys", and metadata
            columns for plotting the lines.
        """
        radius = self.radius3
        x0, y0 = self.xy0

        # 1. Assign start and end positions based on orientation (vectorized)
        c = df["orientation"] == "+"
        df.loc[c, "start"] = df.loc[c, "t_start"].values
        df.loc[c, "end"] = df.loc[c, "t_end"].values
        df.loc[~c, "start"] = df.loc[~c, "t_end"].values
        df.loc[~c, "end"] = df.loc[~c, "t_start"].values

        # 2. Sort once
        df = df.sort_values(by=["q_id", "q_start", "q_end"]).reset_index(drop=True)

        # 3. Create reference mapping for all unique t_ids at once
        all_t_ids = df["t_id"].unique()
        lets = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
        refs = {tid: lets[i % len(lets)] for i, tid in enumerate(all_t_ids)}
        revref = str({refs[k]: k for k in refs.keys()})

        # 4. Map references to the dataframe
        df["ref_letter"] = df["t_id"].map(refs)

        # 5. Build location and read_loc string parts using vectorized Series concatenation
        t_start_str = df["t_start"].astype(str)
        t_end_str = df["t_end"].astype(str)
        q_start_str = df["q_start"].astype(str)
        q_end_str = df["q_end"].astype(str)
        
        df["loc_part"] = df["ref_letter"]+"[" + t_start_str + ":" + t_end_str + "]" + df["orientation"]
        df["read_loc_part"] = df["ref_letter"]+"[" + q_start_str + ":" + q_end_str + "]" + df["orientation"]
        
        # 6. Format numeric columns to strings vectorized (no lambdas)
        df["match_score_str"] = df["match_score"].map("{:.3f}".format)
        df["coverage_str"] = df["coverage"].map("{:.3f}".format)

        # 7. Group by q_id and aggregate strings using built-in methods (no lambdas)
        grouped = df.groupby("q_id")
        locations = grouped["loc_part"].agg("".join)
        read_locs = grouped["read_loc_part"].agg("".join)
        match_scores = grouped["match_score_str"].agg(", ".join)
        coverages = grouped["coverage_str"].agg(", ".join)

        # 8. Map aggregated strings back to the original rows
        df["location"] = df["q_id"].map(locations)
        df["read_loc"] = df["q_id"].map(read_locs)
        df["coverage"] = df["q_id"].map(coverages)
        df["match_score"] = df["q_id"].map(match_scores)
        df["location_key"] = revref

        # 9. Compute xy coordinates for all segments in a single vectorized pass
        starts = df["start"].to_numpy()
        ends = df["end"].to_numpy()

        # Interleave start and end positions: [s1, e1, s2, e2, ...]
        positions = np.empty(len(starts) * 2)
        positions[0::2] = starts
        positions[1::2] = ends

        theta = positions / self.seqlen * 2 * np.pi
        radii = np.full(len(positions), radius)
        xy_all = polar2xy(theta, radii)

        # Shift coordinates
        xy_all[:, 0] += x0
        xy_all[:, 1] += y0

        # Reshape to (N_segments, 2, 2) to extract xs and ys for each segment
        xy_segments = xy_all.reshape(-1, 2, 2)
        xs = list(xy_segments[:, :, 0])
        ys = list(xy_segments[:, :, 1])

        # 10. Construct the final DataFrame
        out_df = pd.DataFrame({
            "q_id": df["q_id"],
            "xs": xs,
            "ys": ys,
            "location": df["location"],
            "location_key": df["location_key"],
            "read_loc": df["read_loc"],
            "coverage": df["coverage"],
            "match_score": df["match_score"]
        })
        return out_df

    def get_arrow(self):
        """Get x,y coordinates for arc arrow patches.
        
        Generates the polygon coordinates for drawing sequence features as
        arrows along the circular backbone, handling both forward and reverse
        strands, as well as features that wrap around the origin.

        Returns:
            pandas.DataFrame: The feature DataFrame updated with "xs" and "ys"
            columns containing the arrow patch coordinates.
        """
        df = self.features
        w = self.arrow_width / 2
        h = self.arrow_length
        x0, y0 = self.xy0

        if len(df) == 0:
            df["xs"] = []
            df["ys"] = []
            return df

        # Extract columns once as numpy arrays to avoid repeated df.loc lookups
        starts = df["start"].to_numpy()
        ends = df["end"].to_numpy()
        strands = df["strand"].to_numpy()
        radii = df["radius"].to_numpy()
        n = len(df)

        # Accumulate (theta, radius) pairs for every vertex of every arrow
        t_chunks = []
        r_chunks = []
        sizes = np.zeros(n, dtype=int)

        for i in range(n):
            start = starts[i]
            end = ends[i]
            strand = strands[i]
            R = radii[i]

            # handle features that wrap around the origin
            if end < start:
                end = end + self.seqlen

            # Build the list of arc segments for this feature
            # Each segment is (t0, t1, r0, r1, N) where N=None means "full resolution"
            if strand == 1:
                mid = end - h
                if mid < start:
                    segments = [
                        (start, start, R, R + w, 2),
                        (end, start, R, R - w, 2),
                    ]
                else:
                    segments = [
                        (start, start, R, R + w, 1),
                        (start, mid, R + w, R + w, None),
                        (end, end, R, R - w, 1),
                        (mid, start, R - w, R - w, None),
                    ]
            else:
                mid = start + h
                if mid > end:
                    segments = [
                        (start, end, R, R + w, 2),
                        (end, end, R - w, R, 1),
                    ]
                else:
                    segments = [
                        (start, start, R, R + w, 1),
                        (mid, end, R + w, R + w, None),
                        (end, end, R, R - w, 1),
                        (end, mid, R - w, R - w, None),
                    ]

            # Convert each segment into linspace arrays
            count = 0
            for t0, t1, r0, r1, N in segments:
                if N is None:
                    N = abs(int(t1 - t0))
                t_chunks.append(np.linspace(t0, t1, N))
                r_chunks.append(np.linspace(r0, r1, N))
                count += N
            sizes[i] = count

        # Single polar2xy call for every vertex of every arrow
        if t_chunks:
            t_all = np.concatenate(t_chunks)
            r_all = np.concatenate(r_chunks)
            theta = t_all / self.seqlen * 2 * np.pi
            xy = polar2xy(theta, r_all)
        else:
            xy = np.empty((0, 2))

        # Split the flat xy array back into per-feature polygons
        cuts = np.cumsum(sizes)[:-1]
        xs_parts = np.split(xy[:, 0], cuts)
        ys_parts = np.split(xy[:, 1], cuts)

        # Shift by the plot origin and assign
        df["xs"] = [x + x0 for x in xs_parts]
        df["ys"] = [y + y0 for y in ys_parts]
        return df