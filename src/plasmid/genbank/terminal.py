#!/usr/bin/env python
"""Class functions related to terminal sequence viewing and colorization."""

from typing import List, Dict, Optional, Tuple, Any, Protocol
import numpy as np
import matplotlib.colors as mcolors

# ==========================================
# 1. STANDALONE TERMINAL COLOR UTILITIES
# ==========================================

def get_ansi_fg(cname: Optional[str]) -> List[str]:
    """Get ANSI foreground escape codes."""
    if cname is None:
        return ['', '']
    try:
        # mcolors.to_rgb natively handles hex, named colors, and tuples
        r, g, b = mcolors.to_rgb(cname)
        r, g, b = int(r * 255), int(g * 255), int(b * 255)
        return [f"\033[38;2;{r};{g};{b}m", "\033[39m"]
    except (ValueError, TypeError):
        return ['', '']

def get_ansi_bg(cname: Optional[str]) -> List[str]:
    """Get ANSI background escape codes."""
    if cname is None:
        return ['', '']
    try:
        r, g, b = mcolors.to_rgb(cname)
        r, g, b = int(r * 255), int(g * 255), int(b * 255)
        return [f"\033[48;2;{r};{g};{b}m", "\033[49m"]
    except (ValueError, TypeError):
        return ['', '']

def get_colored_text(text: str, fg: Optional[str] = None, bg: Optional[str] = None) -> str:
    """Apply ANSI colors to a string."""
    sfg = get_ansi_fg(fg)
    sbg = get_ansi_bg(bg)
    return f"{sfg[0]}{sbg[0]}{text}{sfg[1]}{sbg[1]}"

def split_carr(fgarr: np.ndarray, bgarr: np.ndarray) -> List[Tuple[int, int]]:
    """Find contiguous blocks of identical colors."""
    if len(fgarr) == 0:
        return []
    
    changes = np.where((fgarr[1:] != fgarr[:-1]) | (bgarr[1:] != bgarr[:-1]))[0] + 1
    starts = np.insert(changes, 0, 0)
    stops = np.append(changes, len(fgarr))
    
    return list(zip(starts.tolist(), stops.tolist()))

def color_by_array(text: str, colordict: Dict[int, str], 
                   fgarr: Optional[np.ndarray] = None, 
                   bgarr: Optional[np.ndarray] = None) -> str:
    """Apply color arrays to format a text string."""
    if fgarr is None or len(fgarr) == 0:
        fgarr = np.zeros(len(text), dtype=int)
    if bgarr is None or len(bgarr) == 0:
        bgarr = np.zeros(len(text), dtype=int)
        
    out = []
    cuts = split_carr(fgarr, bgarr)
    
    for start, stop in cuts:
        fg = colordict.get(fgarr[start])
        bg = colordict.get(bgarr[start])
        
        if fg is None and bg is None:
            out.append(text[start:stop])
        else:
            out.append(get_colored_text(text[start:stop], fg, bg))
            
    return "".join(out)


# ==========================================
# 2. STREAMLINED GRAPHIC CLASS
# ==========================================

class Viewer:
    """Renders plasmid features as colored terminal text.
    
    Attributes:
        features: DataFrame containing plasmid feature data.
        seq: String representation of the plasmid sequence.
    """
    
    def __init__(self, rec):
        """Initialize with genbank record data."""
        self.features = rec.to_dataframe(str_only=False)
        self.seq = str(rec)
    
    def __str__(self) -> str:
        """Return the colorized plasmid string."""
        return 'Colorized Plasmid \n' + self.colorize()
    
    def colorize(self, width: int = 60) -> str:
        """Return the sequence with ANSI color formatting."""
        carr, cdict = self.get_color_array()
        seq_len = len(self.seq)
        
        if width == 0 or seq_len == 0:
            return color_by_array(self.seq, cdict, carr)
        
        start_indices = np.arange(0, seq_len, width)
        out = []
        
        for start in start_indices:
            # Cap stop at seq_len to prevent incorrect coordinates on the last chunk
            stop = min(start + width, seq_len)
            text_chunk = self.seq[start:stop]
            fgarr_chunk = carr[start:stop]
            
            out.append(color_by_array(text_chunk, cdict, fgarr_chunk))
            # Print coordinate on a new line to prevent ANSI terminal wrapping issues
            out.append(f"\n{stop:>6}\n") 
            
        return "".join(out)
    
    def get_color_array(self) -> Tuple[np.ndarray, Dict[int, str]]:
        """Map feature locations to a color index array."""
        seq_len = len(self.seq)
        if seq_len == 0:
            return np.array([], dtype=int), {}
        
        df = self.features.sort_values(by=['start'])
        cdict: Dict[int, str] = {}
        idict: Dict[str, int] = {}
        j = 1
        
        carr = np.zeros(seq_len, dtype=int)
        
        for loc, color in df[['location', 'color']].to_numpy():
            if color not in idict:
                cdict[j] = color
                idict[color] = j
                j += 1
                
            val = idict[color]
                        
            for start, end, strand in loc.positions:
                # Handle circular plasmids and out-of-bounds coordinates
                start = start % seq_len
                end = end % seq_len
                
                if start < end:
                    carr[start:end] = val
                else:
                    # Feature wraps around the origin
                    carr[start:] = val
                    carr[:end] = val
                
        return carr, cdict