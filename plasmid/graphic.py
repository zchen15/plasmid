#!/usr/bin/env python
# Class functions related to plotting

# computing libraries
import numpy as np

# plotting libraries
import matplotlib as mpl
import matplotlib.pyplot as plt
import bokeh
import bokeh.plotting
import bokeh.palettes
import sty

from .misc import unpack_FeatureLocation

class Graphic:
    '''
    The class contains functions for drawing plasmids and plotting things
    '''
    params = {'colorlib':[mpl._color_data.XKCD_COLORS, mpl._color_data.BASE_COLORS, mpl._color_data.CSS4_COLORS],
              'plots':{'R':1000,
                       'W_ratio':4,
                       'x':0,
                       'y':0,
                       'res':100,
                       'level':2,
                       'alpha':0.5,
                       'add_index':False}}

    def __init__(self, data=None):
        '''
        Initializes graphic object
        data = plasmid dataframe containing sequence and features
        R = Radius of the plasmid. Default none
        '''
        if data!=None:
            self.Plasmid = data
            self.features = data.to_dataframe(str_only=False)
            self.params['nt'] = len(data.__str__())  # total length of the plasmid
            self.params['plots']['R'] = self.params['nt']/(2*np.pi)
            self.auto_level()

    def auto_level(self):
        '''
        Auto feature level sort and assign by overlap
        # debug preserve feature order
        '''
        N = len(self.features)
        # go in order of largest features
        self.features['len'] = [len(i) for i in self.features['location']]
        self.features = self.features.sort_values(by=['len'], ascending=False)
        level = [0]*N
        for i in range(N):
            loc1 = self.features.iloc[i]['location']
            lvl1 = level[i]
            for j in range(i+1,N):
                loc2 = self.features.iloc[j]['location']
                lvl2 = level[j]
                if Graphic.isOverlapping(loc1, loc2) and level[j] <= level[i]:
                    level[j] = level[i]+1
        self.features['level'] = level
        self.features = self.features.sort_values(by=['start','level'])

    def get_bp_overlap(loc1, loc2):
        '''
        Get bp overlap between two features
        loc1 = feature location of feature1
        loc2 = feature location of feature2
        '''
        loc1 = unpack_FeatureLocation(loc1)
        loc2 = unpack_FeatureLocation(loc2)
        bp = 0
        for a in loc1:
            for b in loc2:
                c1 = max(a.start,b.start)
                c2 = min(a.end,b.end)
                L = c2 - c1
                if L > 0:
                    bp+=L
        return bp

    def isOverlapping(loc1, loc2):
        '''
        Check if features are overlapping
        '''
        loc2 = unpack_FeatureLocation(loc2)
        for L in loc2:
            if loc1.__contains__(L.start) or loc1.__contains__(L.end):
                return True
        return False

    def arc_x(self, t):
        return np.cos(t/self.nt*2*np.pi)
    
    def arc_y(self, t):
        return np.sin(t/self.nt*2*np.pi)

    def arc_xy(self, t):
        return np.transpose([self.arc_x(t), self.arc_y(t)])

    def arc_arrow(self, R, W, start, end, strand):
        '''
        Compute arc arrow polygons
        R = distance from center to draw the arc
        W = width of arc
        L = length of the arc in nucleotides
        strand = direction of arrow
        start = starting nucleotide position
        end = end nucleotide position
        return set of points defining arc arrow polygon
        '''
        # handle coordinates that wrap around the origin
        if end < start:
            end+= self.nt
            t = np.arange(0, end, self.res)
        else:
            t = np.arange(0, self.nt, self.res)
        # compute start and end lines
        D = self.nt*W/2/(R+W/2)/2/np.pi
        if strand > 0: 
            if end - D < start:
                D = end - start
            p_start = np.stack([R*self.arc_xy(start), (R+W)*self.arc_xy(start)])
            p_end = np.stack([(R+W)*self.arc_xy(end-D), (R+W/2)*self.arc_xy(end), R*self.arc_xy(end-D)])
        # compute for reverse direction
        else:
            if start+D > end:
                D = end - start
            p_end = np.stack([(R+W)*self.arc_xy(end), R*self.arc_xy(end)])
            p_start = np.stack((R*self.arc_xy(start+D), (R+W/2)*self.arc_xy(start), (R+W)*self.arc_xy(start+D)))
        # get inner and out arcs
        if strand > 0:
            t = t[(t > start) & (t < end-D)]
        else:
            t = t[(t > start+D) & (t < end)]
        p_outer = (R+W)*self.arc_xy(t)
        p_inner = R*self.arc_xy(t[::-1])
        # return the arc arrow polygon
        out = [p_start, p_outer, p_end, p_inner]
        return np.concatenate(out, axis=0)

    def add_arc_arrow(self, p):
        '''
        Adds arc arrow glyph to current plot
        '''
        data = self.features
        xs = []
        ys = []
        for start, end, strand, R, W in data[['start','end','strand','R','W']].values:
            # compute arch
            xy = self.arc_arrow(R, W, start, end, strand)
            xs.append(xy[:,0]+self.x)
            ys.append(xy[:,1]+self.y)
        data['xs'] = xs
        data['ys'] = ys
        source = bokeh.models.ColumnDataSource(data)
        # add glyph
        glyph = bokeh.models.Patches(xs='xs', ys='ys', fill_color='color', fill_alpha='alpha')
        p.add_glyph(source, glyph)
        return p
    
    def add_ticks(self):
        '''
        Add tick marks
        '''
    
    def annotate(self, p):
        '''
        annotate features
        '''
        data = self.features
        xs = []
        ys = []
        xy = []
        ts = []
        for start, end, R, W, level in data[['start','end','R','W','level']].values:
            # compute the line
            if end < start:
                end+= self.nt
            t = (start + end)/2
            pos1 = (R+W)*self.arc_xy(t)
            pos2 = (R+level*W)*self.arc_xy(t)
            xy.append(pos2)
            xs.append([pos1[0], pos2[0]])
            ys.append([pos1[1], pos2[1]])
            ts.append(t)
        xy = np.array(xy)
        data['x'] = xy[:,0] + self.x
        data['y'] = xy[:,1] + self.y
        data['xs'] = xs
        data['ys'] = ys
        data['ts'] = ts
        # add the line
        glyph = bokeh.models.MultiLine(xs='xs', ys='ys', line_color='black', line_width=1, line_alpha=1)
        # add index to labels 
        if self.add_index:
            for i in range(len(data)):
                data.loc[i,'locus_tag'] = str(i) + ' ' + data.iloc[i]['locus_tag']
        source = bokeh.models.ColumnDataSource(data)
        p.add_glyph(source, glyph)
        # add text
        t = (data['ts']%self.nt)/self.nt
        c = (t > 0.25) & (t < 0.75)
        source = bokeh.models.ColumnDataSource(data[~c])
        glyph = bokeh.models.Text(x='x', y='y', text='locus_tag', angle=0, text_color='color',
                            text_align='left', text_baseline='middle')
        p.add_glyph(source, glyph)
        source = bokeh.models.ColumnDataSource(data[c])
        glyph = bokeh.models.Text(x='x', y='y', text='locus_tag', angle=0, text_color='color',
                            text_align='right', text_baseline='middle')
        p.add_glyph(source, glyph)

    def bokeh_plot(self, p=None, width=400, height=350):
        '''
        Plot features using bokeh
        p = bokeh plot reference
        '''
        tooltips = [('locus_tag', '@label'),
                    ('type', '@type'),
                    ('location', '@location'),
                    ('length','@length'),
                    ('color','@color')]

        # create figure object if it does not exist
        if p==None:
            p = bokeh.plotting.figure(tooltips=tooltips)
            p.plot_width = width
            p.plot_height = height
        self.add_arc_arrow(p)
        self.annotate(p)
        return p 

    def plot(self, p=None):
        '''
        Plot features using matplotlib
        p = matplotlib figure handle
        '''

    def __repr__(self):
        loc = str(type(self))+' at '+hex(id(self))+'\n'
        loc+= self.__str__()
        return str(loc)

    def __str__(self):
        '''
        Format graphic record text as string
        '''
        out = 'Colorized Plasmid \n'
        out+= self.colorize_plasmid()
        return out

    def colorize_plasmid(self, width=0):
        '''
        Return text denoting colored plasmid
        width = chars per line
        '''
        carr, cdict = self.get_color_array()
        seq = self.Plasmid.__str__()
        if width==0:
            return self.color_by_array(seq, cdict, carr)
        else:
            start = np.arange(0, len(seq), width)
            stop = start + width
            stop[stop > len(seq)] = len(seq)
        # print out plasmid
        out = ''
        for i in range(len(start)):
            text = seq[start[i]:stop[i]]
            fgarr = carr[start[i]:stop[i]]
            out+= self.color_by_array(text, cdict, fgarr)
            out+= '   ' + str(stop[i]) + '\n'
        return out

    def get_color_array(self):
        '''
        Returns integer array and dictionary mapping integer to color values
        '''
        df = self.features.sort_values(by=['level','start'])
        cdict = {}
        idict = {}
        j = 1
        seq = self.Plasmid.__str__()
        carr = np.zeros(len(seq))
        for loc, color in df[['location','color']].values:
            if (color in idict.keys())==False:
                cdict[j] = color
                idict[color] = j
                j+= 1
            val = idict[color]
            loc = unpack_FeatureLocation(loc)
            for L in loc:
                carr[L.start:L.end] = val
        return carr, cdict

    def color_by_array(self, text, colordict, fgarr=[], bgarr=[]):
        '''
        Plot color text marked by an array
        text = text to be colorized
        colordict = dictionary mapping integers to colors
        fgarr = integer array of text color
        bgarr = integer array of text background color
        '''
        if len(fgarr)==0:
            fgarr = np.zeros(len(text))
        if len(bgarr)==0:
            bgarr = np.zeros(len(text))
        out = ''
        cuts = self.split_carr(fgarr, bgarr)
        for start, stop in cuts:
            fg = None
            bg = None
            if fgarr[start] in colordict.keys():
                fg = colordict[fgarr[start]]
            if bgarr[start] in colordict.keys():
                bg = colordict[bgarr[start]]
            out+= self.get_colored_text(text[start:stop], fg, bg)
        return out
    
    def split_carr(self, fgarr, bgarr):
        '''
        Split up blocks of text by their coloring scheme
        '''
        start = 0
        val = (fgarr[1:]!=fgarr[:-1]) | (bgarr[1:]!=bgarr[:-1])
        val = np.arange(len(val))[val] + 1
        cuts = []
        start = 0
        for i in range(len(val)):
            cuts.append([start, val[i]]) 
            start = val[i]
        cuts.append([start, len(fgarr)])
        return cuts

    def get_colored_text(self, text, fg=None, bg=None):
        '''
        Format string into colored text
        fg = text color
        bg = background color
        '''
        sfg = self.get_sty_fg(fg)
        sbg = self.get_sty_bg(bg)
        return sfg[0] + sbg[0] + text + sfg[1] + sbg[1]
    
    def get_sty_fg(self, cname):
        '''
        Get sty fg markup
        cname = name of color
        '''
        x = self.color_to_int(cname)
        if x==None:
            return ['', '']
        else:
            return [sty.fg(x[0],x[1],x[2]), sty.fg.rs]

    def get_sty_bg(self, cname):
        '''
        Get sty bg markup
        cname = name of color
        '''
        x = self.color_to_int(cname)
        if x==None:
            return ['','']
        else:
            return [sty.bg(x[0],x[1],x[2]), sty.bg.rs]

    def color_to_int(self, cname):
        '''
        Get color name to 24 bit values
        cname = name of color
        return integer array
        '''
        if type(cname)!=str:
            return None
        if cname[0]!='#':
            cname = self.color_to_hex(cname)
        return self.hex_to_int(cname)

    def color_to_hex(self, cname):
        '''
        Get hex color from color name
        '''
        colorlib = self.params['colorlib']
        for clib in colorlib:
            if cname in clib.keys():
                return clib[cname] 

    def hex_to_int(self, cnum):
        '''
        Convert hex color number to integers
        '''
        a = cnum[1:3]
        b = cnum[3:5]
        c = cnum[5:7]
        return [int(i,16) for i in [a,b,c]]
   
    def get_colormap(name, palette):
        '''
        Generates a colormap for a set of dna features using colors from a palette
        name = list of names of the features
        palette = list of colors
        return dictionary of {name:color}
        '''
        if type(palette)==str:
            palette = bokeh.palettes.__dict__[palette]
        N = len(palette)/2
        cmap = {}
        for i in range(len(name)):
            fwd_color = palette[(i*2)%len(palette)]
            rev_color = palette[(i*2+1)%len(palette)]
            cmap[name[i]] = [fwd_color, rev_color]
        return cmap