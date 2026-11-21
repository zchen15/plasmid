# Custom class for storing feature locations
import numpy as np
import copy

class Location:
    """
    This class holds function related to merging and storing locations

    Attributes:
        positions: unique positions stored by the location.
    """
    positions = np.array([], dtype=int)

    def __init__(self, pos: list):
        """
        Initializes the Location object

        Args:
            pos: list of positions.
            kwargs: settings that override the class defaults.

        Returns:
            Location object
        """
        self.positions = np.unique(pos)
                
    def get_slices(self):
        """
        Transform list of indices into list of slices

        Returns:
            list of slices
        """        
        c0 = self.positions[:-1] != self.positions[1:] - 1
        c1 = np.concatenate(([True], c0))
        c2 = np.concatenate((c0, [True]))
        s1 = self.positions[c1]
        s2 = self.positions[c2]
        
        out = []
        for i in range(len(s1)):
            out.append(slice(int(s1[i]), int(s2[i])+1))
        return out
        
    def remap(self, pos:dict):
        """
        Remap atom_id
        """
        positions = [pos[i] for i in self.positions]
        return Location(positions)
        
    def contains(self, pos:list):
        """
        check if location contains all positions
        """
        return all(np.isin(self.positions, pos))
        
    def intersection(self, pos:list):
        """
        keep only positions that intersect
        """
        pos = np.intersect1d(self.positions, pos)
        return Location(pos)

    def get_positions(self):
        return self.positions
    
    def __len__(self):
        """
        Returns length of position indicies
        """
        return len(self.positions)
        
    def __str__(self):
        """
        Format list of indices into a string
        """
        out = [f"[{s.start}:{s.stop}]" for s in self.get_slices()]
        out = "+".join(out)
        return out

    def __repr__(self):
        """
        Returns the representation of the string
        """
        return self.__str__()
        
    def __add__(self, other):
        """
        Performs self + other
        
        Args:
            other: Record object
        
        Returns:
            modified copy of object
        """
        # check if this should be inplace or not
        if isinstance(other, Location):
            pos = np.concatenate((self.positions, other.positions))
            pos = np.unique(pos)
        elif isinstance(other, int):
            pos = self.positions + other
        return Location(pos)
        
    def __radd__(self, other):
        """
        Performs other + self
        
        Args:
            other: Record object
        
        Returns:
            modified copy of object
        """
        out = self.__add__(other)
        return out
        
    def __sub__(self, other):
        if isinstance(other, Location):
            pos = np.setdiff1d(self.positions, other.positions)
        elif isinstance(other, int):
            pos = self.positions - other
        return Location(pos)

    def __eq__(self, other):
        if isinstance(other, Location) and len(self.positions) == len(other.positions):
            return all(self.positions == other.positions)
        return False
