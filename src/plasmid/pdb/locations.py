# Custom class for storing feature locations
import numpy as np
import copy

class Location:
    """
    This class holds function related to merging and storing locations
    """
    
    def __init__(self, pos:list):
        """
        Initializes the Record object
        
        Args:
            pos: list of positions
        
        Returns:
            Location object
        """
        self.positions = np.unique(pos)
        
    def copy(self):
        """
        Returns a copy of self
        """
        return copy.deepcopy(self)
        
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
        out = ", ".join(out)
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
        out = self.copy()
        x1 = list(out.positions)
        x2 = list(other.positions)
        out.positions = np.unique(x1 + x2)
        return out

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