# Module for storing feature locations as arrays of [start, stop] pairs.

import numpy as np


class Locations:
    """
    This class holds function related to merging and storing locations.

    Unlike the atom-index based :class:`plasmid.pdb.locations.Location`, this
    version stores ``self.positions`` as an array of ``[start, stop]`` pairs,
    one pair per contiguous range of positions.
    """

    positions = np.array([], dtype=int)

    def __init__(self, pos: list | np.ndarray):
        """
        Initializes the Locations object.

        Args:
            pos: list or array of positions.

        Returns:
            Locations: Locations object.

        Raises:
            ValueError: if ``pos`` is empty.
        """

        arr = np.unique(np.asarray(pos, dtype=int))
        if len(arr) == 0:
            raise ValueError("Locations requires at least one position")
        if len(arr) == 1:
            self.positions = np.array([[int(arr[0]), int(arr[0]) + 1]], dtype=int)
            return
        # Find contiguous run boundaries:
        # a break occurs when next_val != current_val + 1
        # i.e., arr[i+1] != arr[i] + 1  ->  arr[1:] != arr[:-1] + 1
        c1 = np.concatenate(([True], arr[1:] != arr[:-1] + 1))  # start markers
        c2 = np.concatenate((arr[1:] != arr[:-1] + 1, [True]))  # end markers
        starts = arr[c1]
        ends = arr[c2]
        self.positions = np.column_stack((starts, ends)).astype(int)

    def get_positions(self):
        """
        Flatten the ``[start, stop]`` pairs into a flat array of positions.

        Uses vectorized numpy to avoid slow Python loops.

        Returns:
            np.ndarray: flat array of positions covered by this location.
        """

        return np.concatenate([np.arange(s, e) for s, e in self.positions])

    def get_slices(self):
        """
        Transform list of indices into list of slices

        Returns:
            list of slices
        """

        return [slice(int(s), int(e)) for s, e in self.positions]

    def remap(self, pos: dict):
        """
        Remap atom_id

        Args:
            pos: dict mapping old atom ids to new atom ids.

        Returns:
            Locations: New Locations object with remapped positions.
        """

        positions = [pos[int(i)] for i in self.get_positions()]
        return Locations(positions)

    def contains(self, pos: list):
        """
        check if location contains all positions

        Args:
            pos: list of positions to check.

        Returns:
            bool: True if all positions are contained in this location.
        """

        return bool(np.all(np.isin(self.get_positions(), np.asarray(pos, dtype=int))))

    def intersection(self, pos: list):
        """
        keep only positions that intersect

        Args:
            pos: list of positions to intersect with.

        Returns:
            Locations: New Locations object containing only the intersecting positions.
        """

        positions = np.intersect1d(self.get_positions(), np.asarray(pos, dtype=int))
        return Locations(positions)

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
            other: Locations object or integer offset.

        Returns:
            Locations: New Locations object resulting from the addition.
        """

        if isinstance(other, Locations):
            flat = np.concatenate((self.get_positions(), other.get_positions()))
        elif isinstance(other, int):
            flat = self.get_positions() + other
        return Locations(flat)

    def __radd__(self, other):
        """
        Performs other + self
        """

        return self.__add__(self)

    def __iadd__(self, other):
        """
        Performs in-place self + other
        """

        if isinstance(other, Locations):
            self.positions = Locations(
                np.concatenate((self.get_positions(), other.get_positions()))
            ).positions
        elif isinstance(other, int):
            self.positions = self.positions + other
        return self

    def __sub__(self, other):
        """
        Performs self - other
        """

        if isinstance(other, Locations):
            flat = np.setdiff1d(self.get_positions(), other.get_positions())
        elif isinstance(other, int):
            flat = self.get_positions() - other
        return Locations(flat)

    def __isub__(self, other):
        """
        Performs in-place self - other
        """

        if isinstance(other, Locations):
            self.positions = Locations(
                np.setdiff1d(self.get_positions(), other.get_positions())
            ).positions
        elif isinstance(other, int):
            self.positions = self.positions - other
        return self

    def __eq__(self, other):
        """
        Performs self == other
        """

        if isinstance(other, Locations) and len(self.positions) == len(other.positions):
            return np.array_equal(self.positions, other.positions)
        return False