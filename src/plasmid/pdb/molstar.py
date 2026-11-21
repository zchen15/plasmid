# Module for molview viewer functionality
#
# science libs
import tempfile
from typing import TYPE_CHECKING

import molviewspec as mvs

from .locations import Location

if TYPE_CHECKING:
    from .record import Record

# Mapping of feature representation names (nglview/py3Dmol style) to the
# representation types supported by molviewspec.
REPRESENTATION_TYPES = {
    "cartoon": "cartoon",
    "backbone": "backbone",
    "ball+stick": "ball_and_stick",
    "ball-and-stick": "ball_and_stick",
    "ball_and_stick": "ball_and_stick",
    "stick": "ball_and_stick",
    "licorice": "ball_and_stick",
    "spacefill": "spacefill",
    "sphere": "spacefill",
    "vdw": "spacefill",
    "line": "line",
    "lines": "line",
    "wireframe": "line",
    "surface": "surface",
    "carbohydrate": "carbohydrate",
}

# Representation used for a feature when none is specified.
DEFAULT_REPRESENTATIONS = {
    "amino_acid": "cartoon",
    "nucleic_acid": "cartoon",
    "ligand": "ball_and_stick",
    "ion": "ball_and_stick",
    "solvent": "ball_and_stick",
}


def feature_selectors(feature):
    """Builds molviewspec selectors for the residues covered by a feature.

    Args:
        feature: Row of a Record features table.

    Returns:
        List of :class:`molviewspec.ComponentExpression` covering the
        feature's residue ranges.
    """
    res_id = feature["res_id"]
    if not isinstance(res_id, Location):
        res_id = Location(res_id)

    chain_id = str(feature["chain_id"])
    return [
        mvs.ComponentExpression(
            auth_asym_id=chain_id,
            beg_auth_seq_id=int(s.start),
            end_auth_seq_id=int(s.stop) - 1,
        )
        for s in res_id.get_slices()
    ]


def feature_representation(feature):
    """Returns the molviewspec representation type for a feature."""
    representation = feature.get("representation")
    if isinstance(representation, str):
        representation = REPRESENTATION_TYPES.get(representation.strip().lower())
    else:
        # None is stored as NaN in pandas, so treat all non-strings as unset
        representation = None
    if representation is None:
        representation = DEFAULT_REPRESENTATIONS.get(
            feature.get("entity_type"), "ball_and_stick"
        )
    return representation


def format_feature_color(representation, color, selectors=None):
    """Applies a feature color to a representation, falling back to defaults."""
    if isinstance(color, (list, tuple)):
        color = color[0] if len(color) > 0 else None
    if isinstance(color, str):
        color = color.strip()

    # unrecognized colors (including representation names accidentally stored
    # as colors) are replaced with the default element/residue coloring
    valid_color = isinstance(color, str) and color != "" and color.lower() not in REPRESENTATION_TYPES
    if not valid_color:
        representation.color(custom={"molstar_use_default_coloring": True})
        return representation

    if selectors is None or selectors == "all":
        representation.color(color=color)
    else:
        representation.color(selector=selectors, color=color)
    return representation


class Viewer:
    """
    This class holds functions pertaining to drawing structures with molviewspec.

    Attributes:
        verbose: print run settings when True.
        height: height of the viewer.
        width: width of the viewer.
        background: background color.
        record: Record object being drawn.
    """
    verbose = False
    height = 500
    width = 700
    background = "black"
    record = None

    def __init__(self, rec: "Record", **kwargs):
        """Initializes the Viewer object for the record.

        Args:
            rec: Record object.
            kwargs: Settings that override the class defaults.
        """
        for key, value in kwargs.items():
            if hasattr(self, key) and value is not None:
                setattr(self, key, value)
        self.record = rec
        self.work_dir = tempfile.TemporaryDirectory(prefix="molview_")

    def format_view(self, structure):
        """
        Apply the default and feature representations to a structure.

        Args:
            structure: molviewspec structure builder.

        Returns:
            The formatted molviewspec structure builder.
        """
        # default representations
        polymer = structure.component(selector="polymer")
        polymer_cartoon = polymer.representation(
            type="cartoon"
        ).color(custom={"molstar_use_default_coloring": True})
        ligand = structure.component(selector="ligand")
        ligand_ball_and_stick = ligand.representation(
            type="ball_and_stick"
        ).color(custom={"molstar_use_default_coloring": True})

        if self.record.features is not None and len(self.record.features) > 0:
            if self.verbose:
                print(f"applying the following representations:\n{self.record.features}")
            for _, feature in self.record.features.iterrows():
                selectors = feature_selectors(feature)
                if len(selectors) == 0:
                    continue

                representation_type = feature_representation(feature)
                representation = None
                # reuse the default representations where possible so that
                # feature colors are applied to the structure in place
                if representation_type == "cartoon" and feature.get("entity_type") in (
                    "amino_acid", "nucleic_acid",
                ):
                    representation = polymer_cartoon
                    representation_selectors = selectors
                elif representation_type == "ball_and_stick" and feature.get("entity_type") in (
                    "ligand", "ion", "solvent",
                ):
                    representation = ligand_ball_and_stick
                    representation_selectors = selectors
                else:
                    representation = structure.component(
                        selector=selectors
                    ).representation(type=representation_type)
                    representation_selectors = None

                format_feature_color(
                    representation,
                    feature.get("color"),
                    selectors=representation_selectors,
                )
        return structure

    def draw(self):
        """
        Draw the record with molviewspec.

        The Mol* viewer is displayed in the active notebook.

        Returns:
            None. The viewer is displayed as a side effect.
        """
        # generate the molviewspec
        builder = mvs.create_builder()
        builder.canvas(background_color=self.background)

        # format the structure
        structure = builder.download(url="mystructure.cif").parse(format="mmcif").model_structure()
        structure = self.format_view(structure)

        # write the structure to file
        filename = f"{self.work_dir.name}/record.cif"
        self.record.save(filename)

        # render it with molstar
        with open(filename) as f:
            cif_data = f.read()
        return builder.molstar_notebook(
            data={"mystructure.cif": cif_data},
            width=self.width,
            height=self.height,
        )

    def draw_from_url(self, url: str, file_format: str):
        """
        Draw a structure from a URL with molviewspec.

        Args:
            url: URL of the structure file.
            file_format: Format of the structure file (e.g. ``"mmcif"``).

        Returns:
            The formatted molviewspec state.
        """
        # generate the molviewspec
        builder = mvs.create_builder()
        builder.canvas(background_color=self.background)

        # fetch and format the structure
        structure = builder.download(url=url).parse(format=file_format).model_structure()
        structure = self.format_view(structure)

        # render it
        viewer = builder.get_state()
        return viewer
