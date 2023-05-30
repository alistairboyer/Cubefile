from typing import Iterator, Union, Tuple, Optional, Any, List, Dict
from numpy.typing import NDArray

import numpy
import os
import re

float_type = numpy.float64
"""The :attr:`dtype` for all new numpy arrays. [numpy.float64]"""

BOHR_TO_ANGSTROM: float = 5.29177210903 / 10.0
"""Conversion from Bohr to Angstrom.
`2018 CODATA <https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0>`_."""


class Cubefile:
    """
    Cubefile information.

    :param data_source: Passed to :meth:`.read` if truthy.
    """

    filename: Optional[str]
    """The file name of the data file or ``None``."""
    header: str
    """Header content from the loaded data."""
    origin: NDArray[float_type]
    """Origin of coordinate system used (Å)."""
    voxel_shape: NDArray[float_type]
    """Shape of each voxel per dimension."""
    unit_conversion: NDArray[float_type]
    """Scaling of units that has been applied relative to Å."""
    scale: NDArray[float_type]
    """Scale."""
    atoms: List[Dict[str, Union[int, float, NDArray[float_type]]]]
    """
    List of atom information.

    Atom information is stored as ``dict`` with keys:
        -   "element": atomic number [int]
        -   "charge": charge [float]
        -   "xyz": atomic coordinates (Å) [NDArray[float_type]]
    """
    voxels: NDArray[float_type]
    """Voxel data. A 3-dimensional numpy array."""

    def __init__(self, data_source: Any = None) -> None:
        self.reset()
        if data_source:
            self.read(data_source=data_source)

    def reset(self) -> None:
        """Initialize object variables."""
        self.filename = None
        self.header = ""
        self.origin = numpy.zeros((3,), dtype=float_type)
        self.voxel_shape = numpy.ones((3, 3), dtype=float_type)
        self.unit_conversion = numpy.ones((3,), dtype=float_type)
        self.scale = numpy.asarray([1.0, 1.0, 1.0], dtype=float_type)
        self.atoms = list()
        self.voxels = numpy.zeros((0,), dtype=float_type)

    @property
    def voxel_count(self) -> Tuple[int, ...]:
        """Number of voxels in each dimension. Alias of :attr:`.voxels.shape`"""
        return self.voxels.shape

    @property
    def voxel_total(self) -> int:
        """Total number of voxels. Alias of :attr:`.voxels.size`"""
        return self.voxels.size

    @property
    def atom_count(self) -> int:
        """Number of atoms."""
        return len(self.atoms)

    @property
    def max_voxel_val(self) -> float:
        """Maximum absolute voxel value."""
        return float(numpy.abs(self.voxels).max())

    def read(self, data_source: Any) -> None:
        """
        Read a cubefile using :meth:`.read_iterator`.

        :param:
            #.  If :attr:`os.path.isfile(data_source)`
                then the file is opened with :func:`open` and loaded.
                This sets the :attr:`filename` attribute.

            #.  If :attr:`data_source` is a ``str``
                then the :attr:`data_source` is loaded using :func:`splitlines`.

            #.  If :attr:`data_source` is an ``Iterator``
                then its contents are loaded.

        :raises ValueError: if the :attr:`data_source` could not be loaded.

        """
        # Path to file
        try:
            if os.path.isfile(data_source):
                with open(data_source, "rt") as iterator:
                    self.filename = data_source
                    return self.read_iterator(iterator)
        except Exception:
            pass

        # Newline delimited str
        if isinstance(data_source, str):
            return self.read_iterator(iter(data_source.splitlines()))

        # Iterator
        try:
            return self.read_iterator(iter(data_source))
        except Exception:
            pass
        # Other passed
        raise ValueError("Could not read", data_source)

    def read_iterator(self, iterator: Iterator[str]) -> None:
        """
        Read cube data from an iterator.
        See also: :meth:`.read`.

        File format reference: `http://paulbourke.net/dataformats/cube/ <http://paulbourke.net/dataformats/cube/>`_.

        :param iterator: An iterator that yields cubefile data line by line as `str`.
        :type iterator: Iterator[str]

        :raises ValueError: if the amount of voxel data is incorrect.
        :raises ValueError: for parsing errors.
        :raises ValueError: if non-square voxels are encountered.

        """
        i: int
        split_line: List[str]
        try:
            # Lines 1-2 = header
            self.header = next(iterator) + next(iterator)

            # Line  3   = atom numbers and origin
            split_line = next(iterator).split()
            atom_count: int = int(split_line[0])
            self.origin = numpy.asarray(list(map(float, split_line[1:4])))

            # Lines 4-6 = voxel count, dimensions and units
            voxel_count: List[int] = [0, 0, 0]
            for i in range(0, 3):
                split_line = next(iterator).split()
                voxel_count[i] = int(split_line[0])
                # Negative values means Angstroms, Positive means Bohr
                if voxel_count[i] < 0:
                    voxel_count[i] = -voxel_count[i]
                    self.unit_conversion[i] = 1.0  # default is Angstrom
                else:
                    self.unit_conversion[i] = BOHR_TO_ANGSTROM
                self.voxel_shape[i] = numpy.asarray(split_line[1:4], dtype=float_type)

            # Check for non-square voxels
            # Multiply the voxel shape by 1-identity matrix, and check for any noon-zero values
            if numpy.multiply(
                self.voxel_shape, 1.0 - numpy.identity(3, dtype=float_type)
            ).any():
                raise ValueError("Non square voxel shape.")

            # Convert origin now we know about units
            self.origin = numpy.multiply(self.origin, self.unit_conversion)

            # Set up voxels and dimensions
            self.voxels = numpy.zeros(voxel_count, dtype=float_type)
            self.scale = numpy.multiply(
                numpy.linalg.norm(self.voxel_shape, axis=1),
                self.unit_conversion,
            )

            # Lines 7-n_atoms+7 = atom type, charge and position
            for _ in range(atom_count):
                split_line = next(iterator).split()
                self.atoms.append(
                    {
                        "element": int(split_line[0]),
                        "charge": float(split_line[1]),
                        "xyz": numpy.multiply(
                            numpy.asarray(split_line[2:5], dtype=float_type),
                            self.unit_conversion,
                        ),
                    }
                )

            # Collect remaining data
            self.voxels = numpy.asarray(
                re.findall(r"\S+", "".join(iterator)),
                dtype=float_type,
                order="C",
            ).reshape(self.voxels.shape)
            if not self.voxels.shape == tuple(voxel_count):
                raise ValueError("Could not read the correct number of voxels")

        # If the file is incorrectly formatted or contains the wrong number of voxels
        # Then reset the object and reraise any exception as a ValueError
        except Exception as e:
            self.reset()
            raise ValueError("Error reading file", *e.args)

    def __str__(self) -> str:
        if self.filename:
            return f"Cubefile.Cubefile with {'×'.join(map(str, self.voxels.shape))} voxels, loaded from {self.filename}."
        return f"Cubefile.Cubefile with {'×'.join(map(str, self.voxels.shape))} voxels."

    def __repr__(self) -> str:
        return str(self)


if __name__ == "__main__":
    cf = Cubefile("/var/www/python/Cubefile/_testfiles/caffeine_54.cube")
    print(cf)
