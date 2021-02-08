# Revised from cirq.google.line.placement.greedy
from typing import Dict, List, Optional, Set, TYPE_CHECKING

import abc
import collections

from cirq.devices import GridQubit
from cirq.google.line.placement import place_strategy
from cirq.google.line.placement.chip import chip_as_adjacency_list
from cirq.google.line.placement.sequence import GridQubitLineTuple

if TYPE_CHECKING:
    from cirq.google.line.placement.sequence import LineSequence
    import cirq.google

class NotFoundError(Exception):
    pass

class ClosestSequenceSearch:

    def __init__(self, device: 'cirq.google.XmonDevice', start: GridQubit) -> None:
        """Level 2 Hilbert curve search constructor.
        Args:
            device: Chip description.
            start: Starting qubit.
        Raises:
            ValueError: When start qubit is not part of a chip.
        """
        if start not in device.qubits:
            raise ValueError('Starting qubit must be a qubit on the chip')

        self._c = device.qubits
        self._c_adj = chip_as_adjacency_list(device)
        self._start = start
        self._sequence = None  # type: Optional[List[GridQubit]]

    def get_or_search(self, curve: dict, err_qubits:List[GridQubit]) -> List[GridQubit]:
        """Starts the search or gives previously calculated sequence.
        Returns:
            The linear qubit sequence found.
        """
        if not self._sequence:
            self._sequence = self._find_sequence(curve, err_qubits)
        return self._sequence

    @abc.abstractmethod
    def _choose_next_qubit(self, qubit: GridQubit, curve: dict, err_qubits:List[GridQubit]) -> Optional[GridQubit]:
        """Selects next qubit on the linear sequence.
        Args:
            qubit: Last qubit which is already present on the linear sequence
                   of qubits.
            curve: dictionary of qubits that map to each other along the curve.
        Returns:
            Next qubit to be appended to the linear sequence, chosen according
            to the hilbert curve. The returned qubit will be the one
            passed to the next invocation of this method. Returns None if no
            more qubits are available and search should stop.
        """

    def _find_sequence(self, curve: dict, err_qubits:List[GridQubit]) -> List[GridQubit]:
        """Looks for a sequence starting at a given qubit.
        Returns:
            The sequence found by this method.
        """

        return self._sequence_search(self._start, curve, err_qubits)

    def _sequence_search(self, start: GridQubit, curve: dict, err_qubits:List[GridQubit]) -> List[GridQubit]:
        """Search for the continuous linear sequence from the given qubit.
        Args:
            start: The first qubit, where search should be triggered from.
            curve: Previously found level 2 hilbert sequence, which qubits are
                     picked from along the curve during the search.
        Returns:
            Continuous linear sequence that begins with the starting qubit.
        """
        seq = []
        n = start  # type: Optional[GridQubit]
        while n is not None:
            # Append qubit n to the sequence
            if n not in err_qubits:
                seq.append(n)
            else:
                seq.clear()
            # Advance search to the next qubit.
            n = self._choose_next_qubit(n, curve)
        return seq


class _PickClosestNeighbors(ClosestSequenceSearch):
    """Pick Next Qubit along the hilbert curve"""

    def _choose_next_qubit(self, qubit: GridQubit, curve: dict) -> Optional[GridQubit]:
        return curve.get(qubit)


class ClosestSequenceSearchStrategy(place_strategy.LinePlacementStrategy):
    """closest search method for linear sequence of qubits on a chip."""

    def __init__(self, start: GridQubit, curve: dict, err_qubits:List[GridQubit]) -> None:
        """Initializes closest sequence search strategy.
        Args:
            start: GridQubit to start
        """
        self.start = start
        self.curve = curve
        self.err_qubits = err_qubits

    def place_line(self, device: 'cirq.google.XmonDevice', length: int) -> GridQubitLineTuple:
        """Runs line sequence search.
        Args:
            device: Chip description.
            length: Required line length.
        Returns:
            Linear sequences found on the chip.
        Raises:
            ValueError: If search algorithm passed on initialization is not
                        recognized.
        """
        if not device.qubits:
            return GridQubitLineTuple()

        if self.start is None:
            raise NotFoundError("No qubit to start")

        sequence = []  # type: LineSequence

        sequence.append(_PickClosestNeighbors(device, self.start).get_or_search(self.curve, self.err_qubits))

        # return GridQubitLineTuple.best_of(sequence[:length]), sequence[length] if length < len(sequence) else None
        return GridQubitLineTuple.best_of(sequence[:length], length)

        # return GridQubitLineTuple.best_of(sequences, length), start, step