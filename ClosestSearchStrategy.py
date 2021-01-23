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

        prod_rule = [1,2,-1,2,2,1,-2,1,2,1,-2,-2,-1,-2,1]
        steps = {1:(0,1),2:(-1,0),-1:(0,-1),-2:(1,0)}
        self._curve = self._get_curve(prod_rule,steps) # compute the curve

    def get_or_search(self) -> List[GridQubit]:
        """Starts the search or gives previously calculated sequence.
        Returns:
            The linear qubit sequence found.
        """
        if not self._sequence:
            self._sequence = self._find_sequence()
        return self._sequence

    @abc.abstractmethod
    def _choose_next_qubit(self, qubit: GridQubit, used: Set[GridQubit]) -> Optional[GridQubit]:
        """Selects next qubit on the linear sequence.
        Args:
            qubit: Last qubit which is already present on the linear sequence
                   of qubits.
            used: Set of forbidden qubits which can not be used.
        Returns:
            Next qubit to be appended to the linear sequence, chosen according
            to the greedy heuristic method. The returned qubit will be the one
            passed to the next invocation of this method. Returns None if no
            more qubits are available and search should stop.
        """

    def _find_sequence(self) -> List[GridQubit]:
        """Looks for a sequence starting at a given qubit.
        Search is issued twice from the starting qubit, so that longest possible
        sequence is found. Starting qubit might not be the first qubit on the
        returned sequence.
        Returns:
            The longest sequence found by this method.
        """
        # # Run the first pass and drop starting qubit from the found sequence.
        # tail = self._sequence_search(self._start, [])
        # tail.pop(0)

        # # Run the second pass and keep the starting qubit.
        # head = self._sequence_search(self._start, tail)
        # head.reverse()

        return self._sequence_search(self._start, [])

    def _sequence_search(self, start: GridQubit, current: List[GridQubit]) -> List[GridQubit]:
        """Search for the continuous linear sequence from the given qubit.
        This method is called twice for the same starting qubit, so that
        sequences that begin and end on this qubit are searched for.
        Args:
            start: The first qubit, where search should be triggered from.
            current: Previously found linear sequence, which qubits are
                     forbidden to use during the search.
        Returns:
            Continuous linear sequence that begins with the starting qubit and
            does not contain any qubits from the current list.
        """
        used = set(current)
        seq = []
        n = start  # type: Optional[GridQubit]
        while n is not None:
            # Append qubit n to the sequence and mark it is as visited.
            seq.append(n)
            used.add(n)
            # Advance search to the next qubit.
            n = self._choose_next_qubit(n, used)
        return seq

    def _get_curve(self, rule:List[int], steps: dict) -> dict:
        trace = {}
        start = min(self._c)+(5,-2) #hardcoded for now
        for r in rule:
            trace[start] = start + steps[r]
            start = trace[start]
        return trace


class _PickClosestNeighbors(ClosestSequenceSearch):
    """Pick Neighbor along the hilbert curve"""

    def _choose_next_qubit(self, qubit: GridQubit, used: Set[GridQubit]) -> Optional[GridQubit]:

        level_2_curve = self._curve

        return level_2_curve.get(qubit)


class ClosestSequenceSearchStrategy(place_strategy.LinePlacementStrategy):
    """closest search method for linear sequence of qubits on a chip."""

    def __init__(self, start: GridQubit) -> None:
        """Initializes closest sequence search strategy.
        Args:
            start: GridQubit to start
        """
        self.start = start

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

        sequence.append(_PickClosestNeighbors(device, self.start).get_or_search())

        # return GridQubitLineTuple.best_of(sequence[:length]), sequence[length] if length < len(sequence) else None
        return GridQubitLineTuple.best_of(sequence[:length], length)

        # return GridQubitLineTuple.best_of(sequences, length), start, step