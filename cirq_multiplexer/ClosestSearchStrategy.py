# Revised from cirq.google.line.placement.greedy
from typing import Dict, List, Optional, Set, Tuple, TYPE_CHECKING

import abc
import collections

# from cirq.devices import GridQubit
import cirq
from cirq.google.line.placement import place_strategy
from cirq.google.line.placement.chip import chip_as_adjacency_list
from cirq.google.line.placement.sequence import GridQubitLineTuple

if TYPE_CHECKING:
    from cirq.google.line.placement.sequence import LineSequence
    import cirq.google

class NotFoundError(Exception):
    pass

class ClosestSequenceSearch:

    def __init__(self, device: 'cirq.google.XmonDevice', start: cirq.GridQubit) -> None:
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
        self._sequence = None  # type: Optional[List[cirq.GridQubit]]

    def get_or_search(self, curve: dict, err_qubits:List[cirq.GridQubit], connectivity:List[Tuple[cirq.GridQubit,...]]) -> List[cirq.GridQubit]:
        """Starts the search or gives previously calculated sequence.
        Returns:
            The linear qubit sequence found.
        """
        if not self._sequence:
            self._sequence = self._find_sequence(curve, err_qubits, connectivity)
        return self._sequence

    @abc.abstractmethod
    def _choose_next_qubit(self, qubit: cirq.GridQubit, curve: dict, err_qubits:List[cirq.GridQubit]) -> Optional[cirq.GridQubit]:
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

    def _find_sequence(self, curve: dict, err_qubits:List[cirq.GridQubit], connectivity:List[Tuple[cirq.GridQubit,...]]) -> List[cirq.GridQubit]:
        """Looks for a sequence starting at a given qubit.
        Returns:
            The sequence found by this method.
        """

        return self._sequence_search(self._start, curve, err_qubits, connectivity)

    def _sequence_search(self, start: cirq.GridQubit, curve: dict, err_qubits:List[cirq.GridQubit], connectivity:List[Tuple[cirq.GridQubit,...]]) -> List[cirq.GridQubit]:
        """Search for the continuous linear sequence from the given qubit.
        Args:
            start: The first qubit, where search should be triggered from.
            curve: Previously found level 2 hilbert sequence, which qubits are
                     picked from along the curve during the search.
        Returns:
            Continuous linear sequence that begins with the starting qubit.
        """
        seq = []
        n = start  # type: Optional[cirq.GridQubit]
        curr_neighbors = []
        while n is not None:
            # print('curr_qubit:', n)
            # print('curr_seq:', seq)
            # print('err_qubits:', err_qubits)
            if n in err_qubits:
                # print('found error qubit:', n)
                n = self._choose_next_qubit(n, curve, err_qubits)
                continue
            # Append first qubit n to the sequence
            if not seq:
                seq.append(n)
            else:
            # Append qubit n to the sequence if connectivity still holds
                curr_neighbors = [self._neighbors(q, connectivity) for q in seq]
                curr_neighbors = sum(curr_neighbors, [])
                if n in curr_neighbors:
                    seq.append(n)
            # Advance search to the next qubit.
            n = self._choose_next_qubit(n, curve, err_qubits)
        return seq

    def _neighbors(self, qubit:cirq.GridQubit, connectivity:List[Tuple[cirq.GridQubit,...]]):
        possibles = [
           cirq.GridQubit(qubit.row+1, qubit.col),
           cirq.GridQubit(qubit.row, qubit.col+1),
           cirq.GridQubit(qubit.row-1, qubit.col),
           cirq.GridQubit(qubit.row, qubit.col-1),
        ]

        return [q for q in possibles if (q in self._c and self._check_connectivity(qubit, q, connectivity))]

    def _check_connectivity(self, qubit1:cirq.GridQubit, qubit2:cirq.GridQubit, connectivity:List[Tuple[cirq.GridQubit,...]]):
        return True if (qubit1, qubit2) or (qubit2,qubit1) in connectivity else False



class _PickClosestNeighbors(ClosestSequenceSearch):
    """Pick Next Qubit along the hilbert curve"""

    def _choose_next_qubit(self, qubit: cirq.GridQubit, curve: dict, err_qubits:List[cirq.GridQubit]) -> Optional[cirq.GridQubit]:
        return curve.get(qubit)


class ClosestSequenceSearchStrategy(place_strategy.LinePlacementStrategy):
    """closest search method for linear sequence of qubits on a chip."""

    def __init__(self, start: cirq.GridQubit, curve: dict, err_qubits:List[cirq.GridQubit], connectivity:List[Tuple[cirq.GridQubit,...]]) -> None:
        """Initializes closest sequence search strategy.
        Args:
            start: cirq.GridQubit to start
        """
        self.start = start
        self.curve = curve
        self.err_qubits = err_qubits
        self.connectivity = connectivity

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

        sequence.append(_PickClosestNeighbors(device, self.start).get_or_search(self.curve, self.err_qubits, self.connectivity))

        # return GridQubitLineTuple.best_of(sequence[:length]), sequence[length] if length < len(sequence) else None
        return GridQubitLineTuple.best_of(sequence[:length], length)

        # return GridQubitLineTuple.best_of(sequences, length), start, step