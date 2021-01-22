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

        hard_coded_rule = [1,2,-1,2,2,1,-2,1,2,1,-2,-2,-1,-2,1]
        steps = {1:(0,1),2:(-1,0),-1:(0,-1),-2:(1,0)}
        self._curve = self._get_curve(hard_coded_rule,steps) # compute the curve

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


    def _expand_sequence(self, seq: List[GridQubit]) -> List[GridQubit]:
        """Tries to expand given sequence with more qubits.
        Args:
            seq: Linear sequence of qubits.
        Returns:
            New continuous linear sequence which contains all the qubits from
            seq and possibly new qubits inserted in between.
        """
        i = 1
        while i < len(seq):
            path = self._find_path_between(seq[i - 1], seq[i], set(seq))
            if path:
                seq = seq[:i] + path + seq[i:]
            else:
                i += 1
        return seq

    def _find_path_between(
        self, p: GridQubit, q: GridQubit, used: Set[GridQubit]
    ) -> Optional[List[GridQubit]]:
        """Searches for continuous sequence between two qubits.
        This method runs two BFS algorithms in parallel (alternating variable s
        in each iteration); the first one starting from qubit p, and the second
        one starting from qubit q. If at some point a qubit reachable from p is
        found to be on the set of qubits already reached from q (or vice versa),
        the search is stopped and new path returned.
        Args:
            p: The first qubit, start of the sequence.
            q: The second qubit, end of the sequence.
            used: Set of forbidden qubits which cannot appear on the sequence.
        Returns:
            Continues sequence of qubits with new path between p and q, or None
            if no path was found.
        """

        def assemble_path(n: GridQubit, parent: Dict[GridQubit, GridQubit]):
            path = [n]
            while n in parent:
                n = parent[n]
                path.append(n)
            return path

        other = {p: q, q: p}
        parents = {p: dict(), q: dict()}  # type: Dict[GridQubit, Dict[GridQubit, GridQubit]]
        visited = {p: set(), q: set()}  # type: Dict[GridQubit, Set[GridQubit]]

        queue = collections.deque([(p, p), (q, q)])

        # Run two BFSs simultaneously.
        while queue:
            n, s = queue.popleft()
            for n_adj in self._c_adj[n]:
                if n_adj in visited[other[s]]:
                    # Connection has been found, construct the path and return.
                    path_s = assemble_path(n, parents[s])[-2::-1]
                    path_other = assemble_path(n_adj, parents[other[s]])[:-1]
                    path = path_s + path_other
                    if s == q:
                        path.reverse()
                    return path
                if n_adj not in used and n_adj not in visited[s]:
                    # Append n_adj to the end of queue of qubit s.
                    queue.append((n_adj, s))
                    visited[s].add(n_adj)
                    parents[s][n_adj] = n

        return None

    def _neighbors_of_excluding(self, qubit: GridQubit, used: Set[GridQubit]) -> List[GridQubit]:
        return [n for n in self._c_adj[qubit] if n not in used]

    def _get_curve(self, rule:List[int], steps: dict) -> dict:
        trace = {}
        start = min(self._c)+(5,-2) #hardcoded for now
        for r in rule:
            trace[start] = start + steps[r]
            start = trace[start]
        return trace


class _PickClosestNeighbors(ClosestSequenceSearch):
    """Minimal qubit connectivity greedy heuristic for linear sequence.
    Traverses the grid by choosing the qubit which has the least number of still
    available neighbours in each step. However, qubits with no remaining
    neighbors at all are avoided if at all possible. The idea is that this will
    cause the search to "hug the walls" and spiral inward, without falling into
    obvious traps.
    """

    def _choose_next_qubit(self, qubit: GridQubit, used: Set[GridQubit]) -> Optional[GridQubit]:

        level_2_curve = self._curve

        return level_2_curve.get(qubit)


class ClosestSequenceSearchStrategy(place_strategy.LinePlacementStrategy):
    """closest search method for linear sequence of qubits on a chip."""

    def __init__(self, start: GridQubit) -> None:
        """Initializes closest sequence search strategy.
        Args:
            algorithm: Greedy algorithm to be used. Available options are:
                best - runs all heuristics and chooses the best result,
                largest_area - on every step takes the qubit which has
                connection with the largest number of unassigned qubits, and
                minimal_connectivity - on every step takes the qubit with
                minimal number of unassigned neighbouring qubits.
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