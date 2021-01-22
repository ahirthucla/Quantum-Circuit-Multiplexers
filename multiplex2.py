import cirq
from cirq import devices, Circuit
from cirq.devices import GridQubit
from cirq.google.line.placement.sequence import GridQubitLineTuple
from ClosestSearchStrategy import ClosestSequenceSearchStrategy, NotFoundError
from functools import reduce

from typing import Union, Sequence, Iterable, Callable, Optional, List
from numbers import Number

# import qsimcirq

def _get_curve() -> dict:
    hard_coded_rule = [1,2,-1,2,2,1,-2,1,2,1,-2,-2,-1,-2,1]
    steps = {1:(0,1),2:(-1,0),-1:(0,-1),-2:(1,0)}
    trace = {}
    start = GridQubit(5,3) #hardcoded for now
    for r in hard_coded_rule:
        trace[start] = start + steps[r]
        start = trace[start]
    return trace

def mult_qubit_opcount_cost(circuit:'cirq.Circuit'): 
    """ Returns tuple of the number of operations applied to more than one qubit, and the total number of operations """

    cnot_count, all_count = reduce(lambda xy,op:(xy[0] + (len(op.qubits) > 1), xy[1] + 1), circuit.all_operations(), (0,0))
    return cnot_count, all_count

def naive_line_mapping(circuit: 'cirq.Circuit', 
                        device: 'cirq.google.XmonDevice', 
                        start : GridQubit,
                        exclude: Optional[Iterable[devices.LineQubit]] = None):
    """ Wraps cirq.google.line_on_device to support excluding a set of qubits from the search. 
     Returns cirq.devices.LineQubit -> cirq.devices.GridQubit mapping function
    """

    width = len(circuit.all_qubits())

    level_2_curve = _get_curve()

    if exclude: 
        # Absolutely revolting hack to search over a subset of the device's qubits
        # Ideally exclude will be an argument to the search function (in this case line_on_device)
        restore = device.qubits
        # device.qubits = device.qubits.difference(exclude)
        device.qubits = [q for q in device.qubits if q not in exclude]
    
    try:
        line = cirq.google.line_on_device(device, length=width, method=ClosestSequenceSearchStrategy(start))
        print(line)
        # print()
    except cirq.google.line.placement.sequence.NotFoundError as e:
        if exclude: device.qubits = restore
        raise e
        # should have generic sequence not found exception here (not line specific)
    
    if exclude: 
        # Fix the device and hope nothing is broken :)
        device.qubits = restore

    return lambda q: line[q.x], level_2_curve[line[-1]]


def multiplex_onto_xmon(circuits:Iterable['cirq.Circuit'],
        device: 'cirq.google.XmonDevice',
        cost_fun: Callable[['cirq.Circuit'], Union[Number, Sequence[Number]]] = mult_qubit_opcount_cost,
        mapping_function: Callable[['cirq.Circuit', 'cirq.google.XmonDevice', Optional[Iterable['cirq.Qid']]],Callable[[devices.LineQubit], devices.GridQubit]] = naive_line_mapping):
    """ Combines circuits by placing them on different qubits, according to those available to a specific device """

    # Order circuits by a cost function
    circuits = list(enumerate(sorted(circuits, key=cost_fun, reverse=True)))

    start = GridQubit(5,3)

    while circuits: 
        # Start a new running circuit
        circuit_id, circuit = circuits.pop(0)
        qubit_map, start = mapping_function(circuit, device, start)
        exclude = {qubit_map(qubit) for qubit in circuit.all_qubits()}

        circuit = cirq.google.optimized_for_sycamore(circuit=circuit, new_device=device, qubit_map=qubit_map)

        # Update measurement keys
        measurement_key_map = {key:str(circuit_id)+"_"+key for key in circuit.all_measurement_keys()}
        new_circuit = cirq.with_measurement_key_mapping(circuit, measurement_key_map)

        remaining_circuits = []
        for index, (circuit_id, circuit) in enumerate(circuits):
            # Try to find a place on the remaining qubits of the device to place additional circuit
            try: 
                qubit_map, start = mapping_function(circuit, device, start, exclude=exclude)
            # except cirq.google.line.placement.sequence.NotFoundError:
            except NotFoundError:
                remaining_circuits.append((circuit_id, circuit))
                continue

            new_qubits = {qubit_map(qubit) for qubit in circuit.all_qubits()}
            assert exclude.isdisjoint(new_qubits), "mapping function does not properly exclude qubits"
            exclude |= new_qubits

            # Place the additional circuit on found qubits
            circuit = cirq.google.optimized_for_sycamore(circuit=circuit, new_device=device, qubit_map=qubit_map)
            # Fix it's measurement keys
            measurement_key_map = {key:str(circuit_id)+"_"+key for key in circuit.all_measurement_keys()}
            circuit = cirq.with_measurement_key_mapping(circuit, measurement_key_map)
            # Add to the running circuit
            new_circuit = new_circuit + circuit

            assert new_circuit.all_qubits() == exclude, "new circuit does not have expected qubits"

        circuits = remaining_circuits


        # Ensure that the new circuit is valid for the device and yield it
        device.validate_circuit(new_circuit)
        yield new_circuit

if __name__ == '__main__': 
    # super simple testing of multiple copies of the same circuit
    def gen():
        n = 10
        depth = 2
        circuit = cirq.Circuit(
            cirq.CZ(cirq.LineQubit(j), cirq.LineQubit(j + 1))
            for i in range(depth)
            for j in range(i % 2, n - 1, 2)
        )
        circuit.append(cirq.measure(*cirq.LineQubit.range(n), key='all'))
        #circuit.append(cirq.measure(*cirq.LineQubit.range(n)))
        return circuit

    def gen2():
        n = 5
        depth = 2
        circuit = cirq.Circuit(
            cirq.H(cirq.LineQubit(i))
            for i in range(n)
        )
        circuit.append(cirq.X(cirq.LineQubit(i)) for i in range(n))
        circuit.append(cirq.H(cirq.LineQubit(i)) for i in range(n))
        circuit.append(cirq.CX(cirq.LineQubit(0),cirq.LineQubit(1)))
        circuit.append(cirq.measure(*cirq.LineQubit.range(n), key='all'))
        #circuit.append(cirq.measure(*cirq.LineQubit.range(n)))
        return circuit

    # Initialize Simulator
    s = cirq.Simulator()

    cs = [gen2() for _ in range(2)]
    print(cs[0])
    res = multiplex_onto_xmon(cs, cirq.google.Sycamore)
    print(cirq.google.Sycamore)
    for c in res:
        print(c)
        print(len(c.all_qubits()))
        print('Simulate the circuit:')
        results=s.simulate(c)
        print(results)
        print()

        # print('qsim results:')
        # qsim_simulator = qsimcirq.QSimSimulator()
        # qsim_results = qsim_simulator.run(circuit, repetitions=5)
        # print(qsim_results)
        # print()
