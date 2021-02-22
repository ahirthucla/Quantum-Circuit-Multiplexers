import cirq
from cirq_multiplexer.ClosestSearchStrategy import ClosestSequenceSearchStrategy, NotFoundError
from functools import reduce
import networkx as nx
import time
import sys
from cirq.contrib.routing import route_circuit

from typing import Union, Sequence, Iterable, Callable, Optional, List, Set
from numbers import Number

DEBUG = False

def get_curve() -> dict:
    """Looks for a hilbert curve sequence starting at a given qubit.
    Returns:
        The dictonary of qubits that map to each other along the curve found by this method, 
        which is similar to a LinkedList of qubits.
    """
    hard_coded_rule = [1,2,-1,2,2,1,-2,1,2,1,-2,-2,-1,-2,1]
    steps = {1:(0,1),2:(-1,0),-1:(0,-1),-2:(1,0)}
    trace = {}
    start = cirq.GridQubit(5,3) #hardcoded for now
    for r in hard_coded_rule:
        trace[start] = start + steps[r]
        start = trace[start]
    return trace

def get_error_qubits(project_id, processor_id, threshold):

    # query for the latest calibration
    engine = cirq.google.Engine(project_id=project_id)
    processor = engine.get_processor(processor_id=processor_id)
    latest_calibration = processor.get_current_calibration()

    err_qubits = []
    for metric_name in latest_calibration:
        for qubit_or_pair in latest_calibration[metric_name]:
            metric_value = latest_calibration[metric_name][qubit_or_pair]
            # find the qubits that have higher error probability(above the threshold)
            if metric_value[0] > threshold:
                # get the first qubit of the tuple from a metric key
                err_qubits.append(qubit_or_pair[0])
    return err_qubits

def mult_qubit_opcount_cost(circuit:'cirq.Circuit'): 
    """ Returns tuple of the number of operations applied to more than one qubit, and the total number of operations """
    cnot_count, all_count = reduce(lambda xy,op:(xy[0] + (len(op.qubits) > 1), xy[1] + 1), circuit.all_operations(), (0,0))
    return cnot_count, all_count

def device_connectivity(device:'cirq.Device', limit: Set['cirq.Qid']) -> nx.Graph:
    # Below is how it seems we "should" extract device connectivity based on the functionality available to the device class

    #pairs = {}
    #for pair in itertools.combinations(device.qubit_set(), 2):
    #    op = cirq.CX.on(*pair)
    #    op = device.decompose_operation(op)
    #    try:
    #        device.validate_operation(op)
    #    except ValueError as e:
    #        print(e)
    #        continue
    #    pairs.add(pair)
    #return nx.Graph(pairs)

    # Below is how connectivity is extracted in device's string representation, when being printed
    #   except it's slightly modified to use Sycamore's connectivity, for the qubits available to the device
    #   this is because the google engine device doesn't seem to have any connectivity?

    pairs = {
                pair
    #            for gate_defs in device.gate_definitions.values()
                for gate_defs in cirq.google.Sycamore.gate_definitions.values()
                for gate_def in gate_defs if gate_def.number_of_qubits == 2
    #            for pair in gate_def.target_set if len(pair) == 2
                for pair in gate_def.target_set if len(pair) == 2 and set(pair).issubset(device.qubit_set()) and set(pair).issubset(limit)
            }
    return nx.Graph(pairs)

def naive_line_mapping(circuit: 'cirq.Circuit', 
                        device: 'cirq.google.XmonDevice', 
                        exclude : set = None,
                        context = None):
    """ Wraps cirq.google.line_on_device to support excluding a set of qubits from the search. 
     Returns cirq.devices.LineQubit -> cirq.devices.GridQubit mapping function
    """

    if exclude is None:
        exclude = set()

    if context == None:
        context = cirq.GridQubit(5,3)

    width = len(circuit.all_qubits())

    curve = get_curve()

    method = ClosestSequenceSearchStrategy(context, curve, exclude)

    try:
        line = cirq.google.line_on_device(device, length=width, 
            method=method)
    except cirq.google.line.placement.sequence.NotFoundError as e:
        raise NotFoundError('No line placment found.')

    qubit_map = dict(zip(sorted(circuit.all_qubits()), line))
    qubit_map = dict(zip(circuit.all_qubits(), line))
    return qubit_map, curve[line[-1]]


def multiplex_onto_sycamore(circuits:Iterable['cirq.Circuit'],
        device: 'cirq.google.XmonDevice',
        mapping_function: Callable = naive_line_mapping,
        exclude_always: Set['cirq.google.GridQubit'] = None):
    """ Combines circuits by placing them on different qubits, according to those available to a specific device
    Args: 
        circuits: ordered iterable of circuits to multiplex. circuits are placed first-come,first-servr
        device: google device to multiplex onto
        mapping_function: fun(circuit, device, exclude) -> qubit_map: function that produces mapping from qubits in circuit to qubits on device, while avoiding device qubits in exclude set
        exclude_always: set of device qubits that are always avoided (for example those with low calibration data)
    Returns: 
        generator yielding (cumulative_circuit, circuits_included) pairs
            cumulative circuit is the multiplexed circuit
            circuits_included are the indices of the circuits included in the multiplexed circuit. used to index the measurement keys
    """

    if exclude_always is None:
        exclude_always = set()

    # set of device qubits to exclude from the qubit map
    exclude = set(exclude_always)

    # empty circuit to build upon 
    cumulative_circuit = cirq.Circuit()

    # indices for circuits included in cumulative one
    circuits_included = set()

    # context/meta-information for mapping function
    #   used by mapping function to communicate between different calls
    context = None

    for index, circuit in enumerate(circuits):
        # generate a qubit map according to the mapping_function
        try:
            qubit_map, context = mapping_function(circuit, device, exclude, context)
        except NotFoundError:
            # check if circuit fits with only base qubits excluded
            exclude = set(exclude_always)
            context = None
            qubit_map, context = mapping_function(circuit, device, exclude, context)

            # yield old cumulative circuit and start a new one
            yield cumulative_circuit, circuits_included
            cumulative_circuit = cirq.Circuit()
            circuits_included = set()


        # update measurement keys with the circuit index
        measurement_key_map = {key:str(index)+"."+key for key in circuit.all_measurement_keys()}
        if DEBUG:
            print('measurement_key_map:',measurement_key_map)
            pre = len(list(circuit.all_operations()))
            pre_time = time.process_time()
        circuit = cirq.with_measurement_key_mapping(circuit, measurement_key_map)
        if DEBUG:
            print("measurement_1: ", len(list(circuit.all_operations())) - pre, time.process_time()-pre_time)

        try:
            # compile and map circuit according to qubit_map
            if DEBUG:
                pre = len(list(circuit.all_operations()))
                pre_time = time.process_time()
            circuit = cirq.google.optimized_for_sycamore(circuit=circuit, new_device=device, optimizer_type='sycamore')
            if DEBUG:
                print("optimize: ", len(list(circuit.all_operations())) - pre, time.process_time()-pre_time)
        except ValueError as e:

            # split multiple qubit measurement operations into multiple, single qubit measurement operations
            # a workaround to use ccr
            def split_measure(measure_gate:'cirq.GateOperation') -> 'cirq.GateOperation':
                if not cirq.protocols.is_measurement(measure_gate):
                    yield measure_gate
                    return
                key = cirq.protocols.measurement_key(measure_gate)
                yield cirq.Moment([cirq.measure(qubit, key=key+'.'+str(qubit)) for qubit in measure_gate.qubits])

            if DEBUG:
                pre = len(list(circuit.all_operations()))
                pre_time = time.process_time()
                print('all_measurement_keys', circuit.all_measurement_keys())
            circuit = cirq.Circuit(*map(split_measure, circuit.all_operations()))
            if DEBUG:
                print('all_measurement_keys2', circuit.all_measurement_keys())
                print("measurement_2: ", len(list(circuit.all_operations())) - pre, time.process_time()-pre_time)

            # apply swaps while remapping circuit
            graph = device_connectivity(device, limit=qubit_map.values())
            reverse_map = {p:l for l,p in qubit_map.items()}
            if DEBUG:
                pre = len(list(circuit.all_operations()))
                pre_time = time.process_time()
                print(circuit)
            swap_network = route_circuit(circuit=circuit, device_graph=graph, algo_name='greedy', initial_mapping=reverse_map)
            circuit = swap_network.circuit
            if DEBUG:
                print(circuit)
                print([op.gate for op in circuit.all_operations()])
                print("route_and_map: ", len(list(circuit.all_operations())) - pre, time.process_time()-pre_time)

            # compile and map circuit according to qubit_map AGAIN
            if DEBUG:
                pre = len(list(circuit.all_operations()))
                pre_time = time.process_time()
            circuit = cirq.google.optimized_for_sycamore(circuit=circuit, new_device=device, optimizer_type='sycamore')
            if DEBUG:
                print("optimize: ", len(list(circuit.all_operations())) - pre, time.process_time()-pre_time)

        # validate circuit
        device.validate_circuit(circuit)

        # update set of excluded qubits
        used_qubits = circuit.all_qubits()
        exclude.update(used_qubits)

        # add to cumulative circuit and record index
        cumulative_circuit += circuit
        circuits_included.add(index)

    yield cumulative_circuit, circuits_included
    #return
    
# =========================================== For Testing ===================================================================

if __name__ == '__main__': 
    print('starting')
    # super simple testing of multiple copies of the same circuit
    def gen():
        n = 6
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
        return circuit

    # Initialize Simulator
    s = cirq.Simulator()

    cs = sorted([gen2(),gen()], key=mult_qubit_opcount_cost)
    for c in cs:
        print(c)
    print(cirq.google.Sycamore)
    err_qubits = get_error_qubits(sys.argv[1], sys.argv[2], 25)
    print(err_qubits)
    res = multiplex_onto_sycamore(circuits=cs, device=cirq.google.Sycamore, exclude_always=err_qubits)
    for c, cis in res:
        print('test')
        print(cis)
        print(c)
        print('num_q:',len(c.all_qubits()))
        print('num_op:',len(list(c.all_operations())))
        print('Simulate the circuit:')
        results=s.simulate(c)
        print(results)
        print()

        # print('qsim results:')
        # qsim_simulator = qsimcirq.QSimSimulator()
        # qsim_results = qsim_simulator.run(circuit, repetitions=5)
        # print(qsim_results)
        # print()
