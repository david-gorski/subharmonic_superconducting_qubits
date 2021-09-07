import qutip
def get_all_expectation_values(system_states, basis_states):
    expectation_values = []
    for state in basis_states:
        expectation_values.append(qutip.expect(qutip.ket2dm(state), system_states))
    return expectation_values

def expectation_value(state1, state2):
    return qutip.expect(qutip.ket2dm(state1), state2)