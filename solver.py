from pulp import LpProblem, LpVariable, LpMaximize, value, lpSum, lpDot, PULP_CBC_CMD
from enum import Enum
from math import sqrt

from patient_donor_pairs import generate_patient_donor_pair

# generate patient-donor pairs
number_of_pairs = 500
all_pairs = []
for i in range(number_of_pairs):
    all_pairs.append(generate_patient_donor_pair())

# problem type enum
class ProblemType(Enum):
    SIMPLE = 1
    POTENTIALS = 2
    FAIRNESS = 3

# cycle data structure
class Cycle:
    def __init__(self, pairs):
        self.pairs = pairs # indices, not pair objects
        self.size = len(self.pairs)

# chain data structure
class Chain:
    def __init__(self, altruistic_donor, pairs):
        self.altruistic_donor = altruistic_donor # index, not donor object
        self.pairs = pairs # indices, not pair object
        self.size = len(self.pairs)

# graph data structure
class Graph:
    def __init__(self, pairs, altruistic_donors, problem_type, curr_time):
        self.pairs = pairs
        self.edges = Graph.find_edges(self.pairs)
        self.cycles = Graph.find_cycles(self.pairs, self.edges)
        self.problem_type = problem_type
        self.curr_time = curr_time
        self.cycle_weights = Graph.find_cycle_weights(self.problem_type, self.pairs, self.cycles, self.curr_time)
        self.altruistic_donors = altruistic_donors
        self.chains = Graph.find_chains(self.altruistic_donors, self.pairs, self.edges)
        self.chain_weights = Graph.find_chain_weights(self.problem_type, self.pairs, self.altruistic_donors, self.chains, self.curr_time)

    # create adjacency list representation of graph of pairs
    def find_edges(pairs):
        edges = [set() for _ in range(len(pairs))]

        for i in range(len(pairs)):
            for j in range(i+1, len(pairs)):
                if pairs[i].donor.is_compatible_with_patient(pairs[j].patient):
                    edges[i].add(j)
                if pairs[j].donor.is_compatible_with_patient(pairs[i].patient):
                    edges[j].add(i)
        
        return edges

    # find all 2 and 3 cycles for graph of pairs
    def find_cycles(pairs, edges):
        cycles = set()

        for i in range(len(pairs)):
            for j in range(i+1, len(pairs)):
                # checks for 2-cycles
                if j in edges[i] and i in edges[j]: 
                    cycles.add((i, j))

                # check for 3-cycles in one direction
                if j in edges[i]:
                    for k in range(len(pairs)):
                        if k in edges[j] and i in edges[k]:
                            c = [i, j, k]
                            ind = c.index(min(c))
                            c = c[ind:] + c[:ind]
                            cycles.add(tuple(c))

                # check for 3-cycles in other direction
                if i in edges[j]:
                    for k in range(len(pairs)):
                        if k in edges[i] and j in edges[k]:
                            c = [j, i, k]
                            ind = c.index(min(c))
                            c = c[ind:] + c[:ind]
                            cycles.add(tuple(c))

        cycles = list(cycles)
        cycles = [Cycle(list(c)) for c in cycles]

        # print(f'Number of Cycles in Graph: {len(cycles)}')
        return cycles

    # function that establishes the optimization weights for each cycle based on the problem type
    def find_cycle_weights(problem_type, pairs, cycles, curr_time):
        if problem_type == ProblemType.SIMPLE: # if simple, weights are size of the cycle
            return [c.size for c in cycles]
        elif problem_type == ProblemType.POTENTIALS: # if potentials, weights are size of cycle minus potential of each vertex in cycle
            return [c.size - sum([pairs[p].patient.potential + pairs[p].donor.potential for p in c.pairs]) for c in cycles]
        elif problem_type == ProblemType.FAIRNESS: # if fairness, weights take into account waiting time and time before departure
            return [1 + sum([sqrt(curr_time - pairs[p].arrival_time) + max(0, 10 - (pairs[p].departure_time - curr_time)) for p in c.pairs]) for c in cycles]

    # function that finds all the chains in a graph from a given list of altruistic donors
    def find_chains(altruistic_donors, pairs, edges):
        chains = []

        # loop over all altruistic donors
        for d in range(len(altruistic_donors)):
            donor = altruistic_donors[d]

            # get all elements that could start a chain
            first_elems = [i for i in range(len(pairs)) if donor.is_compatible_with_patient(pairs[i].patient)]

            # function to find all chains of size at most 10
            def get_chains(start, past, found):
                found.add(start)
                chains.append(Chain(d, past + [start])) # add current chain to list

                # stop if size exceeds 10
                if len(past) + 1 == 10:
                    return

                # explore all possible next pairs that haven't been searched yet
                next_pairs = [i for i in edges[start] if not i in found]
                for next_pair in next_pairs:
                    get_chains(next_pair, past + [start], found)

            # search for all possible chains
            for first_elem in first_elems:
                get_chains(first_elem, [], set())

        # print(f'Number of Chains in Graph: {len(chains)}')
        return chains

    def find_chain_weights(problem_type, pairs, altruistic_donors, chains, curr_time):
        if problem_type == ProblemType.SIMPLE: # if simple, weights are size of the cycle
            return [c.size for c in chains]
        elif problem_type == ProblemType.POTENTIALS: # if potentials, weights are size of chain minus potential of each vertex in cycle and minus potential of donor * constant
            return [c.size - sum([pairs[p].patient.potential + pairs[p].donor.potential for p in c.pairs]) - 3*altruistic_donors[c.altruistic_donor].potential for c in chains]
        elif problem_type == ProblemType.FAIRNESS: # if fairness, weights take into account waiting time and time before departure
            return [1 + sum([sqrt(curr_time - pairs[p].arrival_time) + max(0, 10 - (pairs[p].departure_time - curr_time)) for p in c.pairs]) for c in chains]

def solve_kidney_matching(pairs, altruistic_donors, problem_type, curr_time):
    # construct graph
    graph = Graph(pairs, altruistic_donors, problem_type, curr_time)

    # get cycles and chains
    cycles = graph.cycles
    chains = graph.chains

    # initialize problem
    problem = LpProblem('kidney_matching', LpMaximize)

    # create decision variables for each cycle and chain
    cycle_vars = LpVariable.dicts('cycle', range(len(cycles)), cat='Binary')
    chain_vars = LpVariable.dicts('chain', range(len(chains)), cat='Binary')

    cycles_chains_with_vertex = [[] for _ in range(len(pairs))] # cycles_with_vertex[i] is list of cycle and chain variables that contain vertex i
    chains_with_donor = [[] for _ in range(len(altruistic_donors))] # chains_with_donor[i] is list of chain variables that start with donor i

    # update cycles_chains_with_vertex by looping over all cycles
    for i in range(len(cycles)):
        for j in cycles[i].pairs:
            cycles_chains_with_vertex[j].append(cycle_vars[i])

    # update cycles_chains_with_vertex and chains_with_donor by looping over all chains
    for i in range(len(chains)):
        chains_with_donor[chains[i].altruistic_donor].append(chain_vars[i])
        for j in chains[i].pairs:
            cycles_chains_with_vertex[j].append(chain_vars[i])

    # create constraint for each pair and donor 
    for c in cycles_chains_with_vertex + chains_with_donor:
        if len(c) > 0:
            problem += lpSum(c) <= 1 # add constraint that each vertex can be used in at most 1 cycle or chain

    # get weights of each cycle and chain - will vary by problem type
    cycle_weights = graph.cycle_weights
    chain_weights = graph.chain_weights
    
    # add objective function
    problem += lpDot([cycle_vars[i] for i in range(len(cycle_vars))] + [chain_vars[i] for i in range(len(chain_vars))], cycle_weights + chain_weights)

    # solve for optimal solution
    # print('Solving')
    problem.solve(PULP_CBC_CMD(msg=0))

    # prints out optimal objective value
    # print(f'Objective Value: {value(problem.objective)}')

    # gets (indices of) pairs that have been matched
    matched = [pairs[p] for c in range(len(cycle_vars)) for p in cycles[c].pairs if value(cycle_vars[c]) == 1]
    matched = matched + [pairs[p] for c in range(len(chain_vars)) for p in chains[c].pairs if value(chain_vars[c]) == 1]
    # print(f'Number of Matched Pairs: {len(matched)}')

    # get (indices of) altruistic donors that have been used
    used_altruistic_donors = [altruistic_donors[chains[d].altruistic_donor] for d in range(len(chain_vars)) if value(chain_vars[d]) == 1]
    # print(f'Number of Altruistic Donors Used: {len(used_altruistic_donors)}')

    # check to make sure no pair or donor was used twice
    assert len(matched) == len(set(matched))
    assert len(used_altruistic_donors) == len(set(used_altruistic_donors))

    return matched, used_altruistic_donors

if __name__ == "__main__":
    solve_kidney_matching(all_pairs, [generate_patient_donor_pair().donor for _ in range(5)], ProblemType.FAIRNESS, 5)

### code for greedy approach 
# def greedy_solve_kidney_matching(new_pair, existing_pairs, problem_type):
#     pairs = existing_pairs + [new_pair] # new pair is put at the end of existing pairs
#     new_ind = len(pairs) - 1 # index of new pair is at the end of pairs list

#     edges = Graph.find_edges(pairs) # edges for new graph with all pairs
#     cycles = [] # cycles involving new pair

#     # iterate over all pairs the new pair can donate to
#     for p in edges[new_ind]:
#         # check for 2-cycle with p as other element
#         if new_ind in edges[p]:
#             cycles.append(Cycle([new_ind, p]))
#         # check for 3-cycle with p and v as other elements
#         for v in edges[p]:
#             if new_ind in edges[v]:
#                 cycles.append(Cycle([new_ind, p, v]))

#     # get weights of all possible cycles
#     cycle_weights = Graph.find_cycle_weights(problem_type, pairs, cycles)

#     # return cycle with max weight if possible
#     if len(cycles) > 0:
#         selected_cycle_ind = cycle_weights.index(max(cycle_weights))
#         selected = cycles[selected_cycle_ind].pairs
#     else:
#         selected = []

#     # return the pairs that are matched (if possible)
#     print(selected)
#     return selected

# greedy_solve_kidney_matching(generate_patient_donor_pair(), all_pairs, ProblemType.SIMPLE)