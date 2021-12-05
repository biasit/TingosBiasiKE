from pulp import LpProblem, LpVariable, LpMaximize, value, lpSum, lpDot
from enum import Enum

from patient_donor_pairs import generate_patient_donor_pair

# generate patient-donor pairs
number_of_pairs = 2000
all_pairs = []
for i in range(number_of_pairs):
    all_pairs.append(generate_patient_donor_pair())

class ProblemType(Enum):
    SIMPLE = 1
    POTENTIALS = 2
    FAIRNESS = 3

# cycle class
class Cycle:
    def __init__(self, pairs):
        self.pairs = pairs
        self.size = len(self.pairs)

class Graph:
    def __init__(self, pairs, problem_type):
        self.pairs = pairs
        self.edges = self.find_edges(self.pairs)
        self.cycles = self.find_cycles(self.pairs, self.edges)
        self.problem_type = problem_type
        self.cycle_weights = self.find_weights(self.problem_type, self.pairs, self.cycles)

    # create adjacency list representation of graph of pairs
    def find_edges(self, pairs):
        edges = [set() for _ in range(len(pairs))]

        for i in range(len(pairs)):
            for j in range(i+1, len(pairs)):
                if pairs[i][1].is_compatible_with_patient(pairs[j][0]):
                    edges[i].add(j)
                if pairs[j][1].is_compatible_with_patient(pairs[i][0]):
                    edges[j].add(i)
        
        return edges

    # find all 2 and 3 cycles for graph of pairs
    def find_cycles(self, pairs, edges):
        cycles = []

        # vars for finding 3-cycles
        unexplored = 0
        explored = 1
        done = 2
        status = [unexplored for _ in range(len(pairs))]
        parent = {}

        # dfs to find 3-cycles
        def dfs(u):
            if status[u] == done:
                return

            status[u] = explored
            for v in edges[u]:
                # Back edge
                if status[v] == explored:
                    if parent[parent[u]] == v: 
                        cycles.append(Cycle([u, v, parent[u]]))
                # v unseen
                elif status[v] == unexplored:
                    parent[v] = u
                    dfs(v)      
            status[u] = done

        # find all 3-cycles
        for i in range(len(pairs)):
            dfs(i)

        # find all 2-cycles
        for i in range(len(pairs)):
            for j in range(i+1, len(pairs)):
                if j in edges[i] and i in edges[j]: 
                    c = Cycle([i, j])
                    cycles.append(c)

        return cycles

    # function that establishes the optimization weights for each cycle based on the problem type
    def find_weights(self, problem_type, pairs, cycles):
        # if simple, weights are size of the cycle
        if problem_type == ProblemType.SIMPLE:
            return [c.size for c in cycles]
        # if potentials, weights are size of cycle minus potential of each vertex in cycle
        elif problem_type == ProblemType.POTENTIALS:
            return [c.size - sum([pairs[p][0].potential + pairs[p][1].potential for p in c.pairs]) for c in cycles]
        # if fairness, ... TODO
        elif problem_type == ProblemType.FAIRNESS:
            return [c.size for c in cycles] # change once fairness is established

def solve_kidney_matching(pairs, problem_type):
    # construct graph
    graph = Graph(pairs, problem_type)

    # get cycles
    cycles = graph.cycles

    # initialize problem
    problem = LpProblem('kidney_matching', LpMaximize)

    # create decision variables for each cycle
    cycle_vars = LpVariable.dicts('cycle', range(len(cycles)), cat='Binary')

    # create constraint for each pair
    cycles_with_vertex = [[] for _ in range(len(cycles))] # cycles_with_vertex[i] is list of cycles that contain vertex i
    for i in range(len(cycles)):
        for j in cycles[i].pairs:
            cycles_with_vertex[j].append(cycle_vars[i])
    for c in cycles_with_vertex:
        if len(c) > 0:
            problem += lpSum(c) <= 1 # add constraint that each vertex can be used in at most 1 cycle

    # get weights of each cycle - will vary by problem type
    cycle_weights = graph.cycle_weights
    
    # add objective function
    problem += lpDot([cycle_vars[i] for i in range(len(cycle_vars))], cycle_weights)

    # solve for optimal solution
    problem.solve()

    # prints out number of cycles and optimal objective value
    print(f'Number of Cycles in Graph: {len(cycles)}')
    print(f'Objective Value: {value(problem.objective)}')

    # gets all the (indices of) pairs that have been matched
    matched = [p for c in range(len(cycle_vars)) for p in cycles[c].pairs if value(cycle_vars[c]) == 1]
    # print(matched)
    print(f'Number of Matched Pairs: {len(matched)}')

    return matched

solve_kidney_matching(all_pairs, ProblemType.SIMPLE)