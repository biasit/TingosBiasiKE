from pulp import LpProblem, LpVariable, LpMaximize, value, lpSum, lpDot

from patient_donor_pairs import generate_patient_donor_pair

# generate patient-donor pairs
number_of_pairs = 500
all_pairs = []
for i in range(number_of_pairs):
    all_pairs.append(generate_patient_donor_pair())

# cycle class
class Cycle:
    def __init__(self, pairs):
        self.pairs = pairs
        self.size = len(self.pairs)

# create adjacency list representation of graph of pairs
def find_edges(pairs):
    edges = [set() for _ in range(len(pairs))]

    for i in range(len(pairs)):
        for j in range(i+1, len(pairs)):
            if pairs[i][1].is_compatible_with_patient(pairs[j][0]):
                edges[i].add(j)
            if pairs[j][1].is_compatible_with_patient(pairs[i][0]):
                edges[j].add(i)
    
    return edges

# find all 2 and 3 cycles for graph of pairs
def find_cycles(pairs):
    edges = find_edges(pairs)
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

def solve_kidney_matching(pairs):
    # get cycles
    cycles = find_cycles(pairs)

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
        problem += lpSum(c) <= 1 # add constraint that each vertex can be used in at most 1 cycle

    # initialize weights of each cycle
    # adjust this for when we optimize for different things
    cycle_weights = [c.size for c in cycles]
    
    # add objective function
    problem += lpDot([cycle_vars[i] for i in range(len(cycle_vars))], cycle_weights)

    # solve for optimal solution
    problem.solve()

    # prints out number of cycles and optimal objective value
    print(len(cycles))
    print(value(problem.objective))

solve_kidney_matching(all_pairs)