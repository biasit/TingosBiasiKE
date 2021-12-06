# Here we should run all of the experiments that we need to run

from simulator import DynamicSimulator
from solver import ProblemType


# Example usage of simulator (no altruists)
print()
print("Example Simulator with No Altruists, batch size of 10 arrivals (pair or altruist)")
simulator = DynamicSimulator(pair_arrival_rate=100, pair_departure_rate=2, 
                            altruist_arrival_rate=0, altruist_departure_rate=0, 
                            problem_type=ProblemType.SIMPLE, batch_size=10)
_, _, _, _ = simulator.run(10)

# Example usage of simulator (altruists)
print()
print("Example Simulator with Altruists, batch size of 1 (greedy)")
simulator = DynamicSimulator(pair_arrival_rate=100, pair_departure_rate=2, 
                            altruist_arrival_rate=10, altruist_departure_rate=0, 
                            problem_type=ProblemType.SIMPLE, batch_size=1)
_, _, _, _ = simulator.run(10)