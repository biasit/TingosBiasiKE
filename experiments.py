# Here we should run all of the experiments that we need to run

from simulator import DynamicSimulator
from solver import ProblemType
import random
import time

def run_experiment(number_of_repetitions,
                  time_limit,
                  pair_arrival_rate, pair_departure_rate=0,
                  altruist_arrival_rate=0, altruist_departure_rate=0,
                  problem_type=ProblemType.SIMPLE, batch_size=10):
    start_time = time.time()

    averaged_statistics = {}

    random.seed(0)
    for _ in range(number_of_repetitions):
        simulator = DynamicSimulator(pair_arrival_rate=pair_arrival_rate, pair_departure_rate=pair_departure_rate, 
                                    altruist_arrival_rate=altruist_arrival_rate, altruist_departure_rate=altruist_departure_rate, 
                                    problem_type=problem_type, batch_size=batch_size)
        _, _, _, _, statistics = simulator.run(time_limit)
    
        if len(averaged_statistics) == 0:
            averaged_statistics = statistics
        else:
            for key in statistics:
                averaged_statistics[key] += statistics[key]
    
    for key in averaged_statistics:
        averaged_statistics[key] /= number_of_repetitions
    
    print()
    print()
    print("EXPERIMENT HYPERPARAMETERS")
    print("Number of reps:", number_of_repetitions)
    print("Time limit:", time_limit)
    print("Pair arrival rate:", pair_arrival_rate)
    print("Pair departure rate:", pair_departure_rate)
    print("Altruist arrival rate:", altruist_arrival_rate)
    print("Altruist departure rate:", altruist_departure_rate)
    print("Problem Type:", problem_type)
    print("Batch size:", batch_size)

    print()
    print("EXPERIMENT RESULTS")
    print(f"{number_of_repetitions} Repetitions")
    print()
    print("Averaged Statistics:")
    for key in averaged_statistics:
        print(f"Average {key}: {round(averaged_statistics[key], 4)}")
    
    end_time = time.time()
    print()
    print(f"Total Time for Experiment: {round((end_time - start_time) / 60, 3)} minutes")


# Baseline hyperparameters
number_of_repetitions = 5
time_limit = 20                 # run no longer than 10 units of time
pair_arrival_rate = 100         # means that 100 pairs are going to be arriving every one time period
base_pair_departure_rate = 0.4      # means that each pairs is expected to last 2.25 units of time
base_altruist_arrival_rate = 1.0     # means 1 altruist in expectation will show up each unit of time
base_altruist_departure_rate = 0.4   # means altruist last 2.25 units of time in expectation


# Batch size experiments (no altruists)
batch_sizes = [10, 20, 30, 50, 100, 1]
for batch_size in batch_sizes:
    run_experiment(number_of_repetitions, time_limit, pair_arrival_rate, pair_departure_rate=base_pair_departure_rate, batch_size=batch_size)


# Pair departure rates experiments (no altruists)
departure_rates = [0.2, 0.4, 0.6, 0.8]
batch_sizes = [10, 30, 50, 1]
for batch_size in batch_sizes:
    for departure_rate in departure_rates:
        run_experiment(number_of_repetitions, time_limit, pair_arrival_rate, pair_departure_rate=departure_rate, batch_size=batch_size)


# Impact of altruists (departure rate selected)
altruist_arrival_rates = [0.5, 1.0]
batch_sizes = [10, 30, 50, 1]
for batch_size in batch_sizes:
    for altruist_arrival_rate in altruist_arrival_rates:
        run_experiment(number_of_repetitions, time_limit, pair_arrival_rate, pair_departure_rate=base_pair_departure_rate, altruist_arrival_rate=altruist_arrival_rate, altruist_departure_rate=base_altruist_departure_rate, batch_size=batch_size)


# Potential and Fairness Weighted
for solver_type in [ProblemType.POTENTIALS, ProblemType.FAIRNESS]:
    for batch_size in [1, 10, 30]:
        run_experiment(number_of_repetitions, time_limit, pair_arrival_rate, pair_departure_rate=base_pair_departure_rate, 
                        altruist_arrival_rate=base_altruist_arrival_rate, altruist_departure_rate=base_altruist_departure_rate, batch_size=batch_size, problem_type=solver_type)









