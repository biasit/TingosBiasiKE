# Here we should run all of the experiments that we need to run

from simulator import DynamicSimulator
from solver import ProblemType
import random
import time

def run_experiment(number_of_repetitions,
                  pair_arrival_rate, pair_departure_rate,
                  altruist_arrival_rate=0, altruist_departure_rate=0,
                  problem_type=ProblemType.SIMPLE, batch_size=10,
                  time_limit=10):
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
    print("EXPERIMENT RESULTS")
    print(f"{number_of_repetitions} Repetitions")
    print()
    print("Averaged Statistics:")
    for key in averaged_statistics:
        print(f"Average {key}: {round(averaged_statistics[key], 4)}")
    
    end_time = time.time()
    print()
    print(f"Total Time for Experiment: {round((end_time - start_time) / 60, 3)} minutes")


# Example experiment
run_experiment(5, 100, 1, time_limit=20)


