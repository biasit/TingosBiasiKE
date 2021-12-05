"""
This code is heavily based on the dynamic simulator from Dickinson provided here:
https://github.com/JohnDickerson/KidneyExchange
It was very helpful :)
But it was written in Java, and it didn't account for altruistic donors.

"""

from collections import deque # deques are always nice, hopefully will provide a slight speed up
import heapq                  # just priority queue functionality, will allow us to order patients
import random
import math

from patient_donor_pairs import generate_patient_donor_pair

# Will likely want to introduce a seed at some point
class ExponentialDistribution():
    def __init__(self, rate):
        self.rate = rate 
    
    def draw(self):
        return -(math.log(random.random()) / self.rate) # there's a name in Stat 110 for this, but you basically invert the CDF of the exponential distribution



class DynamicSimulator():
    def __init__(self, arrival_rate, survival_rate, matching_algorithm):
        self.arrival_rate = arrival_rate           # Poisson(arrival_rate) number of pairs arriving every time period
        self.survival_rate = survival_rate         # Exp(survival_rate) - lifespan of a pair in the donor pool
        self.matching_algorithm = matching_algorithm

        self.arrival_generator = ExponentialDistribution(self.arrival_rate)
        self.survival_generator = ExponentialDistribution(self.survival_rate)

    def run(self, time_limit):
        """
            Run the simulation.

            time_limit: how long to run the simulation for (this many time periods)
        """

        print("Simulator Starting")

        # Preload the arrival times, departure times
        arrival_times = deque()
        departure_times = deque()

        curr_time = 0.0
        while True:
            # Get the arrival and departure time from exponential distributions
            entry_time = curr_time + self.arrival_generator.draw()
            exit_time = entry_time + self.survival_generator.draw()
            curr_time = entry_time
            
            if curr_time > time_limit:  # we have reached the limit of our time
                break 

            arrival_times.append(entry_time)
            departure_times.append(exit_time)

        print(f"This simulator will involve {len(arrival_times)} pairs")

        # Track the current state of the pool
        pool = set()
        vertices_by_exit_time = []  # this will be a priority queue of tuples (exit_time, (Patient, Donor) pair)
        all_matched_pairs = set()
        all_expired_pairs = set()

        # General statistics about the process
        total_matched = 0
        total_expired = 0
        total_seen = 0                          # total number of pairs that arrive to the exchange

        # Simulate everything!
        curr_time = 0.0 
        while True:
            # Remove vertices between the current time and the next entry time - these go unmatched for now, we can change this...
            if len(arrival_times) == 0:
                next_entry_time = float('inf')
            else:
                next_entry_time = arrival_times.popleft()

            while len(vertices_by_exit_time) != 0 and min(vertices_by_exit_time)[0] <= next_entry_time:  # a vertex expires before next vertex arrives
                # Get the vertex that is leaving and make sure it hasn't been matched already
                critical_vertex = heapq.heappop(vertices_by_exit_time)[1]
                if critical_vertex not in pool:
                    continue
                else:
                    pool.remove(critical_vertex)   # for now, just remove from pool, we will want to probably match these though (can discuss this)
                    all_expired_pairs.add(critical_vertex)
                    total_expired += 1
            
            # If no new vertices to enter, we are finished!
            if len(arrival_times) == 0:
                break

            # Otherwise, simulate arrivals
            curr_time = arrival_times.popleft()     # move current time to the next arrival

            # Determine number of arrivals at this time (may be more than 1)
            new_arrivals = 1
            while len(arrival_times) != 0 and arrival_times[0] == curr_time:
                arrival_times.popleft()
                new_arrivals += 1

            # Update the number of pairs that have arrived
            total_seen += new_arrivals

            # Generate the new vertices
            new_pairs = set()
            for _ in range(new_arrivals):
                curr_pair = generate_patient_donor_pair()
                curr_pair.arrival_time = curr_time 

                # Obtain departure time for pair
                departure_time = departure_times.popleft()
                curr_pair.departure_time = departure_time

                # track when it will be leaving the simulation
                heapq.heappush(vertices_by_exit_time, (departure_time, curr_pair))

                new_pairs.add(curr_pair)

            # Match the new vertices
            matched_pairs = self.matching_algorithm(new_pairs, pool)

            # Add the new vertices to the pool
            pool |= new_pairs

            # Remove any matched pairs
            for pair in matched_pairs:
                pool.remove(pair)
                all_matched_pairs.add(pair)

            # Update stats
            total_matched += len(matched_pairs)
    
        print()
        print("RESULTS")
        print("Total matched:", total_matched)
        print("Total seen:", total_seen)
        print("Total expired:", total_expired)

        return all_matched_pairs, all_expired_pairs



# Example usage of simulator
def matching_algorithm(new_pairs, pool):
    new_pool = pool | new_pairs

    return [random.choice(list(new_pool))]

simulator = DynamicSimulator(100, 1, matching_algorithm)
m_pairs, e_pairs = simulator.run(100)