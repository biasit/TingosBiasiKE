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
from enum import Enum
import time

from patient_donor_pairs import generate_patient_donor_pair, generate_altruistic_donor, Donor, Pair, BloodType
from solver import solve_kidney_matching

# Will likely want to introduce a seed at some point
class ExponentialDistribution():
    def __init__(self, rate):
        self.rate = rate 
    
    def draw(self):
        return -(math.log(random.random()) / self.rate) # there's a name in Stat 110 for this, but you basically invert the CDF of the exponential distribution

# Altruists vs Pair
class Vertex(Enum):
    Pair = 1
    Altruist = 2

# General functions for calculating statistics
def calculate_average_waiting_time(vertices):       # average waiting time for people that were matched
    if len(vertices) == 0:
        return float('nan')
    
    total_wait_time = 0
    for vertex in vertices:
        total_wait_time += vertex.match_time
    
    return total_wait_time / len(vertices)

def calculate_pra_prop_matched(pairs, pra):
    num_high_pra = 0
    num_high_pra_matched = 0
    for pair in pairs:
        patient = pair.patient
        if patient.pra >= pra:
            num_high_pra += 1

            if pair.was_matched:
                num_high_pra_matched += 1
    if num_high_pra == 0:
        return float('nan')

    return num_high_pra_matched / num_high_pra

def calculate_prop_matched_expiration(pairs, expiration_time):
    num_high_expiration = 0
    num_high_expiration_matched = 0
    for pair in pairs:
        if pair.departure_time - pair.arrival_time <= expiration_time:
            num_high_expiration += 1

            if pair.was_matched:
                num_high_expiration_matched += 1
    
    if num_high_expiration == 0:
        return float('nan')
    
    return num_high_expiration_matched / num_high_expiration

def calculate_prop_matched_blood_type(pairs, blood_type):
    num_blood_type = 0
    num_blood_type_matched = 0
    for pair in pairs:
        patient = pair.patient
        if patient.blood_type == blood_type:
            num_blood_type += 1

            if pair.was_matched:
                num_blood_type_matched += 1

    if num_blood_type == 0:
        return float('nan')

    return num_blood_type_matched / num_blood_type

class DynamicSimulator():
    def __init__(self, pair_arrival_rate, pair_departure_rate, altruist_arrival_rate, altruist_departure_rate, 
                    problem_type, batch_size=1):
        self.pair_arrival_rate = pair_arrival_rate           # Poisson(arrival_rate) number of pairs arriving every time period
        self.pair_departure_rate = pair_departure_rate         # Exp(survival_rate) - lifespan of a pair in the donor pool
        self.altruist_arrival_rate = altruist_arrival_rate
        self.altruist_departure_rate = altruist_departure_rate

        self.problem_type = problem_type                # solver problem type
        self.batch_size = batch_size                    # If Batch frequency, use batch_size


        self.pair_arrival_generator = ExponentialDistribution(self.pair_arrival_rate)
        self.pair_survival_generator = ExponentialDistribution(self.pair_departure_rate)
        self.altruist_arrival_generator = ExponentialDistribution(self.altruist_arrival_rate)
        self.altruist_departure_generator = ExponentialDistribution(self.altruist_departure_rate)

    def run(self, time_limit):
        """
            Run the simulation.

            time_limit: how long to run the simulation for (this many time periods)
        """
        start_time = time.time()

        entry_count = 0               # heapq breaks without the entry_count for ties
        print()
        print()
        print("Simulator Starting")

        # Preload the pair arrival times, departure times
        pair_arrival_times = deque()
        pair_departure_times = deque()

        curr_time = 0.0
        while True:
            # Get the arrival and departure time from exponential distributions
            entry_time = curr_time + self.pair_arrival_generator.draw()
            if self.pair_departure_rate > 0:
                exit_time = entry_time + self.pair_survival_generator.draw()
            else:
                exit_time = float('inf')

            curr_time = entry_time
            
            if curr_time > time_limit:  # we have reached the limit of our time
                break 

            pair_arrival_times.append(entry_time)
            pair_departure_times.append(exit_time)

        # Preload the altruistic donor arrival times, departure times (if necessary)
        altruist_arrival_times = deque()
        altruist_departure_times = deque()

        if self.altruist_arrival_rate > 0:
            curr_time = 0.0
            while True:
                # Get the arrival and departure time of altruist from exponential distribution
                entry_time = curr_time + self.altruist_arrival_generator.draw()
                if self.altruist_departure_rate > 0:
                    exit_time = entry_time + self.altruist_departure_generator.draw()
                else:
                    exit_time = float('inf')
                curr_time = entry_time
                
                if curr_time > time_limit:  # reach limit of time
                    break

                altruist_arrival_times.append(entry_time)
                altruist_departure_times.append(exit_time)

        print(f"This simulator will involve {len(pair_arrival_times)} pairs")
        print(f"This simulator will involve {len(altruist_departure_times)} altruistic donors")
        
        # Combine the arrival times together
        arrival_times = deque()
        departure_times = deque()
        while len(pair_arrival_times) > 0 and len(altruist_arrival_times) > 0:
            if pair_arrival_times[0] <= altruist_arrival_times[0]:
                arrival_times.append((Vertex.Pair, pair_arrival_times[0]))
                departure_times.append((Vertex.Pair, pair_departure_times[0]))      # not that in case tie with altruist, departure times of pairs come first
                pair_arrival_times.popleft()
                pair_departure_times.popleft()
            else:
                arrival_times.append((Vertex.Altruist, altruist_arrival_times[0]))
                departure_times.append((Vertex.Altruist, altruist_departure_times[0]))
                altruist_arrival_times.popleft()
                altruist_departure_times.popleft()
        while len(pair_arrival_times) > 0:
            arrival_times.append((Vertex.Pair, pair_arrival_times[0]))
            departure_times.append((Vertex.Pair, pair_departure_times[0]))      # not that in case tie with altruist, departure times of pairs come first
            pair_arrival_times.popleft()
            pair_departure_times.popleft()
        while len(altruist_arrival_times) > 0:
            arrival_times.append((Vertex.Altruist, altruist_arrival_times[0]))
            departure_times.append((Vertex.Altruist, altruist_departure_times[0]))
            altruist_arrival_times.popleft()
            altruist_departure_times.popleft()

        # Track the current state of the pool
        pair_pool = set()
        altruist_pool = set()

        vertices_by_exit_time = []  # this will be a priority queue of tuples (exit_time, entry_count, (Patient, Donor) pair or altruistic donor)
        all_matched_pairs = set()
        all_matched_altruists = set()
        
        # For stats purposes
        all_expired_pairs = set()
        all_expired_altruists = set()

        # General statistics about the process
        total_pairs_matched = 0
        total_altruists_matched = 0
        total_pairs_expired = 0
        total_altruists_expired = 0
        total_pairs_seen = 0                          # total number of pairs that arrive to the exchange
        total_altruists_seen = 0

        # Simulate everything!
        curr_time = 0.0 
        curr_batch = 0              # if matching with batches, matches whenever curr_batch >= self.batch_size
        while True:
            # Remove vertices between the current time and the next entry time - these go unmatched for now, we can change this...
            if len(arrival_times) > 0 :
                next_entry_time = arrival_times[0][1]
                while len(vertices_by_exit_time) != 0 and min(vertices_by_exit_time)[0] <= next_entry_time:  # a vertex expires before next vertex arrives
                    # Get the vertex that is leaving and make sure it hasn't been matched already
                    critical_vertex = heapq.heappop(vertices_by_exit_time)[2]
                    
                    if type(critical_vertex) == Pair:
                        if critical_vertex not in pair_pool:
                            continue
                        else:
                            pair_pool.remove(critical_vertex)   # for now, just remove from pool, we will want to probably match these though (can discuss this)
                            all_expired_pairs.add(critical_vertex)
                            total_pairs_expired += 1
                    elif type(critical_vertex) == Donor:
                        if critical_vertex not in altruist_pool:
                            continue
                        else:
                            altruist_pool.remove(critical_vertex)
                            all_expired_altruists.add(critical_vertex)
                            total_altruists_expired += 1

            
            # If no new vertices to enter, we are finished!
            if len(arrival_times) == 0:
                break

            # Otherwise, simulate arrivals
            curr_time = next_entry_time

            # Determine number of arrivals at this time (may be more than 1)
            new_pair_arrivals = 0
            new_altruist_arrivals = 0
            while len(arrival_times) != 0 and arrival_times[0][1] == curr_time:
                curr_arrival_type = arrival_times.popleft()[0]
                if curr_arrival_type == Vertex.Pair:
                    new_pair_arrivals += 1
                elif curr_arrival_type == Vertex.Altruist:
                    new_altruist_arrivals += 1

            # Update the number of pairs that have arrived
            total_pairs_seen += new_pair_arrivals
            total_altruists_seen += new_altruist_arrivals
            curr_batch += (new_pair_arrivals + new_altruist_arrivals)

            # Generate the new pairs
            new_pairs = set()
            for _ in range(new_pair_arrivals):
                curr_pair = generate_patient_donor_pair()
                curr_pair.arrival_time = curr_time 

                # Obtain departure time for pair
                departure_entry = departure_times.popleft()
                assert(departure_entry[0] == Vertex.Pair)       # quick sanity check
                departure_time = departure_entry[1]
                curr_pair.departure_time = departure_time

                # track when it will be leaving the simulation
                heapq.heappush(vertices_by_exit_time, (departure_time, entry_count, curr_pair))
                entry_count += 1

                new_pairs.add(curr_pair)

            # Generate the new altruistic donors
            new_altruists = set()
            for _ in range(new_altruist_arrivals):
                curr_donor = generate_altruistic_donor()
                curr_donor.arrival_time = curr_time

                # Obtain departure time
                departure_entry = departure_times.popleft()  # in case of tie between pair and donor, departure_time of donor is after pair...
                assert(departure_entry[0] == Vertex.Altruist)       # quick sanity check
                departure_time = departure_entry[1]
                curr_donor.departure_time = departure_time

                # track when it will be leaving the simulation
                heapq.heappush(vertices_by_exit_time, (departure_time, entry_count, curr_donor))
                entry_count += 1

                new_altruists.add(curr_donor)

            # Add the new vertices to the pools
            pair_pool |= new_pairs
            altruist_pool |= new_altruists

            # Undergo matching algorithm if necessary
            matched_pairs = None
            matched_donors = None
            if curr_batch >= self.batch_size:
                matched_pairs, matched_donors = solve_kidney_matching(list(pair_pool), list(altruist_pool), self.problem_type, curr_time)
                curr_batch = 0


            # Remove any matched pairs
            if matched_pairs is not None:
                for pair in matched_pairs:
                    pair.was_matched = True
                    pair.match_time = curr_time
                    pair_pool.remove(pair)
                    all_matched_pairs.add(pair)
                total_pairs_matched += len(matched_pairs)
                
            
            # Remove any matched donors
            if matched_donors is not None:
                for donor in matched_donors:
                    donor.was_matched = True
                    donor.match_time = curr_time
                    altruist_pool.remove(donor)
                    all_matched_altruists.add(donor)
                total_altruists_matched += len(matched_donors)

        end_time = time.time()
        print()
        print(f"Total time of simulation: {round((end_time - start_time) / 60, 3)} minutes")
        print()

        # Collect helpful statistics in dictionary
        original_pair_pool = pair_pool | all_expired_pairs | all_matched_pairs
        original_altruist_pool = altruist_pool | all_expired_altruists | all_matched_altruists

        statistics = {}
        statistics["Number of Pairs Matched"] = total_pairs_matched
        statistics["Number of Pairs Seen"] = total_pairs_seen
        statistics["Number of Pairs Expired"] = total_pairs_expired
        statistics["Proportion of Pairs Matched"] = total_pairs_matched / len(original_pair_pool)
        statistics["Proportion of Pairs Expired"] = total_pairs_expired / len(original_pair_pool)
        statistics["Proportion of Pairs Left At End"] = len(pair_pool) / len(original_pair_pool)
        statistics["Pair Average Wait Time"] = calculate_average_waiting_time(all_matched_pairs)

        statistics["Number of Altruists Matched"] = total_altruists_matched
        statistics["Number of Altruists Seen"] = total_altruists_seen
        statistics["Number of Altruists Expired"] = total_altruists_expired

        if len(original_altruist_pool) > 0:
            statistics["Proportion of Altruists Matched"] = total_altruists_matched / len(original_altruist_pool)
            statistics["Proportion of Altruists Left At End"] = len(altruist_pool) / len(original_altruist_pool)
            statistics["Proportion of Altruists Expired"] = total_altruists_expired / len(original_altruist_pool)
        statistics["Altruist Average Wait Time"] = calculate_average_waiting_time(all_matched_altruists)

        # Calculate fairness statistics
        statistics["Proportion with PRA > 0.05 Matched"] = calculate_pra_prop_matched(original_pair_pool, 0.05)
        statistics["Proportion with PRA > 0.2 Matched"] = calculate_pra_prop_matched(original_pair_pool, 0.2)
        statistics["Proportion with PRA > 0.4 Matched"] = calculate_pra_prop_matched(original_pair_pool, 0.4)
        statistics["Proportion with PRA > 0.6 Matched"] = calculate_pra_prop_matched(original_pair_pool, 0.6)
        statistics["Proportion with PRA > 0.8 Matched"] = calculate_pra_prop_matched(original_pair_pool, 0.8)
        statistics["Proportion with PRA > 0.9 Matched"] = calculate_pra_prop_matched(original_pair_pool, 0.9)

        statistics["Proportion of Type O Matched"] = calculate_prop_matched_blood_type(original_pair_pool, BloodType.O)
        statistics["Proportion of Type A Matched"] = calculate_prop_matched_blood_type(original_pair_pool, BloodType.A)
        statistics["Proportion of Type B Matched"] = calculate_prop_matched_blood_type(original_pair_pool, BloodType.B)
        statistics["Proportion of Type AB Matched"] = calculate_prop_matched_blood_type(original_pair_pool, BloodType.AB)

        statistics["Proportion with expiration 0.01 Matched"] = calculate_prop_matched_expiration(original_pair_pool, 0.01)
        statistics["Proportion with expiration 0.05 Matched"] = calculate_prop_matched_expiration(original_pair_pool, 0.05)
        statistics["Proportion with expiration 0.1 Matched"] = calculate_prop_matched_expiration(original_pair_pool, 0.1)
        statistics["Proportion with expiration 0.25 Matched"] = calculate_prop_matched_expiration(original_pair_pool, 0.2)
        statistics["Proportion with expiration 0.5 Matched"] = calculate_prop_matched_expiration(original_pair_pool, 0.5)
        
        print()
        print("RESULTS")
        print()
        for key in statistics:
            print(f"{key}: {round(statistics[key], 3)}")
        
        return all_matched_pairs, all_expired_pairs, all_matched_altruists, all_expired_altruists, statistics

