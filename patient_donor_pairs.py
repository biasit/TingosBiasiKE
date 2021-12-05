from enum import Enum
import random

# Global variable for storing the composition of the NKR pool
pra_intervals = [0, 0.05, 0.3, 0.65, 0.875, 0.97, 0.995]
nkr_pool_composition = {}

class BloodType(Enum):
    O = 1
    A = 2
    B = 3
    AB = 4

    # I should probably check this lol
    def can_donor_donate_to_patient(donor_type, patient_type):
        if patient_type == BloodType.O:
            return donor_type == BloodType.O
        elif patient_type == BloodType.A:
            return donor_type == BloodType.A or donor_type == BloodType.O
        elif patient_type == BloodType.B:
            return donor_type == BloodType.B or donor_type == BloodType.O
        else:
            return True # AB is universal receiver I believe

    # Get blood type from string class method
    def get_blood_type_from_string(type_str):
        if type_str == "O":
            return BloodType.O
        elif type_str == "A":
            return BloodType.A
        elif type_str == "B":
            return BloodType.B
        else:
            return BloodType.AB


# We will probably want to add further details such as the longetivity of the patient or even the patient's race
class Patient():
    def __init__(self, blood_type, pra):
        self.blood_type = blood_type
        self.pra = pra
    
    # Donor's virtual pra must be higher than the patient's PRA to be compatible - reflects tissue-type compatibility
    def is_compatible_with_donor(self, donor):
        return BloodType.can_donor_donate_to_patient(donor.blood_type, self.blood_type) and (donor.virtual_pra > self.pra)

class Donor():
    def __init__(self, blood_type, virtual_pra):
        self.blood_type = blood_type
        self.virtual_pra = virtual_pra
    
    # Once again, Donor's virtual pra must be higher than the patient's PRA to be compatible
    def is_compatible_with_patient(self, patient):
        return BloodType.can_donor_donate_to_patient(self.blood_type, patient.blood_type) and (self.virtual_pra > patient.pra)


# A function for generating a new patient and donor pair
def generate_patient_donor_pair():
    # Continue until a valid pair is generated
    while True:
        # Generate ABO of patient donor pair based on NKR Pool Composition
        abo_pair_ran = random.random()
        abo_pair = None
        running_sum = 0.
        for pair in nkr_pool_composition.keys():
            running_sum += nkr_pool_composition[pair]["probability_of_pair"]
            if abo_pair_ran <= running_sum:
                abo_pair = pair
                break

        if abo_pair is None:
            raise Exception(f"Invalid ABO found, running sum is {running_sum}")
        
        # Generate PRA of Patient
        pra_ran = random.random()
        pra = None
        running_sum = 0.
        for i_c, cond_prob in enumerate(nkr_pool_composition[abo_pair]["cond_probability_pra"]):
            running_sum += cond_prob 
            if pra_ran <= running_sum:
                pra = pra_intervals[i_c]
                break
        if not pra:
            pra = pra_intervals[-1]
        
        # Generate virtual pra of donor
        virtual_pra = random.random() # uniform between 0 and 1
        
        # If virtual_pra is not less than pra, then we need to redraw
        if not (virtual_pra < pra) and BloodType.can_donor_donate_to_patient(abo_pair[1], abo_pair[0]):
            continue

        # Return the patient and donor pair
        patient = Patient(abo_pair[0], pra)
        donor = Donor(abo_pair[1], virtual_pra)
        break



    return patient, donor


# Read in the distribution of pairs from the NKR Pool Composition (2010-2014)'
print("Generating NKR Pool Composition")
with open("distributions.txt", "r") as f:
    lines = f.readlines()
    for line in lines:
        line = line.split()

        patient_type, donor_type = line[0].split('-')[0], line[0].split('-')[1]
        patient_type = BloodType.get_blood_type_from_string(patient_type)
        donor_type = BloodType.get_blood_type_from_string(donor_type)

        # Get data for current pair
        curr_pair_data = {}
        curr_pair_data["probability_of_pair"] = float(line[1]) / 100
        curr_pair_data["cond_probability_pra"] = [float(i) / 100 for i in line[2:]]

        # Add data for current pair to the table
        nkr_pool_composition[(patient_type, donor_type)] = curr_pair_data





# Legacy code - I realized I couldn't figure out Baye's rule for Model 2 - big sad
   # # Gives you the blood type given the proteins someone has
    # def determine_blood_type_from_proteins(proteins):
    #     protein1, protein2 = proteins

    #     if (protein1, protein2) in [(BloodProteins.O, BloodProteins.O)]:
    #         return BloodType.O
    #     elif (protein1, protein2) in [(BloodProteins.O, BloodProteins.A), (BloodProteins.A, BloodProteins.O), (BloodProteins.A, BloodProteins.A)]:
    #         return BloodType.A
    #     elif (protein1, protein2) in [(BloodProteins.O, BloodProteins.B), (BloodProteins.B, BloodProteins.O), (BloodProteins.B, BloodProteins.B)]:
    #         return BloodType.B
    #     else:
    #         return BloodType.AB


# # everyone should have two proteins essentially (O implies absence of A or B - this really is just the genes - O is recessive)
# class BloodProteins(Enum):
#     O = 1
#     A = 2
#     B = 3

#     # We could consider certain racial distributions as well!
#     # Current distributions of proteins based on table A.4 from https://web.stanford.edu/~iashlagi/papers/MS-kidney_exchange.pdf
#     def generate_general_distribution_protein():
#         ran = random.random()

#         if ran <= 0.7022:
#             return BloodProteins.O
#         elif ran <= 0.7022 + 0.1998:
#             return BloodProteins.A
#         else:
#             return BloodProteins.B

# # Another generation function for generating PRA based on population distributions
# def generate_general_distribution_pra():
#         ran = random.random()

#         # Distribution
#         l1 = 0.6256
#         l2 = 0.1648
#         l3 = 0.069
#         l4 = 0.0506
#         l5 = 0.0274
#         l6 = 0.0626
        
#         if ran <= l1:
#             return 0
#         elif ran <= l1 + l2:
#             return 30
#         elif ran <= l1 + l2 + l3:
#             return 65
#         elif ran <= l1 + l2 + l3 + l4:
#             return 87
#         elif ran <= l1 + l2 + l3 + l4 + l5:
#             return 97
#         else:
#             return 99.5