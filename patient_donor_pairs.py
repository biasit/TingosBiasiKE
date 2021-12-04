from enum import Enum
import random

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


class BloodProteins(Enum):
    O = 1
    A = 2
    B = 3

    # We could consider certain racial distributions as well!
    # Current distributions of proteins based on table A.4 from https://web.stanford.edu/~iashlagi/papers/MS-kidney_exchange.pdf
    def generate_general_distribution_protein():
        ran = random.random()

        if ran <= 0.7022:
            return BloodProteins.O
        elif ran <= 0.7022 + 0.1998:
            return BloodProteins.A
        else:
            return BloodProteins.B


# Another generation function for generating PRA based on population distributions
def generate_general_distribution_pra():
        ran = random.random()

        # Distribution
        l1 = 0.6256
        l2 = 0.1648
        l3 = 0.069
        l4 = 0.0506
        l5 = 0.0274
        l6 = 0.0626
        
        if ran <= l1:
            return 0
        elif ran <= l1 + l2:
            return 30
        elif ran <= l1 + l2 + l3:
            return 65
        elif ran <= l1 + l2 + l3 + l4:
            return 87
        elif ran <= l1 + l2 + l3 + l4 + l5:
            return 97
        else:
            return 99.5


class Patient():
    def __init__(self, blood_type, pra):
        self.blood_type = blood_type
        self.pra = pra
    
    def is_compatible_with_donor(self, donor):
        return BloodType.can_donor_donate_to_patient(donor.blood_type, self.blood_type) and (donor.virtual_pra > self.pra)

class Donor():
    def __init__(self, blood_type):
        self.blood_type = blood_type
        self.virtual_pra = random.random()
    
    def is_compatible_with_patient(self, patient):
        return BloodType.can_donor_donate_to_patient(self.blood_type, patient.blood_type) and (self.virtual_pra > patient.pra)


# A function for generating a new patient and donor pair
def generate_patient_donor_pair():
    
