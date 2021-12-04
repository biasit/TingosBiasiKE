# This is an example file for how to use patient-donor pair generation code
from patient_donor_pairs import generate_patient_donor_pair

number_of_pairs = 2000
all_pairs = []

for i in range(number_of_pairs):
    all_pairs.append(generate_patient_donor_pair())

# Our compatibility graph
for patient, donor in all_pairs:
    print(f"Patient {patient.blood_type}, Donor {donor.blood_type}")