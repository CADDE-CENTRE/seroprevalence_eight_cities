# seroprevalence_eight_cities
Severe Acute Respiratory Infection (SARI) data from SIVEP-Gripe (files INFLUD20-14-02-2022.csv and INFLUD21-14-02-2022.csv) can be downloaded at https://opendatasus.saude.gov.br/dataset/srag-2021-e-2022 and https://opendatasus.saude.gov.br/dataset/srag-2020.

Data dictionary for Bloodbank.rds:
Each line of this file represents a donation. The columns are:
- donor_id: unique id for each donor
- sample_id: unique id for each sample
- blood_center: Blood center where sample was collected. The blood centers are: HEMOAM (Manaus), FPS (São Paulo), HEMOCE (Fortaleza), HEMOMINAS (Belo Horizonte), HEMOPE (Recife), HEMEPAR (Curitiba), HEMOBA (Salvador).
- sex: Sex of the donor
- idade: Age of the donor (-1 if unknown)
- donation_date: Date of sample collection
- month: Number of months passed since January 1st 2020. For instance, samples collected in March 2020 and March 2021 have respectively month = 3 and month = 15.
- result: Signal-to-cutoff obtained with the anti-N assay. A test is considered positive if result > threshold. In the paper, we use threshold = 0.49.
