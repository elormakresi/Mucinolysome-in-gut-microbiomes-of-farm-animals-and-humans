#!/bin/bash

# run checkm2 for MAG quality check
bin/checkm2 predict --threads 30 --input ./nasal_fecal_wisc_farmers --output-directory ./nasal_fecal_wisc_farmers/checkm2_results_2 --database_path ./checkm2/CheckM2_database/uniref100.KO.1.dmnd
