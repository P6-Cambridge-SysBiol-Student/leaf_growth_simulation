# Provides methods for the artificial evolution
# Francois Nedelec January 2022
# Copyright Sainsbury Laboratory, Cambridge University, UK

try:
    import os, sys, math, re, csv
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n" % str(e))
    sys.exit()

# simulation executable must be in current working directory:
simex = os.path.abspath('/home/finley/SLCU/frapc_with_triangles/GA/frap')

# name of file where simulation provides its results, or 'stdout':
simout = 'outputFourierCoeffs.csv'

# template file used to generate simulation configuration files
template = 'config.cym.tpl'

# Parameters of the genetics: parameter names and ranges
genetic_config = 'genetics.config'

#-------------------------------------------------------------------------------


def load_target():
    # maximum value that the fitness calculation can output is 1 (when the fourth coefficient is the only one)
    return 1


#-------------------------------------------------------------------------------

def calculate_fitness(file_path, target):
    """
    Calculate fitness expressing the relative contribution of the fourth Fourier coefficient.
    Higher fourth coefficient means higher fitness.
    """
    sum_coeff_values = 0  # Initialize the sum of all coefficient magnitudes
    fourth_coeff_value = 0  # Initialize the fourth coefficient value

    # Open the CSV file in read mode
    with open(file_path, 'r') as csvfile:
        # Create a CSV reader object to read the CSV file line by line
        reader = csv.reader(csvfile)

        # Iterate through each row in the CSV file
        for row in reader:
            if row:  # Check if the row is not empty
                try:
                    coeff_index = int(row[0])  # Parse the first column as an integer (coefficient index)
                    coeff_value = float(row[4])  # Parse the fourth column as a float (coefficient magnitude)
                    if coeff_index != 0:
                        sum_coeff_values += coeff_value  # Add the coefficient value to the sum of coeff magnitudes
                    # Check if the current row corresponds to the fourth Fourier coefficient (index 4)
                    if coeff_index == 4:
                        fourth_coeff_value = coeff_value  # Set the fourth coefficient value

                except ValueError:
                    pass  # Ignore lines that cannot be parsed (e.g., headers, non-numeric data)

    # Calculate the fitness value by dividing the fourth coefficient value by the sum of all coefficient magnitudes
    if sum_coeff_values > 0:
        fit = fourth_coeff_value / sum_coeff_values
    else:
        fit = 0
    return fit


