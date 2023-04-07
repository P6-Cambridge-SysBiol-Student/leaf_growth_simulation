# Provides methods for the artificial evolution
# Francois Nedelec January 2022
# Copyright Sainsbury Laboratory, Cambridge University, UK

try:
    import os, sys, math, re, csv
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n" % str(e))
    sys.exit()

# simulation executable must be in current working directory:
simex = os.path.abspath('frap')

# name of file where simulation provides its results, or 'stdout':
simout = 'outputFourierCoeffs.csv'

# template file used to generate simulation configuration files
template = 'config.cym.tpl'

# Parameters of the genetics: parameter names and ranges
genetic_config = 'genetics.config'

#-------------------------------------------------------------------------------


def load_target():
    """return target for optimization, which can be a complex data"""
    filename = 'target.txt'
    return 0


#-------------------------------------------------------------------------------

def calculate_fitness(file_path):
    """
    Calculate fitness expressing the relative contribution of the fourth Fourier coefficient.
    Higher fourth coefficient means higher fitness, while non-4 coefficients negatively affect the score.
    """
    fit = 0  # Initialize the fitness value

    # Open the CSV file in read mode using a context manager
    with open(file_path, 'r') as csvfile:
        # Create a CSV reader object to read the CSV file line by line
        reader = csv.reader(csvfile)

        # Iterate through each row in the CSV file
        for row in reader:
            if row:  # Check if the row is not empty
                try:
                    coeff_index = int(row[0])  # Parse the first column as an integer (coefficient index)
                    coeff_value = float(row[4])  # Parse the fourth column as a float (coefficient magnitude)

                    # Check if the current row corresponds to the fourth Fourier coefficient (index 4)
                    if coeff_index == 4:
                        fit += coeff_value  # Add the coefficient value to the fitness value
                    elif (coeff_index !=4) and (coeff_value!=0):
                        fit -= coeff_value  # Subtract the coefficient value from the fitness value for non-4 coefficients
                except ValueError:
                    pass  # Ignore lines that cannot be parsed (e.g., headers, non-numeric data)
    return fit  # Return the calculated fitness value

