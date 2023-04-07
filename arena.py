# Provides methods for the artificial evolution
# Francois Nedelec January 2022
# Copyright Sainsbury Laboratory, Cambridge University, UK

try:
    import os, sys, math, re
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n"%str(e))
    sys.exit()

# simulation executable must be in current working directory:
simex = os.path.abspath('sim3')

# name of file where simulation provides its results, or 'stdout':
simout = 'result.txt'

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

def calculate_fitness(data, target):
    """
    Calculate fitness expressing how close `data` is to `target`
    Higher fitness is better.
    This implements fitness for the plasmid partitionning problem.
    """
    fit = 0
    cnt = 0
    if isinstance(data, str):
        raise  Exception("Cannot calculate fitness from string")
    #print("calculate_fitness from %i data groups" % len(data))
    for out in data:
        for line in out.split('\n'):
            if len(line) > 0 and line[0] != '%':
                fit += float(line)
                cnt += 1
    if fit:
        fit = cnt / abs(fit)
    else:
        fit = 1000000
    #print("fitness = %.4f from %i datalines" % (fit, x[0]))
    return fit

