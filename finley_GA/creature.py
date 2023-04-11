# The 'Creature' class is used in 'evolve.py' to run a genetic algorithm.
# Maud Formanek, 2020, Francois Nedelec October 2021
# Copyright Sainsbury Laboratory, Cambridge University, UK

try:
    import sys, math, random, os
    from genetics import LocusSet
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n"%str(e))
    sys.exit()


class Creature(object):
    """
    Class representing an individual for the genetic algorithm
    It is defined by its genome (a list of integers) encoding parameter values
    The composition of this genome is defined by the 'code'
    """
    
    def __init__(self, cod, gen, fit = math.nan):
        """ create Creature """
        if isinstance(gen, str):
            #sys.stderr.write(f'genome : {gen}\n')
            s = gen.split("'")
            gen = []
            for i in s:
                gen.append(int(i, 2))
        if not isinstance(gen, list):
            sys.stderr.write(f'Unexpected genome type : {gen}\n')
        if len(gen) != cod.num_genes():
            sys.stderr.write(f'Error: Genome should be of length {cod.num_genes()}!\n')
            sys.exit()
        # member variables:
        self.dic = {}
        self.code = cod
        self.genome = gen
        self.fitness = fit
        self.target = 0
        self.home = ''
    
    def update_home(self, root, id):
        """ Set home path to root/XXXX where XXXX is `id` on 4 digits"""
        self.home = os.path.join(root, str(id).rjust(4, '0'))
    
    def update_dictionary(self, pam):
        """ Set dictionary according to genome, adding values in `pam` """
        self.dic = self.values()
        for k, v in pam.items():
            self.dic[k] = v
    
    def sequence(self):
        """ Return a compact representation of the genome """
        return self.code.genome_to_hex(self.genome)
    
    def __str__(self):
        """ return string representation of bit genomes"""
        return self.code.genome_to_string(self.genome)
    
    def values(self):
        """ return dictionary of parameter values corresponding to genome"""
        return self.code.genome_to_values(self.genome)

    def mate_uniform(self, other):
        """ 
        Produce offspring by randomly selecting bits from each `self` and `other`
        """
        return self.code.random_select(self.genome, other.genome)
    
    def mate_crossover(self, other):
        """ 
        Combine self and other genome as if they are circular, with 1 crossover
        """
        return self.code.mate_crossover(self.genome, other.genome)
    
    def mutated(self, cnt: float):
        """ Return new sequence with approximately 'cnt' bits flipped """
        return self.code.mutated(self.genome, cnt)
    
    def mutate(self, cnt: float):
        """
            Mutate approximately 'cnt' bits
        """
        self.genome = self.mutated(cnt)


