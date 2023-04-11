# `Locus` and `LocusSet` are used by 'evolve.py' to run a genetic algorithm.
# Maud Formanek, 2020, Francois Nedelec October 2021 -- Jan 2022, 09.2022
# Copyright Sainsbury Laboratory, Cambridge University, UK

try:
    import sys, os, random
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n"%str(e))
    sys.exit()

err = sys.stderr

#-------------------------------------------------------------------------------

class Locus(object):
    """
        The Locus defines a parameter subject to optimization.
        It includes the parameter name, its range (lowerbound, upperbound)
        and the number of bits used to encode the value
    """
    def __init__(self, name, low, up, nb):
        self.name = name
        self.nbits = nb
        self.lowerbound = low
        self.upperbound = up
    
    def value_to_bits(self, val): # not sure this works or is needed
        """
            return sequence of bits corresponding to given value
        """
        S = float(max(1, 2**self.nbits-1)) / ( self.upperbound - self.lowerbound )
        v = ( val - self.lowerbound ) * S
        return round(v)
    
    def value_from_bits(self, val):
        """
            return value corresponding to given sequence of bits:
            lowerbound + val * ( upperbound - lowerbound ) / (2**nbits-1)
        """
        if self.nbits > 0:
            S = ( self.upperbound - self.lowerbound ) / float(max(1, 2**self.nbits-1))
            v = self.lowerbound + val * S
        else:
            v = self.lowerbound
        i = int(v)
        if i == v:
            return i
        return v

#-------------------------------------------------------------------------------

class LocusSet(object):
    """ A list of Locuses is used to interpret the complete genomic sequence """
    def __init__(self):
        self.Locuses = []
    
    def add(self, name, lower_bound, upper_bound, bitlength):
        """ add a new parameter """
        P = Locus(name, lower_bound, upper_bound, bitlength)
        if not isinstance(P, Locus):
            err.write("Error: added parameter ", P, "is not of class 'Locus'\n")
            sys.exit()
        else:
            self.Locuses.append(P)
    
    def num_genes(self):
        """number of genes"""
        return len(self.Locuses)
    
    def num_bits(self):
        """length in bits"""
        s = 0
        for L in self.Locuses:
            s += L.nbits
        return s
    
    def read(self, filename):
        read_params = False
        vals = {}
        with open(filename, 'r') as f:
            for s in f:
                try:
                    exec(s, {}, vals)
                except Exception as e:
                    sys.stderr.write("Error: in `%s': %s" %(s, str(e)))
        if not vals: 
            sys.stderr.write("Error: `%s` could not be read correctly!"%filename)
        for key, val in vals.items():
            if isinstance(val, float):
                self.add(key, val, val, 0)
            elif len(val) == 2:
                self.add(key, val[0][0], val[0][1], val[1])
            else:
                sys.stderr.write("Unexpected genetics: `%s = %s`" %(key, str(val)))
    
    def print(self):
        for L in self.Locuses:
            n = L.nbits
            l = L.value_from_bits(0)
            if n > 0:
                u = L.value_from_bits(2**n-1)
                print(" %20s in [ %f %f ] (%i bits)"%(L.name, l, u, L.nbits))
            else:
                print(" %20s = %f"%(L.name, l))
    
    def save(self, filename):
        """ print parameter names along with lower and upper bounds """
        with open(filename, 'w') as f:
            for L in self.Locuses:
                f.write("%s %f %f\n"%(L.name, L.lowerbound, L.upperbound))
    
    def random_genome(self):
        """ create random bit sequence of appropriate length """
        res = []
        for L in self.Locuses:
            x = random.getrandbits(L.nbits)
            #x = random.randrange(0, 1<<L.nbits)
            res.append(x)
        return res
    
    def random_select(self, mom, dad):
        """ select value from self or other randomly """
        res = []
        for m, d, L in zip(mom, dad, self.Locuses):
            x = random.getrandbits(L.nbits)
            res.append((m&x)|(d&(~x)))
        return res
    
    def print_mutations(self, guy, res):
        """ This uses ANSI escape codes to print in color on terminal """
        str = ''
        gs = self.genome_to_string(guy)
        rs = self.genome_to_string(res)
        for g, r in zip(gs, rs):
            if g == r:
                str += g
            else:
                str += '\033[1;35m'+g+'\033[0m'
        print(str)
    
    def one_mutated(self, guy):
        """ Return copy of `guy` with 1 bit flipped """
        res = []
        m = random.randrange(0, self.num_bits())
        for g, L in zip(guy, self.Locuses):
            if 0 <= m and m < L.nbits:
                g = g ^ ( 1 << m )
            m -= L.nbits
            res.append(g)
        return res

    def mutated(self, guy, cnt):
        """ Return copy of `guy` with `cnt` bits flipped on average """
        res = []
        P = float(cnt) / self.num_bits()
        for g, L in zip(guy, self.Locuses):
            N = L.nbits
            s = ['0'] * N
            for i in range(N):
                if random.random() < P:
                    s[i] = '1'
            m = int(''.join(s), 2)
            #print('mutated :', s, bin(m))
            res.append(g^m)
        #self.print_mutations(guy, res)
        return res
    
    @staticmethod
    def print_crossover(res, dad):
        """ debug printout to check the crossover function below """
        str = ''
        for i in range(len(res)):
            if res[i] == dad[i]:
                str += '-'
            else:
                str += 'X'
        print(f' {str} : {res}')
    
    def mate_crossover(self, dad, mom):
        """
        Combine `dad` and `mom` as circular sequences, using 1 crossover
        The crossover cuts between the genes, avoiding any intra-gene split
        """
        N = len(dad)
        A = random.randint(0, N)
        B = random.randint(0, N)
        res = [ i for i in dad ] # need to make a deep copy!
        if B > A:
            for i in range(A, B):
                res[i] = mom[i]
        elif B < A:
            if A < N:
                for i in range(A, N):
                    res[i] = mom[i]
            for i in range(0, B):
                res[i] = mom[i]
        #self.print_crossover(res, dad)
        return res
    
    def genome_to_values(self, guy):
        """ return a dictionary containing all the numeric values """
        res = {}
        for g, L in zip(guy, self.Locuses):
            val = L.value_from_bits(g)
            res[L.name] = val
        #print(guy, "--->", res)
        return res
    
    def genome_to_string(self, guy):
        """ return a string with 0s & 1s separated by single quotes """
        res = ''
        for g, L in zip(guy, self.Locuses):
            res += "'" + format(g, 'b').rjust(L.nbits, '0')
        return res[1:]
    
    def genome_to_hex(self, guy):
        """ return a compact representation of the genome """
        res = ''
        for g, L in zip(guy, self.Locuses):
            nb = (( L.nbits + 3 ) >> 2 )  # two hex characters per byte!
            res += format(g, 'X').rjust(nb, '0')
        return res

#-------------------------------------------------------------------------------

def load_genetics(filename):
    g = LocusSet()
    g.read(filename)
    g.print()
    return g
