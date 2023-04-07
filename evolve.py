#!/usr/bin/env python3
# A script to optimize simulation parameters via a genetic algorithm.
# Maud Formanek, 2020, Francois Nedelec 10.2021, 12.2021--03.2022; 09--10.2022
# Copyright Sainsbury Laboratory, Cambridge University, UK


"""
    This performs a parameter search using a genetic algorithm. A parameter set
    is associated to an individual for which fitness is defined by running simulations.
    The simulation results are analysed on the fly and compared to the 'target' data.
 
    This will use 'preconfig.py' to vary parameters in a template file.
    Thus `preconfig.py' must be in the working directory
    
    evolve.py relies on a module called 'arena' defining the fitness function,
    and to the executable and template files. Check for 'arena.py' in the CWD.
    You can reprogram the genetic algorithm by providing a new version of 'arena.py'
    with the same methods, changing the function calculation of fitness.
    
Syntax:

    evolve.py [RESTART] [evolve.config]
    
    The RESTART argument can be used to adjust initialization.
    RESTART can be a file ('generation.txt'), or a directory ('gen0010') created
    during a previous evolution.
    
    The parameters of the evolution are read from `evolve.config', or another
    file if specified.

Maud Formanek, 2020, Francois J. Nedelec 2021--2022
With contributions and testing by Minghua Yin, 2022
Copyright Sainsbury Laboratory, Cambridge University
"""

# Loading modules 
try:
    import os, sys, re, random, resource
    from multiprocessing import Process, Queue
    from genetics import load_genetics
    from genetools import *
    from creature import *
    import preconfig
    import arena
except ImportError as e:
    sys.stderr.write("Error loading module: %s\n"%str(e))
    sys.exit(11)

#-------------------------------------------------------------------------------

def save_generation(filename, pop):
    """ save summary of all individuals to text file """
    f = open(filename, 'w', buffering=1)
    for i in pop:
        seq = str(i)
        val = tidy_dictionary(i.values())
        f.write(f'{i.home} {i.fitness:10.6f} {seq}  {val}\n')
    f.close()


def load_generation(filename, code):
    """ load a multiple genomes from a text file """
    pop = []
    file = open(filename, 'r')
    for line in file:
        s = line.split()
        if re.match('gen[0-9]{4}', s[0]):
            f = float(s[1])
            g = s[2]
        elif s[0] == 'fitness': # old format
            f = float(s[1])
            g = s[2]
        elif s[1] == 'fitness': # older format
            f = float(s[2])
            g = s[4]
        else:
            continue
        pop.append(Creature(code, g, f))
        #print(f'loaded {g} fitness {f:10.6f}')
    file.close()
    return pop


def recover_generation(arg, code, target):
    """ load a population from previously run generation """
    pop = []
    gen = 0
    if not arg.endswith('/'):
        gen = int(arg[-4:])
        arg += '/'
    else:
        gen = int(arg[-5:-1])
    pop = load_generation(arg+'generation.txt', code)
    cnt = 0
    for i in pop:
        i.update_home(arg, cnt)
        cnt += 1
        nam = os.path.join(i.home, results_file)
        if 1:
            all = load_results(nam)
            if all:
                print(f'found {len(all)} dataset for {nam}', end='')
                i.fitness = arena.calculate_fitness(all, target)
                print(f': {i.home}.fitness = {i.fitness}')
            else:
                i.fitness = math.nan
        else:
            print(f'no data for {i.home}')
    return pop, gen


def random_generation(size, code):
    """ return a population of given size with random genomes """
    pop = []
    for i in range(size):
        g = code.random_genome()
        pop.append(Creature(code, g, i))
        #print(g)
    return pop

#-------------------------------------------------------------------------------

def cpu_sec():
    """ return seconds spent by this and child processes in user mode """
    return resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime


def bars(val, cnt):
    """ generate a string of size 'cnt' representing `val`"""
    if val > 0:
        c = int(val/100)
        x = int((val-100*c)/10)
        n = int(val-100*c-10*x)
        c = min(c, cnt)
        x = min(x, cnt-c)
        n = min(n, cnt-x-c)
        txt = 'C'*c + 'X'*x + '|'*n
    else:
        txt = "%.3f"%val
    return txt.ljust(cnt)


def find_soul(heaven, bug):
    """ recover data from previous generation if 'bug' already lived"""
    res = ''
    new = bug.home
    seq = bug.sequence()
    cnt = 0
    while seq in heaven:
        old = heaven[seq]
        # using first part of name to identify the generation:
        if old.split('/')[0] == new.split('/')[0]:
            #print("sequence %s already exists in same generation" %seq)
            bug.mutate(1)
            seq = bug.sequence()
            #print(" mutated -> %s" %seq)
            cnt += 1
            if cnt > 1024:
                break
        else:
            res = old
            break
    heaven[seq] = new
    return res


def run_simulation(bug):
    """ Run all config files found in directory 'bug.home' and return fitness """
    # Attention: in multiprocessing, any modification to argument 'bug' is lost!
    # and global variables values are not those in the parent thread!
    cpu = cpu_sec()
    template = os.path.abspath(arena.template)
    # generate config files in bug's directory:
    files = preconfig.parse(template, bug.dic, 1, bug.home)
    #files = find_files('.cym')
    #print(bug.home+" -------> "+' '.join(files))
    # run simulations and get textual output:
    all, old = run_many(arena.simex, arena.simout, bug.home, files)
    # fitness from fresh data only:
    fit = arena.calculate_fitness(all, bug.target)
    if math.isnan(fit):
        sys.stderr.write(f'Error: invalid fitness value for {bug.home}!\n')
    verbose = ( old ) or ( int(os.path.basename(bug.home)) < 16 )
    if verbose:
        cpu = cpu_sec() - cpu
        print(f'{bars(fit*10, 16)} CPU {cpu:6.1f} {bug.home} {bug!s}', end='')
    if old:
        fff = fit
        ggg = arena.calculate_fitness(old, bug.target)
        # include data computed in previous generations:
        all.extend(old)
        # set fitness from all available data:
        fit = arena.calculate_fitness(all, bug.target)
        if verbose:
            print(f' fit {ggg:10.4f} x{len(old)} {fff:+10.4f}  ---> {fit:10.4f}')
    else:
        if verbose:
            print(f' fit {fit:10.4f} ({len(all)})')
    return fit


def run_simulation_trace(bug):
    """ Simulate one Creature and return its fitness"""
    # Attention: in multiprocessing, any modification to argument 'bug' is lost!
    # and global variables values are not those in the parent thread!
    try:
        return run_simulation(bug)
    except Exception:
        import traceback
        traceback.print_exc()
        sys.exit(7)

#-------------------------------------------------------------------------------

def worker(food, poop):
    """
    This is executed by a child processes in multiprocessing mode.
    run simulation taking argument from `feed`, returning results in `poop`
    """
    while True:
        try:
            index, bug = food.get(True, 1)
        except:
            break
        res = run_simulation(bug)
        poop.put((index, res))


def run_parallel(pop, njobs):
    """
    Run all 'bug' in the group using multiple threads, and set bug.fitness
    """
    food = Queue()
    poop = Queue()
    # fill queue
    for n, i in enumerate(pop):
        food.put((n, i), False)
    jobs = []
    # create threads, and start them running 'worker(food,poop)'
    for n in range(njobs):
        j = Process(target=worker, args=(food,poop))
        jobs.append(j)
        j.start()
    # wait for completion of all jobs:
    for j in jobs:
        j.join()
        j.close()
    # attribute calculated fitness to local data:
    while True:
        try:
            i, f = poop.get_nowait()
            pop[i].fitness = f
            #print(f'   ({i}).fitness = {f}')
        except:
            break


def run_all_simulations(pop, njobs: int):
    """
    Run all 'bug' in the group, and set bug.fitness
    """
    if njobs > 1:
        # ...in parallel on multiple threads:
        run_parallel(pop, njobs)
    else:
        # ...sequentially:
        for i in pop:
            i.fitness = run_simulation(i)


def average_fitness(pop):
    """ Return average, variance and number of datapoint """
    cnt = 0
    avg = 0
    var = 0
    for i in pop:
        if not math.isnan(i.fitness):
            cnt += 1
            avg += i.fitness
            var += i.fitness * i.fitness
        else:
            sys.stderr.write(f'Error: missing fitness data for {i.home}!\n')
            i.fitness = 0 #sys.exit(6)
    if cnt > 0:
        avg /= cnt
        var = var - cnt * avg * avg
    else:
        avg = math.nan
    if cnt > 1:
        var = var / float(cnt-1);
    else:
        var = math.nan
    return avg, var, cnt

#-------------------------------------------------------------------------------

def main(args):
    """
    By default, this reads parameters from 'evolve.config',
    but another file, ending by .config, can be specified on the command line.
    """
    config = 'evolve.config'
    genstart = ''
    for arg in args:
        if os.path.isdir(arg):
            genstart = arg
        elif os.path.isfile(arg):
            if arg.endswith('.config'):
                config = arg
            elif arg.endswith('.txt'):
                genstart = arg
        else:
            sys.stderr.write(f'Error: unexpected argument `{arg}`\n')
            sys.exit(1)
    """
    read parameters from config file
    """
    pam = read_config(config)
    #print(pam)
    try:
        # extract some parameters out
        population_size = pam.pop('population_size')
        generation_max = pam.pop('generation_max')
        njobs = pam.pop('njobs')
        FIT_DIFF_MAX = pam.pop('FIT_DIFF_MAX')
        FIT_MAX = pam.pop('FIT_MAX')
        elitism = pam.pop('elitism')
        crossover = pam.pop('crossover')
        mutation = pam.pop('mutation')
        mutation_cnt = pam.pop('mutation_cnt')
        if not isinstance(mutation_cnt, list):
            sys.stderr.write("Error: `mutation_cnt` should be a list of 2 values\n")
        resurrect = pam.pop('resurrect')
        n_plus = pam.pop('n_plus')
    except KeyError as e:
        sys.stderr.write(f'Error: parameter {e!s} must be defined in `{config}`\n')
        sys.exit(2)
    # Remaining parameters will be passed on to preconfig, generating multiple files
    for k, v in pam.items():
        print(f'Parameter: {k} = {v}')
    # check probabilities:
    if elitism + crossover + mutation > 1:
        sys.stderr.write("Error: elitism, crossover, mutation must add up to 1!\n")
        sys.exit(3)
    # check executable is there:
    if not os.access(arena.simex, os.X_OK):
        sys.stderr.write(f'Error: could not find executable `{arena.simex}`\n')
        sys.exit(4)
    # initialize variables:
    random.seed()
    heaven = {}
    gen_idx = 0
    target = arena.load_target()
    code = load_genetics(arena.genetic_config)
    """
        initialize population from 'genstart' if specified
    """
    if os.path.isdir(genstart):
        # restart from previous generation folder:
        popolo, gen_idx = recover_generation(genstart, code, target)
        # run simulations that have no fitness:
        todo = [ i for i in popolo if math.isnan(i.fitness) ]
        for i in todo:
            i.target = target
            i.update_dictionary(pam)
        run_all_simulations(todo, njobs)
        save_generation(f'gen{gen_idx:04}/generation.txt', popolo)
        print(f'loaded {genstart} with {len(popolo)} genomes')
        gen_idx += 1
    elif os.path.isfile(genstart):
        # restart from genomes specified in a file:
        popolo = load_generation(genstart, code)
        print(f'loaded {len(popolo)} genomes from {genstart}')
    else:
        popolo = random_generation(population_size, code)
        print(f'initialized with {len(popolo)} random genomes')
    fitness_old = popolo[0].fitness
    """
        start evolution loop
    """
    while gen_idx < generation_max:
        cpu = cpu_sec()
        # create directory for new generation:
        genpath = f'gen{gen_idx:04}/'
        make_new_dir(genpath)
        gen_idx += 1
        # copy template and genetics file:
        copy_file(arena.template, genpath)
        copy_file(arena.genetic_config, genpath)
        # prepare population:
        cnt = 0
        for i in popolo:
            i.target = target
            i.fitness = math.nan
            i.update_home(genpath, cnt)
            cnt += 1
            # Attention: find_soul() can mutate 'i'!
            old = find_soul(heaven, i)
            i.update_dictionary(pam)
            if old and resurrect:
                copy_files(old, i.home)
                #move_files(old, i.home)
            else:
                os.mkdir(i.home)
        # save generation without fitness values:
        save_generation(genpath+"generation.txt", popolo)
        # optional verification:
        # report_duplicates(popolo)
        # Calculate fitness values for all...
        run_all_simulations(popolo, njobs)
        # check that all fitness values are set and valid
        avg, var, cnt = average_fitness(popolo)
        print(f'{genpath} average fitness for {cnt} bugs : {avg:.3f} +/- {var:.3f}')
        if cnt < 1:
            print("Fatal error: no fitness in whole generation!")
            return
        # sort population by decreasing fitness:
        popolo = sorted(popolo, key = lambda x:-x.fitness)
        # best is now first in list:
        best = popolo[0]
        print(f'Best {best.fitness:.4f} {best.home}', end=': ')
        print(tidy_dictionary(best.values()))
        # save generation with fitness values:
        save_generation(genpath+"generation.txt", popolo)
        # save tested genomes:
        save_dictionary(genpath+"heaven.txt", heaven)
        print(f'CPU {cpu_sec()-cpu:8.1f}', end=' ')
        print(f'{genpath}heaven.txt has {len(heaven)} souls\n')
        """
            Check for convergence, given the best fitness reached in this generation
        """
        if abs(best.fitness - fitness_old) < FIT_DIFF_MAX:
            print("FIT_DIFF_MAX is achieved")
            break
        if best.fitness > FIT_MAX:
            print("FIT_MAX is exceeded")
            break
        if gen_idx >= generation_max:
            break
        """
            Creation of a new generation
        """
        fitness_old = best.fitness
        #elders = [ i for i in popolo if i.fitness > 0 ]
        elders = popolo
        popolo = []

        S, P, Q = partition(len(elders), elitism, crossover, mutation)
        #print(f'elitism {S}, crossover {P}, mutation {Q} ')
        """
            Elitism: transfer to the next generation without change
        """
        for i in elders[:S]:
            popolo.append(i)
            #print(f'elite: {i.home}  {str(i)}  fitness {i.fitness:10.4f}")
        print(f'{S} elite', end='')
        """
            Crossover: mating of two parents to produce one child
        """
        # selection of parent's indices:
        parents = fitness_rank_selection(len(elders), 2*P, n_plus, 2.0-n_plus)
        # check that we have at least 2 different possible parents:
        if len(parents) > 1:
            print(f', crossover from {len(parents)>>1} pairs', end='')
            for i in range(P):
                m = parents[2*i]
                d = parents[2*i+1]
                # dad and mom should be different, so take anyone else in last resort:
                while d == m:
                    d = random.randrange(0, len(elders))
                #print(f' {min(d,m)}:{max(d,m)}', end='')
                mom = elders[m]
                #child = mom.mate_uniform(elders[d])
                child = mom.mate_crossover(elders[d])
                popolo.append(Creature(mom.code, child))
            #print("")
        """
            Mutation: apply mutation rate to the rest of the population
        """
        parents = fitness_rank_selection(len(elders), Q, n_plus, 2.0-n_plus)
        #parents = fitness_proportional_selection(elders, Q)
        if parents:
            print(f', mutations for {len(parents)}', end='')
            for i in parents:
                bug = elders[i]
                # use different rate of mutation above/below the mean:
                mut = mutation_cnt[bug.fitness < avg]
                # print(f'mutation({bug.fitness}) = {mut}')
                popolo.append(Creature(bug.code, bug.mutated(mut)))
        """
           Complete population with random genomes
        """
        cnt = 0
        while len(popolo) < population_size:
            g = code.random_genome()
            popolo.append(Creature(code, g))
            cnt += 1
        # final message:
        if cnt > 0:
            print(f', {cnt} new random genomes')
        else:
            print('')
    """
        Print final results at termination
    """
    print(f'Best fitness: {best.fitness}')
    print(best.values())


#------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) > 2:
        print(__doc__)
    elif len(sys.argv) == 2 and sys.argv[1]=='help':
        print(__doc__)
    else:
        main(sys.argv[1:])
