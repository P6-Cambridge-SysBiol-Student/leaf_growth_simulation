# Handy tools for a genetic algorithm.
# Maud Formanek, 2021, FJ Nedelec, 2021--2022, Copyright Cambridge University

try:
    import os, sys, shutil, re, random
    from subprocess import Popen, PIPE
except ImportError as e:
    sys.stderr.write("Error loading modul: %s\n"%str(e))
    sys.exit()

# name of file used to save consolidated output of simulation runs
results_file = 'all_results.txt'

#-------------------------------------------------------------------------------

def fitness_rank_selection(top, cnt, n_plus, n_minus):
    """
    Selection based on rank: returning 'cnt' integers in [0, top-1]
    'n_plus' is the probability of selecting the topmost individual,
    and 'n_minus' the probability of selecting the last ranking one.
    Having 'm_minus = 2 - n_plus' results in 'sum == top'
    """
    if top < 2:
        return []
    sel = []
    sum = 0
    beta = [0] * top
    alpha = ( n_minus - n_plus ) / ( top - 1 )
    # calculate cummulative distribution:
    for j in range(top):
        sum += n_plus + alpha * j
        beta[j] = sum
    #print("rank_selection ", beta, sum, top)
    for _ in range(cnt):
        t = random.random() * sum
        s = 0
        for j in range(top):
            if t < beta[j]:
                s = j
                break
        sel.append(s)
    return sel


def fitness_proportional_selection(pop, cnt):
    """
        Selection of parents based on relative fitness value
        Returns 'cnt' integers in [0, len(pop)-1]
    """
    sel = []
    sum = 0
    beta = [0] * len(pop)
    for i, p in enumerate(pop):
        sum += max(p.fitness, 0)
        beta[i] = sum
    #print("fitness_selection ", beta)
    for _ in range(cnt):
        t = random.random() * sum
        s = 0
        for j, b in enumerate(beta):
            if t < b:
                s = j
                break
        sel.append(s)
    return sel


def partition(TOT, pA, pB, pC):
    """
    Return 3 random integers with sum = TOT * ( pA + pB + pC ),
    distributed according to the probabilities pA, pB and pC
    """
    A = 0
    B = 0
    C = 0
    pAB = pA + pB
    pABC = pA + pB + pC
    for _ in range(TOT):
        x = random.random()
        if x < pA:
            A += 1
        elif x < pAB:
            B += 1
        elif x < pABC:
            C += 1
    return A, B, C


def report_duplicates(pop):
    """ find duplicales in 'pop' """
    seen = set()
    for i in pop:
        g = i.sequence()
        if g in seen:
            print("duplicate genome: " + g)
        else:
            seen.add(g)

#-------------------------------------------------------------------------------

def move_files(old, new):
    """
        Move files between generations
    """
    try:
        shutil.move(old, new)
        #print(old, " ---move---> ", new)
        return 0
    except Exception as e:
        print("move_files failed: " + str(e))
        return 1


def copy_files(old, new):
    """
        Copy directorty and contained files between generations
    """
    #shutil.copytree(old[:-1], self.path[:-1])
    try:
        shutil.copytree(old, new)
        #print(old, " ---copy---> ", new)
        return 0
    except Exception as e:
        print("copy_files failed: " + str(e))
        return 1

#-------------------------------------------------------------------------------

def copy_file(file, path):
    """ make a copy of `file` into directory `path` """
    try:
        if os.path.isdir(path):
            shutil.copyfile(file, os.path.join(path,file))
        else:
            shutil.copyfile(file, path)
    except Exception as e:
        sys.stderr.write(f'Error: could not copy `{file}`: {e!s}\n')
        sys.exit(9)


def find_files(ext):
    """ return list of files with given extension in current working directory """
    files = []
    with os.scandir() as it:
        for e in it:
            if e.name.endswith(ext) and e.is_file():
                files.append(e.name)
    return sorted(files)


def make_new_dir(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        sys.stderr.write(f'Error: directory `{path}` exists already\n')
        sys.exit(5)


def make_dir(path):
    """ Make directory if it does not exist"""
    if not os.path.exists(path):
        os.mkdir(path)


def make_directory(root, cnt):
    dir = root+'%04i' % cnt
    while os.path.exists(dir):
        cnt += 1
        dir = root+'%04i' % cnt
    try:
        os.mkdir(dir)
    except:
        import tempfile
        dir = tempfile.mkdtemp('', 'root', '.')
    return (dir, cnt)


#-------------------------------------------------------------------------------

def save_dictionary(file, dic):
    """write dictionary from file"""
    with open(file, 'w') as f:
        for k, v in dic.items():
            f.write("%s %s\n" %(k, str(v)))


def load_dictionary(file):
    """read values from file, returning a dictionary"""
    res = {}
    f = open(file, "r")
    for s in f:
        if not s or s[0] == '%' or s[0] == '#':
            continue
        k, v = s.split()
        res[k] = v
    f.close()
    return res


def tidy_dictionary(dic):
    """return string representing dictionary's values"""
    res = ''
    for k, v in dic.items():
        try:
            res += f'; {k}="{len(v)} values"'
        except:
            res += f'; {k}={v:.4f}'
    return res[2:]


def read_values(file, unwanted):
    """load parameter values as a dictionary from file"""
    res = {}
    with open(file, 'r') as f:
        for s in f:
            if not s.strip().startswith('%'):
                exec(s, {}, res)
    for u in unwanted:
        if u in res:
            res.pop(u)
    return res


def read_config(file):
    """ read file containing definitions, eg `x=0.3` """
    res = {}
    with open(file, 'r') as f:
        for s in f:
            if not s.strip().startswith('%'):
                for e in s.split(';'):
                    exec(e.strip(), {}, res)
    return res

#-------------------------------------------------------------------------------

def save_results(name, data, txt=''):
    """ add results to file """
    with open(name, 'a') as f:
        if txt:
            f.write(txt)
        for x in data:
            f.write(x)


def load_results(name):
    """ return list of results: one string for each run """
    out = []
    res = ''
    try:
        f = open(name, 'r')
        for s in f:
            if re.match('%gen[0-9]{4}', s): #s.startswith('%gen'):
                if res:
                    out.append(res)
                res = ''
            else:
                res += s
        f.close()
        #print("loaded %i lines of old data"%len(res))
    except FileNotFoundError:
        pass
    if res:
        out.append(res)
    return out


def load_file(name):
    """ return file content as a string """
    res = ''
    try:
        with open(name, 'r') as f:
            for s in f:
                res += s
    except FileNotFoundError:
        pass
    return res

#-------------------------------------------------------------------------------

def start_job(simex, conf, path):
    """
    Start simulation in current working directory as subprocess
    """
    # Option 'cwd=path' sets the current working directory of the subprocess
    # Cytosim's option '-' triggers the silent mode, reducing its output
    sub = Popen([simex, '', os.path.basename(conf)], stdout=PIPE, stderr=PIPE, cwd=path)
    #print(f'started process {sub.pid} in `{os.getcwd()}`')
    return sub


def wait_job(sub, simout, path):
    """
    wait for subprocess to complete and return its output
    """
    try:
        out, err = sub.communicate()
        if err:
            sys.stderr.write(err.decode().partition('\n')[0])
        if sub.returncode:
            sys.stderr.write(f'exitcode {sub.returncode} ')
        #print(f'collected process {sub.pid} in `{os.getcwd()}`')
        if simout == 'stdout':
            return out.decode()
        else:
            return load_file(os.path.join(path, simout))
    except Exception as e:
        sys.stderr.write(f'Failed: {e!s}\n')
    return ''


def run_many(simex, simout, path, files):
    """
    Run all simulations sequentially within directory `path`
    """
    all = []
    tmp = ''
    for conf in files:
        sub = start_job(simex, conf, path)
        res = wait_job(sub, simout, path)
        if res:
            all.append(res)
            #print("  %s/%s : %i bytes" %(path, conf, len(res)))
            if tmp:
                os.remove(tmp)
            tmp = conf
        else:
            sys.stderr.write(f'[{path}/{conf}]')
    # load data from previous generations:
    nam = os.path.join(path, results_file)
    old = load_results(nam)
    if all:
        # save fresh data to file:
        save_results(nam, all, '%'+path+'\n')
    else:
        sys.stderr.write(f'Error: no result found in `{path}`\n')
    return all, old


def run_many_separate(simex, root, files):
    """
    Runs all simulations sequentially in separate working directory
    It is assumed that current directory is set to 'root' already
    """
    all = []
    tmp = []
    """ run each simulation in a sub-directory """
    for i, f in enumerate(files):
        path = "run%04i" % i
        os.mkdir(path)
        [conf, ext] = os.path.splitext(f)
        conf = 'config'+ext
        os.rename(f, os.path.join(path,conf))
        sub = start_job(simex, conf, path)
        res = wait_job(sub, simout, path)
        if res:
            all.append(res)
            #print("  %s/%s : %i bytes" %(root, path, len(res)))
            tmp.append(path)
        else:
            sys.stderr.write(f'[{root}/{conf}]')
    # load data from previous generations:
    nam = os.path.join(path, results_file)
    old = load_results(nam)
    if all:
        # save fresh data to file:
        save_results(nam, all, '%'+root+'\n')
        # delete temporary directories:
        for p in tmp:
            shutil.rmtree(p)
    else:
        sys.stderr.write(f'Error: runs in `{root}` all failed!\n')
        # grab all data into a single file:
        #sys = os.system("cd %s && cat run??/results.txt > all_results.txt"%root)
        #res = load_results(nam)
    return all, old

