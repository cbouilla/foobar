import glob
import subprocess
import multiprocessing
import itertools
import os.path
import math

"""
This script runs all the preprocessing.
"""


# on how many bits are the preimage files split. This has no influence on the resulting computation.
# small values yields too big dictionnaries. high values yields too many dictionnaries.
# With approx 1Gb premage files, 2 is a reasonable choice.
split_bits = 2

cores = 4



KINDS = {'foo': 0, 'bar': 1, 'foobar': 2}
PREIMAGE_DIR = '../data/preimages'
DICT_DIR = '../data/dict'
HASH_DIR = '../data/hash'
SPLITTER = './splitter'
DICT_CHECKER = './dict_checker'
SORTER = './sorter'
MERGER = './merger'
HASH_CHECKER = './hash_checker'

do_split = False
check_split = False
do_sort = False
check_sort = False
do_merge = False
check_merge = True

n_preimages = {}
n_dict = {}
n_hash = {}


def preimage_stats(verbose=True):
    for kind, kind_id in KINDS.items():
        total = 0
        for file in glob.glob('{}/{}.*.work'.format(PREIMAGE_DIR, kind)):
            total += os.path.getsize(file)
        n_preimages[kind] = total // 12
    if verbose:
        print('Preimage stats')
        print('==============')
        for kind in KINDS:
            total = n_preimages[kind]
            print("|{:^6}| = {} (2^{:.2f}, {:.1f}M preimages)".format(kind, total, math.log(total, 2), total / 1024**2))

def dict_stats(verbose=True):
    for kind, kind_id in KINDS.items():
        total = 0
        for file in glob.glob('{}/{}.*'.format(DICT_DIR, kind)):
            total += os.path.getsize(file)
        n_dict[kind] = total // 20
    if verbose:
        print('Dictionary stats')
        print('================')
        for kind in KINDS:
            total = n_dict[kind]
            pre = n_preimages[kind]
            dup = 100 * (pre - total) / pre
            print("|{:^6}| = {} (2^{:.2f}, {:.1f}M dict entries) [invalid rate={:.3f}%]".format(kind, total, math.log(total, 2), total / 1024**2, dup))

def hash_stats(verbose=True):
    for kind, kind_id in KINDS.items():
        total = 0
        for file in glob.glob('{}/{}.*'.format(HASH_DIR, kind)):
            total += os.path.getsize(file)
        n_hash[kind] = total // 8
    if verbose:
        print('Hash files stats')
        print('================')
        for kind in KINDS:
            total = n_hash[kind]
            d = n_dict[kind]
            dup = 100 * (d - total) / d
            print("|{:^6}| = {} (2^{:.2f}, {:.1f}M hashes) [duplicate rate={:.3f}%]".format(kind, total, math.log(total, 2), total / 1024**2, dup))


def splitting():
    """
    split all preimage files using [[split_bits]]
    """
    print("1. Splitting (preimage --> dictionaries), [split_bits={}]".format(split_bits))
    mpi_n_process = 1 + (1 << split_bits) + 2 * cores
    for kind, kind_id in KINDS.items():
        files = glob.glob('{}/{}.*.work'.format(PREIMAGE_DIR, kind))
        for preimage in files:
            args = ['mpirun', '-np', mpi_n_process, SPLITTER, '--kind', kind_id, 
                    '--bits', split_bits, '--output-dir', DICT_DIR, preimage]
            print("    -> {}".format(" ".join(map(str, args))))
            subprocess.run(map(str, args)).check_returncode()

def _do_sort(job):
    _, file = job
    args = [SORTER, file]
    print("    -> {}".format(" ".join(args)))
    subprocess.run(args).check_returncode()

def sorting():
    """
    Sort all dictionnaries.
    """
    print("2. Sorting (dictionaries -> dictionaries)")
    jobs = []
    for kind, kind_id in KINDS.items():
        for file in glob.glob('{}/{}.*'.format(DICT_DIR, kind)):
            jobs.append((kind_id, file))
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(_do_sort, jobs, chunksize=1)

def _do_check_dict(job):
    kind, file = job
    args = [DICT_CHECKER, '--kind', str(kind), file]
    print("    -> {}".format(" ".join(args)))
    subprocess.run(args).check_returncode()

def check_dict():
    """
    Verify all dictionnaries.
    """
    print("X. Checking dictionaries")
    jobs = []
    for kind, kind_id in KINDS.items():
        for file in glob.glob('{}/{}.*'.format(DICT_DIR, kind)):
            jobs.append((kind_id, file))
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(_do_check_dict, jobs, chunksize=1)

            
def _do_merge(job):
    kind, bucket = job
    input_files = glob.glob('{}/{}.*.dict.{}'.format(DICT_DIR, kind, bucket))
    output_file = '{}/{}.hash.{}'.format(HASH_DIR, kind, bucket)
    #print("    -> merging {}/{}".format(kind, bucket))
    args = [MERGER, '--output', output_file] + input_files
    print("    -> {}".format(" ".join(args)))
    subprocess.run(args).check_returncode()            


def merging():
    """
    merge all (sorted) dictionnaries.
    """
    print("3. Merging (dictionaries -> hash files)")
    jobs = itertools.product(KINDS, range(1 << split_bits))
    with multiprocessing.Pool(processes=1) as pool:
        pool.map(_do_merge, jobs, chunksize=1)

def _do_check_hash(job):
    _, file = job
    args = [HASH_CHECKER, file]
    print("    -> {}".format(" ".join(args)))
    subprocess.run(args).check_returncode()

def check_hash():
    """
    Verify all dictionnaries.
    """
    print("Y. Checking hashes")
    jobs = []
    for kind, kind_id in KINDS.items():
        for file in glob.glob('{}/{}.*'.format(HASH_DIR, kind)):
            jobs.append((kind_id, file))
    with multiprocessing.Pool(processes=cores) as pool:
        pool.map(_do_check_hash, jobs, chunksize=1)

        

preimage_stats()
if do_split:
    splitting()
if check_split:
    check_dict()
dict_stats()
if do_sort:
    sorting()
if check_sort:
    check_dict()
if do_merge:
    merging()
hash_stats()
if check_merge:
    check_hash()


print("DONE")