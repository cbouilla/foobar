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
split_bits = 10
dry_run = False

# number of CPU cores of this machine.
cores = 4

L1_cache_size = 16384


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
check_merge = False

n_preimages = {}
n_dict = {}
n_hash = {}

def ensure_dirs():
    for i in range(1 << split_bits):
        d = '{}/{:03x}'.format(DICT_DIR, i)
        if not os.path.exists(d):
            os.mkdir(d)

def preimage_stats(verbose=True):
    for kind, kind_id in KINDS.items():
        total = 0
        for file in glob.glob('{}/{}.*'.format(PREIMAGE_DIR, kind)):
            total += os.path.getsize(file)
        n_preimages[kind] = total // 12
    if verbose:
        print('Preimage stats')
        print('==============')
        for kind in KINDS:
            total = n_preimages[kind]
            print("|{:^6}| = {} (2^{:.2f}, {:.1f}M preimages)".format(kind, total, math.log(total, 2), total / 1024**2))
        print()

def dict_stats(verbose=True):
    for kind in KINDS:
        total = 0
        for file in glob.glob('{}/*/{}.*.unsorted'.format(DICT_DIR, kind)):
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
        print()

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
        n_entries = 1
        for kind in KINDS:
        	n_entries *= n_hash[kind]
        print("Est # (64+k)-bit solutions = {:.1f}".format(n_entries / 2**(64 + split_bits)))
        print()

def splitting():
    """
    split all preimage files using [[split_bits]]
    """
    print("1. Splitting (preimage --> dictionaries), [split_bits={}]".format(split_bits))
    mpi_n_process = 2 + 2 * cores
    for kind, kind_id in KINDS.items():
        files = glob.glob('{}/{}*'.format(PREIMAGE_DIR, kind))
        for preimage in files:
            args = ['mpirun', '-np', mpi_n_process, SPLITTER, 
                    '--partitioning-bits', split_bits, '--output-dir', DICT_DIR, preimage]
            print("    -> {}".format(" ".join(map(str, args))))
            if not dry_run:
            	subprocess.run(map(str, args)).check_returncode()

def _do_sort(file):
    args = [SORTER, file]
    print("    -> {}".format(" ".join(args)))
    sorted_filename = file[:-8] + 'sorted'
    if os.path.exists(sorted_filename) and os.path.getsize(file) == os.path.getsize(sorted_filename):
        print("        [SKIP]")
        return
    if not dry_run:
        subprocess.run(args).check_returncode()

def sorting():
    """
    Sort all dictionnaries.
    """
    print("2. Sorting (dictionaries -> dictionaries)")
    jobs = []
    for i in range(1 << split_bits):
        for file in glob.glob('{}/{:03x}/*.unsorted'.format(DICT_DIR, i)):
            _do_sort(file)

def _do_check_dict(file):
    args = [DICT_CHECKER, '--partitioning-bits', str(split_bits), file]
    print("    -> {}".format(" ".join(args)))
    if not dry_run:
        subprocess.run(args).check_returncode()

def check_dict():
    """
    Verify all dictionnaries.
    """
    print("X. Checking dictionaries")
    jobs = []
    for i in range(1 << split_bits):
        for file in glob.glob('{}/{:03x}/*'.format(DICT_DIR, i)):
            _do_check_dict(file)
    


def merged_hashfile_name(kind, i):
    return '{}/{}.{:03x}'.format(HASH_DIR, kind, i)

def _do_merge(job):
    kind, i = job
    input_files = glob.glob('{}/{:03x}/{}.*.sorted'.format(DICT_DIR, i, kind))
    output_file = merged_hashfile_name(kind, i)
    args = [MERGER, '--output', output_file] + input_files
    print("    -> {}".format(" ".join(args)))
    if not dry_run:
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
    if not dry_run:
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


ensure_dirs()
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