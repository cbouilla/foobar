import sys
import glob
import subprocess
import os
import os.path
import math

import argparse
parser = argparse.ArgumentParser(description="This script runs all the preprocessing")
parser.add_argument("--partitioning-bits", help="number of bits to split", type=int, required=True, dest="k")
parser.add_argument("--stats", help="display stats", action="store_true")
parser.add_argument("--dry-run", help="print commands but don't run them", action="store_true")
parser.add_argument("--verbose", help="display more information", action="store_true")
parser.add_argument("--cores", help="number of cores of this machine", type=int, default=1)
#parser.add_argument("--check", help="run check programs", action="store_true")
parser.add_argument("--slice", help="run the slicer [for joux's algo]", action="store_true")

args = parser.parse_args()

print("Running with k = {}".format(args.k))
print()


L1_cache_size = 16384

#KINDS = {'foo': 0, 'bar': 1, 'foobar': 2}
#PREIMAGE_DIR = '../data/preimages'
#PREIMAGE_DIR = '../foobar'
#DICT_DIR = '../data/dict'
HASH_DIR = '../data/hash_c'
SLICE_DIR = '../data/slice_test'
#SPLITTER = './splitter'
#DICT_CHECKER = './dict_checker'
#SORTER = './sorter'
#MERGER = './merger'
#HASH_CHECKER = './hash_checker'
SLICER = './slicer_test'

SLICE_L = 20

#do_split = True
#check_split = False
#do_sort = True
#check_sort = False
#do_merge = True
#check_merge = False

#n_preimages = {}
n_dict = {}
n_hash = {}





def slicing():
    """
    compute the slice for foobar hash files.
    """
    for i in range(1 << args.k):
        input_file = '{}/foobar.{:03x}'.format(HASH_DIR, i)
        output_file = '{}/{:03x}'.format(SLICE_DIR, i)
        if not os.path.exists(input_file):
            continue
        if os.path.exists(output_file):
            input_mtime = os.stat(input_file).st_mtime
            output_mtime = os.stat(output_file).st_mtime
            if input_mtime < output_mtime:
                if args.verbose:
                    print('# skipping {} [already sliced]'.format(output_file))
                continue
        cmds = [SLICER, '--l', str(SLICE_L), '--target-dir', SLICE_DIR, input_file]
        if args.verbose:
            print(" ".join(cmds))
        if not args.dry_run:
            subprocess.run(cmds, stdout=subprocess.DEVNULL).check_returncode()



if args.slice:
    print("4. Slicing (hash files -> slice files)")
    slicing()
