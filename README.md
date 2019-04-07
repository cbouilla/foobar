Mining
======

Look at the mining/README.md file



Preprocessing the result of the mining phase
============================================

sudo apt-get update
sudo apt-get install git build-essential noweb mpi-default-dev s3fs libm4ri-dev (?)

git clone https://github.com/cbouilla/foobar.git
mkdir foobar/foobar

s3fs 3sum-foobar foobar -o ro,uid=1000,umask=0004,public_bucket=1
sudo sh -c "echo 800000 > /proc/sys/fs/file-max"

[new shell process]

ulimit -n 100000

mkdir data/dict
mkdir data/hash
mkdir data/slice
cd preprocessing
screen python3 preprocess.py --cores 4 --partitioning-bits 15 --verbose


bad data
========

foo.8 vérollé à 92%

---> couper à l'offset 0x62fa9ffc
---> reprendre à l'offset 0x62faa000

Solving the instance of the problem
===================================

Look at solving/



TODO
====

+ données auxiliaires
++ quad

---> constructions de cuckoo hash à l'avance pour chaque chunk, avec biais pour H1 (algo hongrois ?).

+ ijoux

---> matrices

+ quadratic on GPU
+ quadratic on KNL
+ ijoux on BG/Q

Budget 30-40x plus grand / item que pour quad.

++ hash join ? 
+++ 1 step partitioning ? ---> impose petites tailles 
+++ 2 step partitiining ? ---> permet grandes tailles
++ sort join ?
