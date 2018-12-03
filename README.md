Mining
======

Look at the mining/README.md file


Preprocessing the result of the mining phase
============================================

Look at preprocessing/preprocess.py


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
