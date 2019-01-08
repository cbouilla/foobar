To try the quadratic algorithm :

1) Generate an instance of a given size N with known solutions

	# preprocessing/forger N

N = 128000 is reasonable. The 3 files:
- bar.000  
- foo.000
- foobar.000 
are created in the current directory.


2) Solve it

	# ./quad_standalone --partitioning-bits 0 --hash-dir HASHPATH

Set HASHPATH must be the path where the 3 hashfiles are located.
This should report the solutions announced by the forger.