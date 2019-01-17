To try the iterated Joux algorithm :

1) Generate an instance of a given size N with known solutions

	# preprocessing/forger N

N = 128000 is reasonable. The 3 files:
- bar.000  
- foo.000
- foobar.000 
are created in the current directory.

1.5) Generate the slices (in the current directory).

	# slicer foobar.000 --target-dir . --l 19


2) Solve it

	# ./joux_standalone --partitioning-bits 0 --hash-dir HASHPATH --slice-dir SLICEPATH

Set HASHPATH must be the path where the 3 hashfiles are located.
This should report the solutions announced by the forger.