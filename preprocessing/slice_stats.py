import subprocess

k = 15

for i in range(1 << k):
	slice_file = '../data/slice/{:03x}'.format(i)
	cmds = [, ]
	subprocess.run("./slice_stats --slice {f} > {f}.csv".format(f=slice_file), shell=True).check_returncode()
