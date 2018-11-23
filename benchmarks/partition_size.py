from matplotlib import pyplot as plt
import numpy as np
import statistics

def load_file(fname):
    A = {}
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            k, *values = line.strip().split(";")
            k = int(k)
            values = list(map(float, values[:-1]))
            m = min(values[1:])
            M = max(values[1:])
            med = statistics.median(values[1:])
            A[k] = (med, m, M)
    return A


def process_file(filename, cpu_name):
    A_ = load_file(filename + '.txt')
    
    K = max(A_)

    A = np.zeros(shape=(K+1, 3))
    for k, (a,b,c) in A_.items():
        A[k] = (a, b, c)
    
    fig = plt.figure(figsize=(16, 9))
    
    plt.title("Partitioning Efficiency (T=max, p=8)".format(cpu_name))
    plt.grid(zorder=-10)
    
    ax = fig.add_subplot(111)
    ax.set_xlabel("Input size ($K item$)")
    ax.set_ylabel('G item/s')
    
    x = np.arange(0, K+1)
    ax.errorbar(x, A[:, 0], yerr=A[:, 1:2])
    
    #ax.set_xticks(np.arange(1, P+1))
    #ax.legend(loc='lower left',  bbox_to_anchor=(-0.025, -0.35), ncol=2)
    plt.savefig("partition_size_" + filename + ".pdf")
    #plt.show()


process_file('hpac_size', 'Intel Xeon E5-4620 @ 2.20GHz')
process_file('zen2_size', 'AMD EPYC 7301 @ 2.20GHz')
process_file('ppti_size', 'Intel Xeon E5-2695 v4 @ 2.10GHz')
