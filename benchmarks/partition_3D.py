from matplotlib import pyplot as plt
import numpy as np
import statistics

def load_file(fname):
    D = {}
    WC = {}
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            algo, p, T, *values = line.strip().split(";")
            p = int(p)
            T = int(T)
            values = list(map(float, values[:-1]))
            if algo.strip() == 'direct':
                D[p, T] = statistics.median(values[1:])
            else:
                WC[p, T] = statistics.median(values[1:])
    return (D, WC)


def process_file_p(filename, cpu_name):
    D_, WC_ = load_file(filename + '.txt')
    
    T = max(t for (p, t) in D_)
    P = max(p for (p, t) in D_)
    
    D = np.zeros(shape=(P+1, T+1))
    for (p, t), x in D_.items():
        D[p, t] = x
    WC = np.zeros(shape=(P+1, T+1))
    for (p, t), x in WC_.items():
        WC[p, t] = x
    
    
    fig = plt.figure(figsize=(16, 9))
    #T = 48

    plt.title("Partitioning Efficiency ({}, T={})".format(cpu_name, T))
    plt.grid(zorder=-10)
    
    ax = fig.add_subplot(111)
    ax.set_xlabel("Radix bits ($p$)")
    ax.set_ylabel('G item/s')
    
    xs = np.arange(0, P+1)
    ax.bar(xs, D[:, T], width=0.1, label='direct, T={}'.format(T))
    
    xs = np.arange(0, P+1) + 0.1
    ax.bar(xs, WC[:, T], width=0.1, label='WC, T={}'.format(T))
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2*box.height, box.width, box.height * 0.8])
    #ax.set_xticks(L)
    
    ax.set_xticks(np.arange(1, P+1))
    ax.legend(loc='lower left',  bbox_to_anchor=(-0.025, -0.35), ncol=2)
    plt.savefig("partition_p_" + filename + ".pdf")
    #plt.show()


def process_file_T(filename, cpu_name):
    D_, WC_ = load_file(filename + '.txt')
    T = max(t for (p, t) in D_)
    P = max(p for (p, t) in D_)
    D = np.zeros(shape=(P+1, T+1))
    for (p, t), x in D_.items():
        D[p, t] = x
    WC = np.zeros(shape=(P+1, T+1))
    for (p, t), x in WC_.items():
        WC[p, t] = x
    
    fig = plt.figure(figsize=(16, 9))
    
    plt.title("Scalability ({})".format(cpu_name))
    #plt.grid(zorder=-100)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel("#Threads")
    ax.set_ylabel('G item/s')
    xs = np.arange(0, T+1)
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray')

    i = 0
    for p in [8, 10, 12]:
        ax.bar(xs + i * 0.1, D[p, :], width=0.1, label='direct, p={}'.format(p))
        i += 1
        ax.bar(xs + i * 0.1, WC[p, :], width=0.1, label='WC, p={}'.format(p))
        i += 1
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.2*box.height, box.width, box.height * 0.8])
    #ax.set_xticks(L)
    
    ax.set_xticks(np.arange(1, T+1))
    ax.legend(loc='lower left',  bbox_to_anchor=(-0.025, -0.35), ncol=2)
    plt.savefig("partition_T_" + filename + ".pdf")
    #plt.show()

process_file_p('BGQ', 'IBM PowerPC A2 @ 1.6GHz')
process_file_T('BGQ', 'IBM PowerPC A2 @ 1.6GHz')

process_file_p('BGQ_gemv', 'IBM PowerPC A2 @ 1.6GHz, GEMV')
process_file_T('BGQ_gemv', 'IBM PowerPC A2 @ 1.6GHz, GEMV')

process_file_p('BGQ.gemv_sep', 'IBM PowerPC A2 @ 1.6GHz, GEMV [separate]')
process_file_T('BGQ.gemv_sep', 'IBM PowerPC A2 @ 1.6GHz, GEMV [separate]')


# process_file_p('hpac', 'Intel Xeon E5-4620 @ 2.2GHz')
# process_file_T('hpac', 'Intel Xeon E5-4620 @ 2.2GHz')
# process_file_p('hpac_gemv', 'Intel Xeon E5-4620 @ 2.2GHz, GEMV')
# process_file_T('hpac_gemv', 'Intel Xeon E5-4620 @ 2.2GHz, GEMV')


# process_file_p('zen2', 'AMD EPYC 7301 @ 2.2GHz')
# process_file_T('zen2', 'AMD EPYC 7301 @ 2.2GHz')
# process_file_p('zen2_gemv', 'AMD EPYC 7301 @ 2.2GHz, GEMV')
# process_file_T('zen2_gemv', 'AMD EPYC 7301 @ 2.2GHz, GEMV')


# process_file_p('ppti_gpu_3', 'Intel Xeon E5-2695 v4 @ 2.1GHz')
# process_file_T('ppti_gpu_3', 'Intel Xeon E5-2695 v4 @ 2.1GHz')
# process_file_p('ppti_gemv', 'Intel Xeon E5-2695 v4 @ 2.1GHz, GEMV')
# process_file_T('ppti_gemv', 'Intel Xeon E5-2695 v4 @ 2.1GHz, GEMV')