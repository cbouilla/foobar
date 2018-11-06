import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

print("loading")
df = pd.read_csv('dump.txt', delimiter=";", header=None, 
    names=['ID', 'start_time', 'duration', 'hostname_length', 'hostname', 'result_size'], 
    index_col=0, usecols=[0,1,2,4], skipinitialspace=True, 
    dtype={'ID': np.int32, 'start_time': np.uint32, 'duration': np.float64})

print("working")
pq = []
for row in df.itertuples():
    pq.append((row.start_time, 1))
    pq.append((row.start_time + row.duration, -1))

pq.sort()
ts = []
values = []

n = 0
for t, action in pq:
    n += action
    ts.append(t)
    values.append(n)

print("resampling")
ts = np.array(ts, dtype=np.float64)
plop = pd.to_datetime(ts, unit='s')
df = pd.DataFrame(values, columns=['active'], dtype=np.uint32, index=plop)

ts = df.resample('1Min').mean()

print("plotting")
ts.plot()
plt.show()