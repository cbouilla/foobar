import sys
import array
import os.path

filenames = sys.argv[1:]
prefix = array.array('Q')
assert prefix.itemsize == 8

prefix.append(len(filenames))

s = 0
i = 0
for f in filenames:
    size =  os.path.getsize(f)
    prefix.append(s)
    s += size
    i += 1
prefix.append(s)

sys.stdout.buffer.write(prefix.tobytes())

for f in filenames:
    with open(f, 'rb') as g:
        sys.stdout.buffer.write(g.read())
