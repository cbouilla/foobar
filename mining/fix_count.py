import struct

with open('bar.status.alt', 'wb') as f:
    b = struct.pack('<Q', 0x000000074375db1d)
    f.write(b)

