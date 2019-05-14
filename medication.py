f = open("foo.18", 'rb')
#a = f.read(0x62fa9ffc)

# with open('foo.8a', 'wb') as g:
#      g.write(a)
# 
f.seek(4096)
b = f.read()

with open('bar.18a', 'wb') as f:
     f.write(b)