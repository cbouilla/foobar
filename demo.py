from hashlib import sha256 
from binascii import hexlify

inputs = [
  "FOO-0x00000000038AD6921                                                     ".encode() + b"\x33\xf4\x97\x44",
  "BAR-0x0000000001EE3B6BB                                                     ".encode() + b"\xd7\xe5\xbb\x4a",
  "FOOBAR-0x000000000396A78FC                                                  ".encode() + b"\xc4\xae\x89\xd7"
]

hashes = [sha256(sha256(x).digest()).digest() for x in inputs]

xor = bytes([hashes[0][i] ^ hashes[1][i] ^ hashes[2][i] for i in range(32)])


print('  sha256({})'.format(sha256(inputs[0]).hexdigest()))
print('^ sha256({})'.format(sha256(inputs[1]).hexdigest()))
print('^ sha256({})'.format(sha256(inputs[2]).hexdigest()))
print('  ========================================================================')
print('         {}'.format(hexlify(xor).decode()))
print()
N = int.from_bytes(xor, byteorder="big")
i = 0
while N & 1 == 0:
    i += 1
    N >>= 1
print("3SUM on {} bits".format(i))
print("3SUM on {} bits".format(256 - N.bit_length()))
print()
print("hash preimages : ")
print('  sha256d({})'.format(hexlify(inputs[0])))
print('^ sha256d({})'.format(hexlify(inputs[1])))
print('^ sha256d({})'.format(hexlify(inputs[2])))

