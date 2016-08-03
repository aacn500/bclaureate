def int_to_bytes(n):
    s = bin(n)[2:]
    l = [s[:8], s[8:16], s[16:24], s[24:]]
    return [chr(int(c, 2)) for c in l]

def int_to_bytes(n):
    s = bin(n)[2:]
    l = [s[:8], s[8:16], s[16:24], s[24:]]
    return ''.join([chr(int(c, 2)) for c in l])
