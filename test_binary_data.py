import numpy as np
import struct


if __name__ == "__main__":
    counter = 0
    with open("FXLATEST_000_MF002","rb") as f:
        text = f.read()
        data = text[158:]
        for x in range(0,len(data),2):
            sample = data[x:x+2]
            value = struct.unpack('<h', sample)[0]
            if value == 0:
                dbz = np.nan
            elif value == 10692:
                dbz = np.nan
            else:
                dbz = (value/10)/2 - 32.5
                print(counter,dbz)

            counter += 1
    print(counter)