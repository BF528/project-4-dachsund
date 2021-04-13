#!/usr/bin/env python
import sys
import gzip
import pandas as pd

#make sure that file exists
try:
    df = gzip.open(sys.argv[1], "r")
except IndexError as ie:
    print("Call sort.py file_to_open file_to_write\n")

d ={} #dictionary for counts (barcode:count)
i=0 #for debugging fd
if df:
    for line in df:
        i+=1#fd
        if i > 100: #fd
            break #fd
        elif line.startswith(b'@'): #after debugging, make this if
            bc = line.decode("utf-8").split(' ')
            if len(bc) > 2:
                code = bc[1][4:] + bc[2][4:]
                if code not in d:
                    d[code] = 1
                else:
                    d[code]+=1 
                    
                    
                    
for key, value in d.items():
    print(key, value)
print(len(d))