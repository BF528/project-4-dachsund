#!/usr/bin/env python
import sys
import gzip
import csv

#make sure that file exists
try:
    df = gzip.open(sys.argv[1], "r")
except IndexError as ie:
    print("Call sort.py file_to_open name_to_write(optional)\n")

d ={} #dictionary for counts (barcode:count)
i=0 #for debugging fd
if df:
    for line in df:
        i+=1#fd
        if i > 100: #fd
            break #fd
        #if line has barcode, continue. Else skip
        elif line.startswith(b'@'): #after debugging, make this if
            bc = line.decode("utf-8").split(' ')
            #Only use lines w/barcodes
            if len(bc) > 2:
                #concat bc1 and bc2
                code = bc[1][4:] + bc[2][4:]
                #
                if code not in d:
                    d[code] = 1
                else:
                    d[code]+=1 

#will write to csv 
y = False
#Checks if a file name was specified
if sys.argv[2]:
    f = sys.argv[2]
else: #if not, use original file name
    y = True
    f = sys.argv[1]
    #if file path, delete anything but the actual name
    #Only if you use sys.argv[1]
    if "/" in f:
        f = f.split("/")
        f = f[len(f)-1]
        f = ''.join(f)
#Formats the output into a csv no matter what
f = f.split(".") 
#if there's a .{anything} will remove
if len(f) > 1:
    f = f[0]
f = ''.join(f)
if y == True: #if we used the argv[1] argument, add "_output"
    f += "_output"
f += ".csv" #ensures that a .csv

with open(f, "w") as file:
    w = csv.writer(file, delimiter=',')
    for key, value in d.items():
        w.writerow([key, value])