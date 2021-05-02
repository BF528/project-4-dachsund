#!/usr/bin/env python
import csv
import sys

dictArr = [] #Stores an array of dictionary

l = len(sys.argv)# measures the number of arguments given

#Checks if there are more than 2 args given, else print error message
if l > 2:
    #For each argument (sys.argv[0] is this file's name)
    for i in range(1,l):
        #Try to open. If unable, throw error message
        try:
            #open each file
            with open(sys.argv[i], "r") as file:
                d ={} #temp d to hold
                r = csv.reader(file)
                for line in r:
                    #try so if there is a mistake, nothing happens
                    try:
                        d[line[0]] = int(line[1])
                    except:
                        pass
                dictArr.append(d) #append to array
        except:
            print("Cannot open " + sys.argv[i])
else:
    print("Not enough arguments! Call with 2 or more csv files")
    print("Call: 'python3 merge.py file1.csv file2.csv [...] fileN.csv'")
    print("csv files must have 'name', integer \\n for each row")
    
l = len(dictArr)#length of dictArr; l is a re-used variable

#Checks that there are at least 2 sucessful dictionaries
if l > 1:
    #Combine all matching keys to dictArr[0]
    for key, value in dictArr[0].items():
        for i in range(1, l):
            if key in dictArr[i]:
                dictArr[0][key] += dictArr[i][key]
                dictArr[i].pop(key)#pop off matching keys in such dict
    #Combine remaining keys to dictArr[0]
    for i in range(1, l):
            for key, value in dictArr[i].items():
                dictArr[0][key] = value
    #Prints result to file  
    f = "merge_output_" + str(l) + "_files.csv"
    with open(f, "w") as file:
        w = csv.writer(file, delimiter=',')
        for key, value in dictArr[0].items():
            w.writerow([key, value])
elif l == 1:
    print("Only one dictionary sucessful. Check your .csv files again")
else:
    print("No dictionaries sucessful in converting. Check your files again")
