## resource file을 분리하는 작업 (\n 한개)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

with open('resource_50.csv') as f:
    count = 0
    idx = 1
    file=[]
    list_cnt = 0
    while 1:
        data = f.readline()
        #print(data)
        if not data : break

        if data == "\n":
            # write data in a new file
            strr = 'resource_v'+str(idx)+'.csv'
            idx=idx+1
                
            with open(strr, 'w', newline='') as subfile:
                writer =  csv.writer(subfile)
                for i in range(0,len(file)):
                    writer.writerow(file[i])

                        
            file=[]
            count=0
            list_cnt=0
                
        else :
            row = data.split(",")
            row.remove('\n')
            file.append(row)
            list_cnt=list_cnt+1

