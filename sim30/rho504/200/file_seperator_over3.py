## resource file을 분리하는 작업 (\n 한개)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

csv_file = pd.read_csv(r'resource.csv', header=None )

with open('resource.csv') as f:
    count = 0
    idx = 1
    file=[]
    list_cnt = 0
    while 1:
        data = f.readline()
        print(data)
        if not data : break

        if data == "\n":
            if list_cnt > 2 :
                print("=========", idx, "==========")
                print(file)
                
                # write data in a new file
                strr = 'resource_v'+str(idx)+'.csv'
                idx=idx+1
                print("file name is = ", strr)
                
                with open(strr, 'w', newline='') as subfile:
                    writer =  csv.writer(subfile)
                    for i in range(0,len(file)):
                        writer.writerow(file[i])
                        
                #for x in xx:
                    #xx.write(",", join(x)+'\n')
                #xx.close()
                        
            file=[]
            count=0
            list_cnt=0
                
        else :
            row = data.split(",")
            row.remove('\n')
            file.append(row)
            list_cnt=list_cnt+1

