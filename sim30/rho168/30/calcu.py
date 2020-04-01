import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

print("=====start=====")

# seperate file number
file_num = 8
f = open('output.csv','w', newline='')
wr= csv.writer(f)
header = ["version","count_row","overlap count","%"]
wr.writerow(header)

for row in range(1,file_num+1):
    strr = 'resource_v'+str(row)+'.csv'

    csv_file = pd.read_csv(strr, header=None)
    
    a = np.array(csv_file)

    count_row = a.shape[0]
    #print(count_row)
    x_value = a

    count = 0


    for k in range(1,count_row):
        for i in x_value[0]:
            for j in x_value[k]:
                if i==j:
                    count=count+1

                
    #print("총 행 갯수 : ", count_row)
    #print("총 겹치는 갯수 : ",count)
    avg = float(count)/float(count_row*40.0)
    num = 'v'+str(row)
    #print("겹치는 평균 비율", avg)
    llist = [num, count_row, count, avg]
    wr.writerow(llist)
    
f.close()
print("=====finish=====")
    


            


