import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

# 이거는 알아서 바꿔주셈
filelen = 13

for i in range(1,filelen+1):
    x_value=[]
    y_value=[]
    strr = 'resource_v'+str(i)+'.csv'
    name = 'resource_v'+str(i)+'.png'
    csv_file = pd.read_csv(strr, header=None)
    a = np.array(csv_file)

    count_row = a.shape[0]

    x_value = a
    y_value = np.zeros((count_row,40))


    for i in range(0,count_row):
        for j in range(0,40):
            y_value[i][j] = i+1


    plt.figure(figsize=(15,3))
    scatter = plt.scatter(x_value,y_value, s=10)


    plt.xlim(left=0)
    plt.xlim(right=200)
    tick = []
    tick =np.arange(0, count_row+2, step=1)
    plt.yticks(tick)
    plt.xlabel('select resource')
    plt.ylabel('number of snapshot')

    plt.title('selected resource locate')
    #plt.set_size_inches(18, 10)
    #plt.figure(figsize=(9,9))
    
    directory = './graph/'+name
    plt.savefig(directory)
    plt.clf()
   # plt.show()
    


