import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import seaborn as sns

df = pd.read_csv('resource_v10.csv', header=None)

a = np.array(df)

count_row = a.shape[0]
count_col = a.shape[1]

x = a

count = 0
# for indexing location of the resource
idx = 0
overArr = np.zeros((count_col,100))
print(overArr.shape)

y = np.zeros((count_row, count_col))
for j in range(0,count_row):
    for i in range(0,count_col):
        y[j][i] = j+1


            
print("for loop")
for k in range(0,count_row):
    overArr[k] = a[k]

new=[[]]
for i in range(0,count_row):
    new[i] = np.concatenate((overArr[i], y[i]), axis=0)
    new[i] = np.reshape(new[i],(2, count_col))
    new.append([])

print("check the loop result")
'''
for i in range(1, count_row+1):
    file_name = './graph_SW/SW_v'+str(i)+".png"
    #plt.subplot(i,1,1)
    hm = sns.heatmap(new[i-1], annot=False, cmap="YlGnBu")
    #plt.pcolor(new[i-1])
    #plt.grid()
    plt.title('count overlapping resource_rho672_200')
    plt.savefig(file_name)
    plt.clf()
'''
print(new[0])
hm = sns.heatmap(a, annot=False, cmap="YlGnBu")
plt.grid()
plt.title('count overlapping resource_rho672_200')
plt.show()


