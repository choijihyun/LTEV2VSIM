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
overArr = np.zeros((1,100))
print(overArr.shape)

print("for loop")
for k in range(1,count_row):
    for i in x[0]:
        idx = 0
        for j in x[k]:
            if i==j:
                print(i)
                overArr[0][i] = overArr[0][i]+1
                
                count=count+1
            
        idx=idx+1
print("check the loop result")
print(overArr)



hm = sns.heatmap(overArr, annot=False, cmap="YlGnBu")
plt.title('count overlapping resource_rho672_200')
plt.show()
