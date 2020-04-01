import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import seaborn as sns

df = pd.read_csv('resource_v1.csv', header=None)

a = np.array(df)

count_row = a.shape[0]
count_col = a.shape[1]

count = 0
# for indexing location of the resource
idx = 0
overArr = np.zeros((200))
x = np.array(range(1,201))

print(overArr.shape)
print(x.shape)

maxx=0
print("for loop")
for k in range(1,count_row):
    for i in a[0]:
        idx = 0
        for j in a[k]:
            if i==j:
                print(i)
                overArr[i] = overArr[i]+1
                
                count=count+1
                if count> maxx:
                    maxx = int(overArr[i])
            
        idx=idx+1
print("check the loop result")
print(overArr)

print(maxx)
plt.bar(x,overArr)
yy= np.array(range(maxx+2))
plt.yticks(yy)
plt.title('count overlapping resource_rho84_30')
#plt.show()
plt.savefig("overlapping rho84_30")
