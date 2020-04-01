import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import seaborn as sns

df = pd.read_csv('resource_v4.csv', header=None)

a = np.array(df)

count_row = a.shape[0]
count_col = a.shape[1]
print("row", count_row)
print("col", count_col)

count = 0
# for indexing location of the resource
idx = 0
A = np.zeros((100))
B = np.zeros((100))
x = np.array(range(1,101))


for j in range(0, count_col):
    A[a[0][j]-1] = A[a[0][j]-1]+1
for j in range(0, count_col):
    B[a[1][j]-1] = B[a[1][j]-1]+1
plt.bar(x, A, color='b')
plt.bar(x, B, color='r',  bottom=A)


print("check the loop result")
#print(overArr)

#plt.bar(x,overArr)
yy= np.array(range(2+2))
plt.yticks(yy)
plt.title('count overlapping resource_rho168_30')
#plt.show()
plt.savefig("overlapping rho168_30")
