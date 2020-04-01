import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv

#csv_file = pd.read_csv(r'C:\Users\ubuntu\Desktop\YM_DSRC_Based\LTE-V2Vsim_v3.5_0229_low80_50s_to_Jihyun\resource_v1.csv', header = None)

#csv_file = pd.read_csv(r'C:\Users\ubuntu\Desktop\YM_DSRC_Based\LTE-V2Vsim_v3.5_0229_low80_50s_to_Jihyun\test.csv', header=None )

#csv_file = pd.read_csv(r'resource_v1.csv', header=None )
csv_file = pd.read_csv(r'resource_v2.csv', header=None )
#csv_file = pd.read_csv(r'resource_v3.csv', header=None )
#csv_file = pd.read_csv(r'resource_v4.csv', header=None )
#csv_file = pd.read_csv(r'resource_v5.csv', header=None )


a = np.array(csv_file)

count_row = a.shape[0]

x_value = a
y_value = np.zeros((count_row,40))

for i in range(0,count_row):
    for j in range(0,40):
        y_value[i][j] = int(i+1)

print(x_value.shape)
print(y_value.shape)

print(count_row)
x_val_new = []
print(x_value[0])


print(x_val_new)

plt.scatter(x_value,y_value,s=1)

#ax = plt.subplots()
#ax.set_yticks([0, count_row])

plt.xlim(left=0)
plt.xlim(right=200)
#plt.yticks([0,1,2, count_row+1])
plt.xlabel('select resource')
plt.ylabel('number of snapshot')

plt.title('selected resource locate')

plt.show()


