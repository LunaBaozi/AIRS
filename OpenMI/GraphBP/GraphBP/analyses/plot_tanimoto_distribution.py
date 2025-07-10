# import matplotlib.pyplot as plt
# import csv
  
# x = []
# y = []
  
# with open('./results/tanimoto_simil_AURKA.csv','r') as csvfile:
#     plots = csv.reader(csvfile, delimiter = ',')
#     next(plots)
      
#     for row in plots:
#         x.append(row[4])
#         y.append(row[5])
  
# plt.bar(x, y, color = 'g', width = 0.72, label = "Age")
# plt.xlabel('Morgan fps')
# plt.ylabel('Cinfony fps')
# plt.title('Tanimoto similarity')
# plt.legend()
# #plt.show()
# plt.savefig('./results/images/ciao.png')


import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('./results/tanimoto_simil_AURKA.csv', sep=',', header=0, index_col=0)
morgan = data['simil']
cinfony = data['cinfony-simil']

morgan.plot(kind='bar')
plt.ylabel('Frequency')
plt.xlabel('Words')
plt.title('Title')

plt.show()
plt.savefig('./results/images/ciao.png')