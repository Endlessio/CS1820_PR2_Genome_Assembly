import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle

fr = open('./res_store_graph_dict_simple.txt','rb')
graph_dict = pickle.load(fr)
x_axis = [i for i in range(6,26)]
y_2 = graph_dict[2]
y_4 = graph_dict[4]
y_6 = graph_dict[6]
# y_8 = graph_dict[8]
# y_10 = graph_dict[10]
# y_12 = graph_dict[12]

plt.title('Result Analysis')
plt.plot(x_axis, y_2, color='blue', label='t=2')
plt.plot(x_axis, y_4, color='green', label='t=4')
plt.plot(x_axis, y_6, color='red', label='t=6')
# plt.plot(x_axis, y_8, color='blue', label='t=8')
# plt.plot(x_axis, y_10, color='yellow', label='t=10')
# plt.plot(x_axis, y_12, color='black', label='t=12')
plt.legend() # 显示图例

plt.xlabel('k value')
plt.ylabel('score')
plt.show()
