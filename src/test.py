import pickle
import numpy as np
fr = open('./res_store/graph_dict_simple.txt','rb')
graph_dict = pickle.load(fr)
print(graph_dict[2][11])