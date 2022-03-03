import pickle


fr = open('./25_2_simple.txt','rb')
res = pickle.load(fr)
print(len(res), len(res[-52]))

# fd = open('./new_graph_dictsimple.txt','rb')
# graph_dict_one = pickle.load(fd)
# print(graph_dict_one)

# graph_dict_tot[2] = graph_dict_one[2]

# f = open("./graph_dict_simple.txt",'wb')
# pickle.dump(graph_dict_tot,f)
# f.close()