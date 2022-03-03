from collections import defaultdict
import sys
import random
from os.path import exists
import pickle
import math
from tracemalloc import start

from numpy import gradient

from local_alignment import LocalAlignment
from correction import Correction

def extract(reads):
    true_reads = []
    mutated_reads = []

    long_read = "".join(i for i in reads)
    length = len(long_read)

    # # G: genome length
    # # N: sampled reads
    # # L: read length = 50-mer
    sample_time = 1269
    
    with open ('true_reads.txt', 'w') as true:
        for _ in range(sample_time):
            start_idx = random.randint(0, len(long_read)-50)
            true_read = long_read[start_idx:start_idx+50]
            true_reads.append(true_read)
            true.write(true_read)
            true.write('\n')
    return true_reads

def mutation(reads):
    err_idx = []
    true_reads = list(reads)
    err_reads = list(reads)
    read_num = len(true_reads)
    read_len = 50
    tot_ele = read_num*read_len
    mutate_num = int(tot_ele*0.01)
    print(mutate_num)


    for _ in range(mutate_num):
        read_idx = random.randint(0, read_num-1)
        mutate_idx = random.randint(0, 50-1)
        ATCG_set = {"A", "T", "C", "G"}
        mutate_from = err_reads[read_idx][mutate_idx]
        temp_set = ATCG_set-{mutate_from}
        mutate_to = random.sample(temp_set, 1)
        err_idx.append((read_idx, mutate_idx))
        # print(mutate_to)
        err_reads[read_idx] = err_reads[read_idx][:mutate_idx]+mutate_to[0]+err_reads[read_idx][mutate_idx+1:]



    with open('mutated_reads.txt', 'w') as f:
        for read in err_reads:
            f.write(read)
            f.write('\n')

    return err_reads, err_idx


def preprocessing(argv):
    # variable
    reads = []
    matrix = defaultdict(dict)

    # get db
    with open(argv[1]) as f:
        for line in f.readlines():
            cur = line.strip()
            # assert cur!="", "ERROR: empty string in DB"
            reads.append(cur)

    # get matrix
    with open(argv[2]) as f:
        chars = f.readline().strip().split()[1:]
        temp = []
        for _ in range(len(chars)):
            temp.append(f.readline()[1:].strip().split())
        row = len(temp)
        col = len(temp[0])
        for i in range(row):
            for j in range(col):
                matrix[chars[i]][chars[j]] = int(temp[i][j])


    return reads, matrix



def main(argv):
    S_m = None
    reads, matrix = preprocessing(argv)
    if not exists("./true_reads.txt"):
        print("check true read")
        true_reads = extract(reads)
    else:
        true_reads = []
        with open("./true_reads.txt") as f:
            for line in f.readlines():
                cur = line.strip()
                true_reads.append(cur)


    if not exists("./mutated_reads.txt"):
        mutated_reads, err_idx = mutation(true_reads)
        print("check mutated read")
    else:
        mutated_reads = []
        with open("./mutated_reads.txt") as f:
            for line in f.readlines():
                cur = line.strip()
                mutated_reads.append(cur)
    
    # # do the alignment
    # if not exists("./alignment_score.txt"):
    #     print("check alignment score")
    #     seq = []
    #     for i in range(len(mutated_reads)):
    #         seq.append(true_reads[i])
    #         seq.append(mutated_reads[i])
    #     score_list = []
    #     for idx in range(0, len(seq), 2):
    #         align = LocalAlignment(seq, matrix, float('-inf'))
    #         align.local_alignment(idx)
    #         # align.print_output(idx)
    #         score_list.append(align.imax)
    #     with open('alignment_score.txt', 'w') as f:
    #         S_m = sum(score_list)/len(score_list)
    #         f.write(str(S_m))
    #         f.write('\n')
    #         f.write(" ".join(str(i) for i in score_list))
    # else:
    #     with open('alignment_score.txt', 'r') as f:
    #         S_m = float(f.readline())
    #     # S_m = 49.03804347826087

    # # correct
    # k = [i for i in range(6, 26)]
    # t = [4, 6, 8, 10, 12]
    # # k = [7]
    # # t = [4]
    # d = 2
    # res_dict = defaultdict(list)
    # for cur_k in k:
    #     for cur_t in t:
    #         correction = Correction(mutated_reads, cur_k, cur_t, d)
    #         kmer_list, tot_kmer_dict = correction.form_kmer(cur_k)
    #         # print(tot_kmer_dict)
    #         infrequent_kmer_dict, frequent_kmer_dict = correction.find_infrequent(tot_kmer_dict, cur_t)
    #         # print(frequent_kmer_dict)
    #         closest_pair_dict = correction.find_closest(frequent_kmer_dict, infrequent_kmer_dict, d, cur_k)
    #         # print(closest_pair_dict)
    #         if argv[3] == 'simple':
    #             res, res_idx = correction.simple_replace(closest_pair_dict)
    #         if argv[3] == 'naive':
    #             res, res_idx = correction.naive_replace(closest_pair_dict)
    #             # print(res_idx)
    #         if argv[3] == 'merge':
    #             res, res_idx = correction.merge_replace(closest_pair_dict, infrequent_kmer_dict)
    #         if argv[3] == 'opt_merge':
    #             res, res_idx = correction.opt_merge_replace(closest_pair_dict, infrequent_kmer_dict)
    #         if argv[3] == 'stack':
    #             res, res_idx = correction.stack_replace(closest_pair_dict)
    #         # for i in range(len(res)):
    #         #     # print(res[i]==mutated_reads[i])
    #         res_dict[(cur_k, cur_t)] = res
    #         f=open('./'+str(cur_k)+"_"+str(cur_t)+"_"+str(argv[3])+".txt",'wb')
    #         pickle.dump(res,f)
    #         f.close()
                


    # # calculate final alignment and form graph
    # graph_dict = defaultdict(list)
    # for cur_t in t:
    #     for cur_k in k:
    #         f = open('./'+str(cur_k)+"_"+str(cur_t)+"_"+str(argv[3])+".txt",'rb')
    #         corrected_reads = pickle.load(f)
    #         f.close()      
    #         seq = []
    #         for i in range(len(corrected_reads)):
    #             seq.append(corrected_reads[i])
    #             seq.append(true_reads[i])
    #         score_list = []
    #         for idx in range(0, len(seq), 2):
    #             align = LocalAlignment(seq, matrix, float('-inf'))
    #             align.local_alignment(idx)
    #             score_list.append(align.imax)
    #         S_k = sum(score_list)/len(score_list)
    #         res = -math.log((50-S_k)/(50-S_m))
    #         graph_dict[cur_t].append(res)
    #         print(cur_t, cur_k, score_list)
    # print(graph_dict)
    # f = open("./graph_dict"+str(argv[3])+".txt",'wb')
    # pickle.dump(graph_dict,f)
    # f.close()



    fr = open("./graph_dict"+str(argv[3])+".txt",'rb')
    graph_dict = pickle.load(fr)
    t_optimal = 4
    d = 2
    t = [2]
    S_k_optimal = max(graph_dict[t_optimal])
    k = [i for i in range(6, 26)]
    new_res_dict = defaultdict(list)
    for cur_k in k:
        for cur_t in t:
            correction = Correction(mutated_reads, cur_k, cur_t, d)
            kmer_list, tot_kmer_dict = correction.form_kmer(cur_k)
            # print(tot_kmer_dict)
            infrequent_kmer_dict, frequent_kmer_dict = correction.find_infrequent(tot_kmer_dict, cur_t)
            # print(frequent_kmer_dict)
            closest_pair_dict = correction.find_closest(frequent_kmer_dict, infrequent_kmer_dict, d, cur_k)
            # print(closest_pair_dict)
            if argv[3] == 'simple':
                res, res_idx = correction.simple_replace(closest_pair_dict)
            if argv[3] == 'naive':
                res, res_idx = correction.naive_replace(closest_pair_dict)
                # print(res_idx)
            if argv[3] == 'merge':
                res, res_idx = correction.merge_replace(closest_pair_dict, infrequent_kmer_dict)
            if argv[3] == 'opt_merge':
                res, res_idx = correction.opt_merge_replace(closest_pair_dict, infrequent_kmer_dict)
            if argv[3] == 'stack':
                res, res_idx = correction.stack_replace(closest_pair_dict)
            # for i in range(len(res)):
            #     # print(res[i]==mutated_reads[i])
            new_res_dict[(cur_k, cur_t)] = res
            f=open('./'+str(cur_k)+"_"+str(cur_t)+"_"+str(argv[3])+".txt",'wb')
            pickle.dump(res,f)
            f.close()
    # calculate final alignment and form graph
    new_graph_dict = defaultdict(list)
    for cur_t in t:
        for cur_k in k:
            f = open('./'+str(cur_k)+"_"+str(cur_t)+"_"+str(argv[3])+".txt",'rb')
            corrected_reads = pickle.load(f)
            f.close()      
            seq = []
            for i in range(len(corrected_reads)):
                seq.append(corrected_reads[i])
                seq.append(true_reads[i])
            score_list = []
            for idx in range(0, len(seq), 2):
                align = LocalAlignment(seq, matrix, float('-inf'))
                align.local_alignment(idx)
                score_list.append(align.imax)
            S_k = sum(score_list)/len(score_list)
            res = -math.log((50-S_k)/(50-S_k_optimal))
            new_graph_dict[cur_t].append(res)
            print(cur_t, cur_k, score_list)
    print(new_graph_dict)
    f = open("./new_graph_dict"+str(argv[3])+".txt",'wb')
    pickle.dump(new_graph_dict,f)
    f.close()







    
main(sys.argv)