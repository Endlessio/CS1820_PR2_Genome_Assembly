import numpy as np
import sys
import heapq
from collections import defaultdict


class Correction:
    def __init__(self, err_read, k, t, d):
        self.k = k
        self.t = t
        self.d = d
        self.err_read = err_read
        self.kmer_storage = defaultdict(dict)
        self.infreq_storage = defaultdict(dict)
        self.freq_storage = defaultdict(dict)
        self.distance_storage = defaultdict(dict)

    def form_kmer(self, k):
        kmer_list = []
        tot_kmer_dict = defaultdict(list)
        for idx, read in enumerate(self.err_read):
            split_kmer_dict = defaultdict(list)
            for i in range(len(read)-k+1):
                split_kmer_dict[read[i:i+k]].append(i)
                tot_kmer_dict[read[i:i+k]].append((idx, i))
            kmer_list.append(split_kmer_dict)
        assert len(kmer_list) == len(self.err_read), "ERROR: len(kmer_list) == len(self.err_read)"
        return kmer_list, tot_kmer_dict

    def find_infrequent(self, tot_kmer_dict, t):
        infrequent_dict = defaultdict(set)
        frequent_dict = defaultdict(set)
        for key, val in tot_kmer_dict.items():
            if len(val) < t:
                infrequent_dict[key] = val
            else:
                frequent_dict[key] = val
        return infrequent_dict, frequent_dict


    def find_closest(self, frequent_kmer, infrequent_kmer, d, k):
        def cal_distance(infreq, freq):
            cnt = 0
            for i in range(k):
                if infreq[i] != freq[i]:
                    cnt += 1
            return cnt
        closest_pair_dict = defaultdict(list)
        for infreq in infrequent_kmer.keys():
            # print("here", infreq)
            imin = d
            for freq in frequent_kmer.keys():
                distance = cal_distance(infreq, freq)
                # print(freq, distance, imin)
                if distance <= imin:
                    imin = min(distance, imin)
                    if distance == imin:
                        heapq.heappush(closest_pair_dict[infreq], (-len(frequent_kmer[freq]), freq))
                    else:
                        closest_pair_dict[infreq] = [(-len(frequent_kmer[freq]), freq)]
        return closest_pair_dict

    def naive_replace(self, closest_pair_dict):
        res = []
        res_idx = []
        # print(closest_pair_dict)
        for idx, read in enumerate(self.err_read):
            dummy_read = read
            flag = True
            for i in range(len(read)-self.k+1):
                if dummy_read[i:self.k+i] in closest_pair_dict:
                    flag = False
                    # print("here", heapq.heappop(closest_pair_dict[dummy_read[i:self.k+i]]))
                    target = closest_pair_dict[dummy_read[i:self.k+i]][0][1]
                    read = read[:i]+target+read[self.k+i:]
            if not flag:
                res_idx.append(idx)
            res.append(read)
        return res, res_idx

    def merge_find_closest(self, frequent_kmer, infrequent_kmer, d, k):
        def cal_distance(infreq, freq):
            cnt = 0
            for i in range(k):
                if infreq[i] != freq[i]:
                    cnt += 1
            return cnt
        closest_pair_dict = defaultdict(list)
        for infreq in infrequent_kmer.keys():
            for freq in frequent_kmer.keys():
                distance = cal_distance(infreq, freq)
                heapq.heappush(closest_pair_dict[infreq], (-len(frequent_kmer[freq]), distance, freq))
        return closest_pair_dict

    def merge_replace(self, closest_pair_dict, infrequent_kmer_dict):
        def find_target_pos(dp_cnt, dp_check):
            length = len(dp_cnt)
            target_pos = []
            imax = np.max(dp_cnt)
            # no infrequent
            if imax == 0:
                return -1
            # no overlap
            if imax == 1:
                return 0

            idx = 0
            while idx < length:
                ele = dp_cnt[idx]
                if ele == imax:
                    start = idx
                    while dp_cnt[idx] == imax:
                        idx += 1
                    end = idx - 1
                    target_pos.append((start, end, imax))
                else:
                    idx += 1   

            return target_pos


        def replace(target_pos, read, dp_check):
            # print("target pos", target_pos, read)
            res_read = read
            for ele in target_pos:
                k_len = ele[1]-ele[0]+1
                # print(ele, k_len)
                # form new kmer if not exist, identify frequency, distance map
                if k_len not in self.kmer_storage:
                    # print("not exist")
                    kmer_list, tot_kmer_dict = self.form_kmer(k_len)
                    self.kmer_storage[k_len] = tot_kmer_dict

                    infrequent_kmer_dict, frequent_kmer_dict = self.find_infrequent(tot_kmer_dict, self.t)
                    self.infreq_storage[k_len] = infrequent_kmer_dict
                    self.freq_storage[k_len] = frequent_kmer_dict

                    distance_map = self.merge_find_closest(frequent_kmer_dict, infrequent_kmer_dict, self.d-1, k_len)
                    self.distance_storage[k_len] = distance_map
                else:
                    # print("exist")
                    tot_kmer_dict = self.kmer_storage[k_len]
                    infrequent_kmer_dict = self.infreq_storage[k_len]
                    frequent_kmer_dict = self.freq_storage[k_len]
                    distance_map = self.distance_storage[k_len]
                
                start = ele[0]
                end = ele[1]
                imax = float("-inf")
                res = ""
                for freq, dis, ele in distance_map[res_read[start:end+1]]:
                    temp = res_read[:start] + ele + res_read[end+1:]
                    cnt = 0
                    for i in range(max(0,start-self.k+1), min(len(read), end+self.k)):
                        if temp[i:i+self.k] in self.freq_storage[self.k]:
                            cnt += 1
                    if cnt > imax:
                        # print("check", ele, read, temp)
                        res_read = temp
                        imax = cnt
            # print("final", read, res_read)
            return res_read


        # print("pair", closest_pair_dict)
        res = []
        res_idx = []
        for idx, read in enumerate(self.err_read):
            dp_cnt = np.zeros((len(read),), int)
            dp_check = np.zeros((len(read),), int)
            dummy_read = read
            flag = True
            for i in range(len(read)-self.k+1):
                if dummy_read[i:self.k+i] in infrequent_kmer_dict:
                    dp_cnt[i:i+self.k] += 1
                    dp_check[i] = 1
            # print(read, dp_cnt, dp_check)
            target_pos = find_target_pos(dp_cnt, dp_check)
            if target_pos == -1:
                res.append(read)

            # TODO fix this
            elif target_pos == 0:
                return self.naive_replace(closest_pair_dict)
            else:
                res_idx.append(idx)
                res.append(replace(target_pos, read, dp_check))
         

        return res, res_idx
    
    def print_output(self, res, res_idx):
        print(",".join(str(i) for i in res_idx))
        print("-"*20)
        for ele in res:
            print(ele)
            



def preprocessing(argv):
    # variable
    err_read = []

    # get error_read
    with open(argv[1]) as f:
        for line in f.readlines():
            cur = line.strip()
            err_read.append(cur)

    # get k
    k = int(argv[2])

    # get t
    t = int(argv[3])

    # get d
    d = int(argv[4])


    return err_read, k, t, d

def main(argv):
    err_read, k, t, d = preprocessing(argv)
    # print(err_read, k, t, d)
    correction = Correction(err_read, k, t, d)

    kmer_list, tot_kmer_dict = correction.form_kmer(k)
    # print(tot_kmer_dict)
    infrequent_kmer_dict, frequent_kmer_dict = correction.find_infrequent(tot_kmer_dict, t)
    # print(infrequent_kmer_dict)
    closest_pair_dict = correction.find_closest(frequent_kmer_dict, infrequent_kmer_dict, d, k)
    # print(closest_pair_dict)
    # res, res_idx = correction.naive_replace(closest_pair_dict)
    # correction.print_output(res, res_idx)
    # if argv[5] == "naive":
    #     res, res_idx = correction.naive_replace(closest_pair_dict)
    #     correction.print_output(res, res_idx)
    # if argv[5] == "merge":
    res, res_idx = correction.merge_replace(closest_pair_dict, infrequent_kmer_dict)
        # correction.print_output(res, res_idx)

    correction.print_output(res, res_idx)




if __name__ == "__main__":
    main(sys.argv)
