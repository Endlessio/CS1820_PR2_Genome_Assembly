from nbformat import read
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
        self.tot_kmer_dict = None
        self.closest_pair_dict = None
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
        self.tot_kmer_dict = tot_kmer_dict
        return kmer_list, tot_kmer_dict

    def find_infrequent(self, tot_kmer_dict, t):
        infrequent_dict = defaultdict(set)
        frequent_dict = defaultdict(set)
        for key, val in tot_kmer_dict.items():
            if len(val) < t:
                infrequent_dict[key] = val
            else:
                frequent_dict[key] = val
        self.freq_storage[self.k] = frequent_dict
        self.infreq_storage[self.k] = infrequent_dict
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
        self.closest_pair_dict = closest_pair_dict
        return closest_pair_dict

    def simple_replace(self, closest_pair_dict):
        res_idx = set()
        iteration = 3
        read_num = len(self.err_read)
        res = list(self.err_read)
        for _ in range(iteration):
            for read_idx in range(read_num):
                cur_read = res[read_idx]
                for check_idx in range(len(cur_read)-self.k+1):
                    if cur_read[check_idx: check_idx+self.k] in closest_pair_dict:
                        res_idx.add(read_idx)
                        replace_to = closest_pair_dict[cur_read[check_idx: check_idx+self.k]][0][1]
                        cur_read = cur_read[:check_idx]+replace_to+cur_read[check_idx+self.k:]
                res[read_idx] = cur_read
        return res, res_idx

    
    def stack_replace(self, closest_pair_dict):
        def stack_check(read):
            dp = [0]*len(read)
            stack = [float("-inf")]
            for i in range(len(read)-self.k+1):
                if read[i:self.k+i] in self.closest_pair_dict:
                    if i-stack[-1] < self.k:
                        dp[i] = 0
                    else:
                        dp[i] = 1
                        stack.append(i)
            return dp

        res_idx = set()
        iteration = 3
        read_num = len(self.err_read)
        res = list(self.err_read)
        for _ in range(iteration):
            for read_idx in range(read_num):
                cur_read = res[read_idx]
                dp = stack_check(cur_read)
                for check_idx in range(len(cur_read)-self.k+1):
                    if dp[check_idx] == 1:
                        res_idx.add(read_idx)
                        # print("dp", dp)
                        # print(check_idx, cur_read[check_idx: check_idx+self.k], closest_pair_dict[cur_read[check_idx: check_idx+self.k]])
                        replace_to = closest_pair_dict[cur_read[check_idx: check_idx+self.k]][0][1]
                        cur_read = cur_read[:check_idx]+replace_to+cur_read[check_idx+self.k:]
                        # print(len(cur_read))
                res[read_idx] = cur_read
        return res, res_idx



    def naive_replace(self, closest_pair_dict):
        def dp_check(read):
            # print(closest_pair_dict)
            dp = [0]*len(read)
            for i in range(0, len(read)-self.k+1):
                if read[i:self.k+i] in self.closest_pair_dict:
                    if i == 0:
                        dp[i] = 1
                    else:
                        dp[i] = dp[i-1]+1
                        if dp[i] == self.k:
                            start = i-self.k+1
                            end = i
                            mid = (start+end)//2
                            dp[start:end+1] = [0]*self.k
                            dp[mid] = 1
                else:
                    if i == 0:
                        continue
                    else:
                        if dp[i-1] > 1:
                            start = i-dp[i-1]
                            end = i-1
                            mid = (start+end)//2
                            dp[start:end+1] = [0]*dp[i-1]
                            dp[mid] = 1
            # print(dp)
            return dp

        res = list(self.err_read)
        res_idx = set()
        iteration_num = 3
        # print(closest_pair_dict)
        for iter in range(iteration_num):
            # print("first", self.err_read[6])
            for idx in range(len(self.err_read)):
                # print("read", self.err_read[6])
                read = self.err_read[idx]
                dp_check_res = dp_check(read)
                # print(iter, idx, dp_check_res)
                for i in range(len(read)-self.k+1):
                    # print(i, len(dp_check_res))
                    if dp_check_res[i] == 1:
                        # print("bug", self.infreq_storage[self.k][read[i:self.k+i]], self.closest_pair_dict[read[i:self.k+i]])
                        if len(self.closest_pair_dict[read[i:self.k+i]])<=0:
                            continue
                        replace_target = self.closest_pair_dict[read[i:self.k+i]][0][1]
                        new_read = read[:i]+replace_target+read[self.k+i:]
                        # print("fffff", iter, idx, i, len(new_read), len(read))
                        read = new_read
                        res_idx.add(idx)
                        # print("fuck", res_idx)
                res[idx] = read
            self.err_read = list(res)

            # self.form_kmer(self.k)
            # infrequent_dict, frequent_dict = self.find_infrequent(self.tot_kmer_dict, self.t)
            # self.find_closest(frequent_dict, infrequent_dict, self.d, self.k)
        return res, sorted(list(res_idx))

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


    def opt_merge_replace(self, closest_pair_dict, infrequent_kmer_dict):
        def find_target_pos(dp_cnt, dp_check):
            length = len(dp_cnt)
            target_pos = []
            imax = np.max(dp_cnt)
            # no infrequent
            if imax < 0:
                return -1
            #no overlap
            if imax == -1:
                return 0


            idx = 0
            while idx < length:
                ele = dp_cnt[idx]
                if imax == ele:
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
                for _, _, ele in distance_map[res_read[start:end+1]]:
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
                else:
                    dp_cnt[i:i+self.k] -= 1
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
    res, res_idx = correction.stack_replace(closest_pair_dict)
    # correction.print_output(res, res_idx)
    # if argv[5] == "naive":
    #     res, res_idx = correction.naive_replace(closest_pair_dict)
    #     correction.print_output(res, res_idx)
    # if argv[5] == "merge":
    # res, res_idx = correction.merge_replace(closest_pair_dict, infrequent_kmer_dict)
        # correction.print_output(res, res_idx)
    # res, res_idx = correction.opt_merge_replace(closest_pair_dict, infrequent_kmer_dict)
    correction.print_output(res, res_idx)




if __name__ == "__main__":
    main(sys.argv)
