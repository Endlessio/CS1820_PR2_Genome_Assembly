import sys
from collections import defaultdict

class Contamination:
    def __init__(self, contam_read, vector, k):
        self.k = k
        self.contam_read = contam_read
        self.vector = vector

    def vector2kmer(self):
        kmer2idx_dict = defaultdict(list)
        length = len(self.vector)
        for i in range(length-self.k+1):
            kmer2idx_dict[self.vector[i:i+self.k]].append(i)
        return kmer2idx_dict
        # print(kmer2idx_dict)

    def left_end_extend(self, seq, start_idx_list):
        imax = self.k

        for _, ele in enumerate(start_idx_list):
            start = 0
            end = start + self.k 
            v_end = ele + self.k 

            while v_end < len(self.vector) and end < len(seq):
                if self.vector[v_end] == seq[end]:
                    end += 1
                    v_end += 1
                else:
                    break
            
            imax = max(end - start, imax)

        return imax



    def right_end_extend(self, seq, start_idx_list):   
        imax = self.k

        for _, ele in enumerate(start_idx_list):
            end = len(seq)-1
            start = end - self.k
            v_start = ele-1

            while v_start >=0 and start >= 0:
                if self.vector[v_start] == seq[start]:
                    v_start -= 1
                    start -= 1
                else:
                    break
            
            imax = max(end - start, imax)

        return imax



    def end_match(self, kmer2idx_dict):
        res = []
        res_idx = []
        for idx, seq in enumerate(self.contam_read):
            # print(seq, self.k)
            # check left end
            if seq[:self.k] in kmer2idx_dict:
                imax = self.left_end_extend(seq, kmer2idx_dict[seq[:self.k]])
                seq = seq[imax:]

            # check right end
            if seq[-self.k:] in kmer2idx_dict:
                imax = self.right_end_extend(seq, kmer2idx_dict[seq[-self.k:]])
                seq = seq[:-imax]
            res.append(seq)
            res_idx.append(idx)
        # print(res)
        # print(res_idx)
        return res, res_idx

    
    def print_output(self, res, res_idx):
        print(",".join(str(i) for i in res_idx))
        print("-"*20)
        for ele in res:
            print(ele)
            



def preprocessing(argv):
    # variable
    contam_read = []
    vector = None
    k = None

    # get contam_read
    with open(argv[1]) as f:
        for line in f.readlines():
            cur = line.strip()
            contam_read.append(cur)


    # get vector
    with open(argv[2]) as f:
        vector = f.readline().strip()

    # get T
    k = int(argv[3])


    return contam_read, vector, k

def main(argv):
    contam_read, vector, k = preprocessing(argv)
    # print(contam_read, vector, k)
    contam = Contamination(contam_read, vector, k)
    kmer2idx_dict = contam.vector2kmer()
    res, res_idx = contam.end_match(kmer2idx_dict)
    contam.print_output(res, res_idx)




if __name__ == "__main__":
    main(sys.argv)
