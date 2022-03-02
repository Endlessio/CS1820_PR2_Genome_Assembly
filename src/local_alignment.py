from distutils.log import error
import sys
from collections import defaultdict


class LocalAlignment:
    def __init__(self, seq, matrix, gap_penalty):
        self.seq = seq
        self.matrix = matrix
        self.gap_penalty = gap_penalty
        self.val_dp = None
        self.pos_dp = None

        self.imax = 0
        self.idx = (-1, -1)
        self.output1 = ""
        self.output2 = ""

    def local_alignment(self, idx=0):
        if len(self.seq)<idx+1:
            raise Exception("invalid sequence amount for alignment")
        seq1 = self.seq[idx]
        seq2 = self.seq[idx+1]
        row = len(seq1)
        col = len(seq2)

        self.val_dp = [[0]*(col+1) for _ in range(row+1)]
        self.pos_dp = [[0]*(col+1) for _ in range(row+1)]
        for i in range(1,row+1):
            for j in range(1,col+1):
                char1 = seq1[i-1]
                char2 = seq2[j-1] 
                self.val_dp[i][j] = max(
                    self.val_dp[i-1][j-1] + self.matrix[char1][char2],
                    self.val_dp[i-1][j] + self.gap_penalty,
                    self.val_dp[i][j-1] + self.gap_penalty,
                    0
                )
                if self.val_dp[i][j] == self.val_dp[i-1][j-1] + self.matrix[char1][char2]:
                    self.pos_dp[i][j] = 2
                elif self.val_dp[i][j] == self.val_dp[i-1][j] + self.gap_penalty:
                    self.pos_dp[i][j] = 1
                elif self.val_dp[i][j] == self.val_dp[i][j-1] + self.gap_penalty:
                    self.pos_dp[i][j] = 3
                elif self.val_dp[i][j] == 0:
                    self.pos_dp[i][j] = 0

                if self.val_dp[i][j] >= self.imax:
                    if self.val_dp[i][j] == self.imax and i != j:   
                        continue
                    else:                     
                        self.imax = self.val_dp[i][j]
                        self.idx = (i-1,j-1)


    def print_output(self, idx=0):
        i = self.idx[0]
        j = self.idx[1]
        while i+1>0 and j+1>0 and self.val_dp[i+1][j+1]>0:
            if self.pos_dp[i+1][j+1] == 2:
                self.output1 = self.seq[idx][i] + self.output1
                self.output2 = self.seq[idx+1][j] + self.output2
                i -= 1
                j -= 1
            elif self.pos_dp[i+1][j+1] == 1:
                self.output1 = self.seq[idx][i] + self.output1
                self.output2 = "-" + self.output2
                i -= 1
            elif self.pos_dp[i+1][j+1] == 3:
                self.output1 = "-" + self.output1
                self.output2 = self.seq[idx+1][j] + self.output2
                j -= 1
        print(self.output1)
        print(self.output2)
        print(self.imax)

def preprocessing(argv):
    # variable
    matrix = defaultdict(dict)
    gap_penalty = float("-inf")
    seq = []

    # get seq
    with open(argv[1]) as f:
        for line in f.readlines():
            seq.append(line.strip())

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
        
    # get penalty
    if argv[3] != "negInf":
        gap_penalty = int(argv[3])

    return seq, matrix, gap_penalty


def main(argv):
    seq, matrix, gap_penalty = preprocessing(argv)
    
    align = LocalAlignment(seq, matrix, gap_penalty)
    align.local_alignment()
    align.print_output()


if __name__ == "__main__":
    main(sys.argv)
