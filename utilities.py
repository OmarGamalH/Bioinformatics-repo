import pandas as pd
import bisect
import numpy as np
class FASTA:

    def __init__(self , filename):
        self.filename = filename

    def read_all(self):
        file = {}
        current_id = None
        with open(self.filename , "r") as f:
            for line in f:
                if line.startswith(">"):
                    id , description = line.split()[0].lstrip(">") , " ".join(line.split()[1:])
                    current_id = id
                    file[id] = {"description" : description , "sequence" : ""}
                else:
                    file[current_id]["sequence"] = file[current_id]["sequence"] + line.strip("\n").upper()
        return file

    def read_dna(self):
        DNA = ""
        with open(self.filename) as f:
            for line in f:
                if not line.startswith(">"):
                        DNA = DNA + line.strip("\n").upper()
        
        return DNA

    def read_df(self):
        all_data = self.read_all()
        ids = []
        descriptions = []
        sequence = []
        for id in all_data:
            ids.append(id)
            descriptions.append(all_data[id]['description'])
            sequence.append(all_data[id]['sequence'])
        
        return pd.DataFrame({"id" : ids , "description" : descriptions , "sequence" : sequence})




class index:
    
    def __init__(self , t , k):
        self.k = k
        self.t = t
        self.k_mers = []
        
        for i in range(len(t) - k + 1):
           self.k_mers.append((t[i:i+k] , i))
    
        self.k_mers.sort()


    def find_p(self , p):
        """
        a function that takes a pattern(p) and finds left most index of it in the k-mer index table
        
        p: pattern that you want to search for

        output: a hit which is the leftmost index in the k-mer index table
        
        """

        hit = bisect.bisect_left(self.k_mers , (p[:self.k] , -1))
        return hit

    def find_all_p(self , hit):
        
        """
        takes a hit from func find_p and search for all other hits 
        
        """



        processed_hits = [hit]
        i = hit + 1
        while i < len(self.k_mers) and self.k_mers[hit][0] == self.k_mers[i][0]:
                processed_hits.append(i)
                i = i + 1
        return processed_hits
    

    def show_k_mers(self):
        return self.k_mers
    

    def query(self , p):
        """
        takes pattern(p) and search for it exactly inside a given text(t)
        
        output: a list of tuples of (start , end) , where start = start index in the text , end = end index in the text
        """
        found = []
        hit = self.find_p(p)
        hits = self.find_all_p(hit)
        k_mers = self.show_k_mers()

        for k in hits:
            k_mer = k_mers[k]
            t_i = k_mer[1]
            if p[self.k:] == self.t[t_i + self.k : t_i + len(p) ]:
                found.append(t_i)

        return found  
    
    def df_query(self , patterns:list):
        """
        similar as query but takes more than 1 pattten + returns the output as pandas dataframe
        """
        founds = []
        for p in patterns:
            if len(p) >= self.k:
                founds.append(self.query(p))
            else:
                print(f"Sequence : {p} has length less than k")
        starts = []
        ends = []
        sequences = []
        for found in founds:
             for start , end in found:
                 starts.append(start)
                 ends.append(end)
                 sequences.append(self.t[start : end])


        return pd.DataFrame({"sequence" : sequences , "start" : starts , "end" : ends})

    

def pigeonhole(p , n , t):
    """
    p: pattern
    n: number of mismatches

    """
    ind = index(t , 8)
    segment_length = int(round(len(p)/(n + 1)))
    all_matches = []
    all_hits = []
    for i in range(n + 1):
        start = i * segment_length
        end = min(start + segment_length , len(p))
        segment = p[start:end]
        matches = ind.query(segment)
        all_hits.append(matches)
        print(f"segment:{segment} , matches:{matches}")
        print(f"0:start -> {p[0:start]} , start:end->{p[start:end]} , end:len(p)->{p[end:len(p)]}")
        for match in matches:
            n_mismatches = 0
            if match < start or match - start + len(p) > len(t):
                continue
        
            for j in range(0 , start):
                if p[j] != t[match - start + j]:
                    n_mismatches += 1
                    
            for j in range(end , len(p)):
                if p[j] != t[match - start + j]:
                    n_mismatches += 1
            

            if n_mismatches <= n:
                all_matches.append(match - start)

        count_hits = 0 
        for hit in all_hits:
            count_hits += len(hit)
        print(f"hits count : {count_hits}")
    return list(set(all_matches))



def naive_with_counts(p , t):
    number_of_alignments = len(t) - len(p) + 1
    occurences = []
    num_character_comparisons = 0
    for i in range(number_of_alignments):
        matched = True
        
        for j in range(len(p)):
            num_character_comparisons += 1
            if t[i + j] != p[j]:
                matched = False
                break
        
        if matched:
            occurences.append(i)
    return occurences , number_of_alignments , num_character_comparisons
        

def naive_with_nmismatches(p , t , n):
    number_of_alignments = len(t) - len(p) + 1
    occurences = []
    num_character_comparisons = 0
    for i in range(number_of_alignments):
        mismatches = 0
        for j in range(len(p)):
            num_character_comparisons += 1
            if t[i + j] != p[j]:
                mismatches += 1
        
        if mismatches <= n:
            occurences.append(i)
    return occurences


class edit_distance:

    def __init__(self , t_1 , t_2):
        self.t_1 = t_1
        self.t_2 = t_2
        dp = np.zeros((len(t_1) + 1 , len(t_2) + 1))
        dp[: , 0] = [i for i in range(len(t_1) + 1)]
        dp[0 , :] = [i for i in range(len(t_2) + 1)]
        self.rows , self.columns = dp.shape
        self.dp = dp

        for i in range(1 , self.columns):
            char_2 = self.t_2[i - 1]
            for j in range(1 , self.rows):
                char_1 = t_1[j - 1]
                if char_1 == char_2:
                    self.dp[j][i] = self.dp[j - 1][i - 1]
                else:
                    self.dp[j][i] = min(1 + self.dp[j-1][i] , 1 + self.dp[j][i - 1] , 1 + self.dp[j - 1][i - 1])

    def get_dp_table(self):
        return self.dp
    
    def get_final_answer(self):
        return int(self.dp[self.rows - 1][self.columns - 1])

    def __edit(self , t:str, index:int,  edit_type:str, edit = None):
        
        new_text = t
        new_text_array = list(new_text)
        
        if edit_type == "substituation" and edit != None:
            new_text_array[index] = edit
            return "".join(new_text_array)
        elif edit_type == "deletion":
            new_text_array[index] = 'removed'
            new_text_array.remove("removed")
            return "".join(new_text_array)
        elif edit_type == "insertion" :
            new_text = t[:index] + edit + t[index:]
            return new_text
            
        
    def back_tracking(self):
        i = self.columns - 1 
        j = self.rows - 1

        
        answer = ""
        while j >= 0 and i >= 0:
            char_2 = self.t_2[i - 1]
            if self.dp[j][i] == 1 + self.dp[j][i-1]:
                answer = answer + f"1 substitution {j - 1}\n"
                i = i - 1
            elif self.dp[j][i] == 1 + self.dp[j - 1][i]:
                answer = answer + f"1 deletion at {j - 1} \n"
                j = j - 1
            elif self.dp[j][i] == 1 + self.dp[j-1][i-1]:
                answer = answer + f"1 insetion at {j - 1}\n"
                i = i - 1
                j = j - 1
            else:
                j = j - 1
                i = i - 1
        return answer
    


class FASTAQ:
    
    def __init__(self , filename):
        self.filename = filename

    def get_all_data(self):

        with open(self.filename , 'r') as f:
            all_headers   = []
            all_sequences = []
            all_qualities = []
            while True:
                header = f.readline().strip('\n')
                all_headers.append(header)
                sequence = f.readline().strip('\n')
                all_sequences.append(sequence)
                f.readline()
                quality = f.readline().strip('\n')
                all_qualities.append(quality)
                if len(header) == 0:
                    break
        return (all_headers , all_sequences , all_qualities)



class overlap:

    def __init__(self , sequences , k):
        self.sequences = sequences
        self.min_length = k
        self.pre_table = {}
        for seq in sequences:
            for i in range(len(seq) - k + 1):
                start = i
                end = start + self.min_length
                subseq = seq[start : end]
                if subseq not in self.pre_table:
                    self.pre_table[subseq] = set()
                self.pre_table[subseq].add(seq)
        self.table = {}
        for subseq in self.pre_table:
            if len(self.pre_table[subseq]) > 1:
                self.table[subseq] = self.pre_table[subseq]
        


    def get_table(self):
        return self.table

    def is_empty(self , k_mer):
        if k_mer in self.table and len(self.table[k_mer]) > 0:
            return False
        return True
    
    def get_all_sequences(self , k_mer):
        return self.table[k_mer]
    
    def longest_overlap(self , a : str , b : str):
        prefix = b[:self.min_length]
        result = 0
        while True:
            result = a.find(prefix , result)
            if result == -1:
                return 0
            
            if b.startswith(a[result:]):
                return len(a) - result
            
            result += 1


        
    
    def create_graph(self):
        graph = {}
        for read in self.sequences:
            subseq = read[-self.min_length : -1] + read[-1]
            all_prefix_candidates = self.table.get(subseq)
            if all_prefix_candidates is not None and len(all_prefix_candidates) >= 1:
                for prefix_candidate in all_prefix_candidates:
                    overlap_length = self.longest_overlap(read , prefix_candidate)
                    if prefix_candidate != read and overlap_length > 0:
                        graph[(read , prefix_candidate)] = overlap_length



        return graph
    

class FASTQ:
    
    def __init__(self , filename):
        self.filename = filename

    def get_all_data(self):

        with open(self.filename , 'r') as f:
            all_headers   = []
            all_sequences = []
            all_qualities = []
            while True:
                header = f.readline().strip('\n')
                if len(header) == 0:
                    break
                all_headers.append(header)
                sequence = f.readline().strip('\n')
                all_sequences.append(sequence)
                f.readline()
                quality = f.readline().strip('\n')
                all_qualities.append(quality)
        return (all_headers , all_sequences , all_qualities)
    
    


class scs_index:
    
    def __init__(self , sequences , k):
        self.sequences = sequences
        self.min_length = k
        self.table = {}
        for seq in sequences:
            for i in range(len(seq) - k + 1):
                start = i
                end = start + self.min_length
                subseq = seq[start : end]
                if subseq not in self.table:
                    self.table[subseq] = set()
                self.table[subseq].add(seq)

    def get_all_sequeces(self , k_mer):
        return self.table[k_mer]
    
    def longest_overlap(self , a : str , b : str):
        prefix = b[:self.min_length]
        result = 0
        while True:
            result = a.find(prefix , result)
            if result == -1:
                return 0
            
            if b.startswith(a[result:]):
                return len(a) - result
            
            result += 1
        
    def search_index(self , read):
        suffix = read[-self.min_length:-1] + read[-1]
        sequences = self.get_all_sequeces(suffix)
        a = read
        b = None
        max_length = 0

        for sequence in sequences:
            curr_length = self.longest_overlap(a , sequence)
            if curr_length > max_length and a != sequence:
                max_length = curr_length
                b = sequence
        
        return a , b , max_length
    
    def compute_max_lengths(self):
        result_a = None
        result_b = None
        result_max_length = 0
        for read in self.sequences:
            a , b , max_length = self.search_index(read)
            if max_length > result_max_length:
                result_a , result_b , result_max_length = a , b , max_length
        
        return result_a , result_b , result_max_length
    
def merge(a , b , max_length):
    return a + b[max_length:]



def overlap_f(a , b , min_length = 1):

    start = 0 

    while True:

        start = a.find(b[:min_length] , start)

        if start == -1:
            return 0 
        
        if b.startswith(a[start:]):
            return start
        
        start = start + 1

def greedy_SCS_optimized(reads):
    obj = scs_index(reads , 30)

    a , b , max_length = obj.compute_max_lengths()
    while max_length > 0 or len(reads) > 2:
        print(len(reads))
        if b == None:
            a = reads[0]
            b = reads[1]
            max_length = overlap_f(a , b , 30)
        reads.remove(a)
        reads.remove(b)
        result = merge(a , b, max_length)
        reads.append(result)
        obj = scs_index(reads , 30)
        a , b , max_length = obj.compute_max_lengths()

    return reads