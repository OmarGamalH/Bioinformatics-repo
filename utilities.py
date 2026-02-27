import pandas as pd
import bisect
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