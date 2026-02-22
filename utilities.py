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
                found.append((t_i , t_i + len(p)))

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

    