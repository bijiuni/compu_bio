from dna2vec.multi_k_model import MultiKModel
from scipy import spatial

filepath = 'pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v'
mk_model = MultiKModel(filepath)

def generate_candidates(kmer):
    sequence_list = ['A', 'T', 'G', 'C']
    candidates = ["A", "T", "G", "C"]
    candidates_new = []

    for _ in range(kmer - 1):
        for candidate in candidates:
            for seq in sequence_list:
                candidates_new.append(candidate+str(seq))
        candidates = candidates_new
        candidates_new = []

    #print(candidates)
    #print(len(candidates))
    return candidates


def rank_neighbors(candidates, added_vector, k, kmer_1, kmer_2):
    possibility = {}
    for candidate in candidates:
        candidate_vector = mk_model.vector(candidate)
        similarity = 1 - spatial.distance.cosine(added_vector, candidate_vector)
        possibility[candidate] = similarity

    k_neighbors = sorted(possibility.items(), key = lambda kv:(kv[1], kv[0]), reverse=True)
    
    success = False
    for index in range(k):
        if k_neighbors[index][0] == kmer_1 + kmer_2:
            #print("the valid kmer: {}.".format(kmer_1 + kmer_2))
            success = True
        elif k_neighbors[index][0] == kmer_2 + kmer_1:
            #print("the valid kmer: {}.".format(kmer_2 + kmer_1))
            success = True
    return success



def main():
    kmer_len1 = 3
    kmer_len2 = 3
    kmer = kmer_len1 + kmer_len2
    k = 5

    kmer_list1 = generate_candidates(kmer_len1)
    kmer_list2 = generate_candidates(kmer_len2)
    candidates = generate_candidates(kmer)
    total = 0
    matched = 0
    for kmer_1 in kmer_list1:
        for kmer_2 in kmer_list2:
            total += 1
            added_vector = mk_model.vector(kmer_1) + mk_model.vector(kmer_2)

            success = rank_neighbors(candidates, added_vector, k, kmer_1, kmer_2)
            if success:
                matched += 1
            else:
                print(kmer_1, " ", kmer_2)
    print("The accuracy is {}".format(matched/total))

  
if __name__== "__main__":
    main()