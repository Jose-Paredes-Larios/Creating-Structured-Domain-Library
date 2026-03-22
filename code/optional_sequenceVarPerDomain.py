import pandas as pd
from itertools import combinations
from Bio.Align import PairwiseAligner

def percentIdentity(seq1, seq2, aligner):
    aligned=aligner.align(seq1, seq2)[0]
    aSeq1=aligned[0]
    aSeq2=aligned[1]

    matches=0
    comparedNoGaps=0
    comparedWithGaps=len(aSeq1)

    for res1, res2 in zip(aSeq1, aSeq2):
        if res1 != '-' and res2 != '-':
            comparedNoGaps+=1
            if res1==res2:
                matches += 1

    pidNoGaps=matches/comparedNoGaps if comparedNoGaps > 0 else 0
    pidWithGaps=matches/comparedWithGaps if comparedWithGaps > 0 else 0

    return pidNoGaps, pidWithGaps

input="/Users/joseparedes/Desktop/kappelLab/structuredDomainLibrary/filteredDomainSequences2.tsv"
df=pd.read_csv(input, sep="\t")
aligner=PairwiseAligner()
aligner.mode="global"

domains=df.groupby("interProList")
results=[]
for interproID, group in domains:
    sequences=group["Domain Sequence"].tolist()
    interproName=group["InterProName"].iloc[0]
    print(interproName)
    numSeqs=len(sequences)

    pidNoGapsList=[]
    pidWithGapsList=[]
    for seq1, seq2 in combinations(sequences, 2):
        pidNoGaps, pidWithGaps=percentIdentity(seq1, seq2, aligner)
        pidNoGapsList.append(pidNoGaps)
        pidWithGapsList.append(pidWithGaps)

    results.append({
        "InterProID": interproID,
        "InterProName": interproName,
        "NumSequences": numSeqs,
        "MeanPIDNoGaps": sum(pidNoGapsList)/len(pidNoGapsList),
        "MeanPIDWithGaps": sum(pidWithGapsList)/len(pidWithGapsList)})

output="/Users/joseparedes/Desktop/kappelLab/structuredDomainLibrary/domainLibraryDiversity.tsv"
pd.DataFrame(results).to_csv(output, sep="\t", index=False)
print(f"Results saved to {output}")
