import sys

filename = sys.argv[1]

def removeIntrons(sequence,exons):
    result_sequence = list(sequence)

    # Create a list to represent introns with 'N's
    introns = ['X'] * len(sequence)
    #introns = list(sequence.lower())
    try:
        for exon_start, exon_end in exons:
            exon_start = exon_start.replace("<","")
            # Set intron positions to empty in the intron list
            introns[int(exon_start)-1:int(exon_end)] = result_sequence[int(exon_start)-1:int(exon_end)]
            #print(result_sequence[exon_start:exon_end])
    except TypeError as error:
        print(error)
        print(exons)
        
   
    #print(''.join(introns))
    return ''.join(introns)


def intronRemover():
    sequence = ""
    exons = None
    savefile = "hla_without_introns.fasta"
    with open(savefile,'w') as writefile:
        with open(filename, 'r') as file:
            sequenceFlag = False
            for line in file:
                if line.startswith("AC"):
                    variant = line.split()[1].replace(";","")
                elif line.startswith("DE"):
                    HLAvariant = line.split()[1].replace(",","").replace("HLA-","")
                elif line.startswith("FT") and line.split()[1] == "CDS":
                    exons = line.strip().split()[2].replace("join(",'').replace(")","").split(",")[:-1]
                    #print(exons)
                    exons = [(exon.split("..")[0],exon.split("..")[1]) for exon in exons]
                elif line.startswith("SQ"):
                    sequenceFlag = True
                elif line.startswith("//"):
                    if len(sequence) > 0:
                        sequence = f'{removeIntrons(sequence,exons)}\n'
                        seqlength = len(sequence)
                        variant = f'>HLA:{variant} {HLAvariant} {seqlength} bp\n'
                        writefile.write(variant)
                        writefile.write(sequence)
                    sequence = ""
                    exons = None
                    sequenceFlag = False
                elif sequenceFlag is True:
                    #print(line.strip().split()[:-1])
                    sequenceline = "".join(line.strip().split()[:-1])
                    sequence += sequenceline.upper()


intronRemover()