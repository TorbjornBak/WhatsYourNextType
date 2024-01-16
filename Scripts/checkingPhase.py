import sys, os

blastfile = sys.argv[1]
allelename = sys.argv[2]
assemblyfile = sys.argv[3]

with open(assemblyfile,'r') as file:
    sequence = ""
    contigname = None
    sequences = dict()
    for line in file:
        if line.startswith(">contig"):
            contigname = line.strip().replace(">","")
            sequences[contigname] = ""
        else:
            sequences[contigname] += line.strip()

notpseudogene = list()
flag = False

with open(blastfile,'r') as file:
    for line in file:
        if line.startswith("Query="):
            contigname = line.split()[1]
            flag = True
        elif line.startswith("HLA") and flag:
            print(line.split())
            print(line.split()[1].split("*")[0])
            if line.split()[1].split("*")[0] in ["A","B","C","DRB1","DQB1","DQA1","DPB1"]:
                notpseudogene.append(contigname)
            flag = False
            
    
print(len(notpseudogene),notpseudogene)
if len(notpseudogene) == 1:
    dir = f"{allelename}_wrong_assembly"
    os.mkdir(dir)
    outfile = open(f"{dir}/assembly.fasta",'w')
    outfile.write(f'>{notpseudogene[0]}\n')
    outfile.write(f'{sequences[notpseudogene[0]]}\n')
    outfile.close()
else:
    dir = f"{allelename}_assembly"
    os.mkdir(dir)
    outfile = open(f"{dir}/{allelename}_assembly.fasta",'w')
    for gene in sequences:
        outfile.write(f'>{gene}\n')
        outfile.write(f'{sequences[gene]}\n')
    outfile.close()
