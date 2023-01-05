import time

# Start the timer
start_time = time.time()
from Bio import SeqIO
chromolist=[]
with open('input.vcf', 'r') as f:
        for line in f:
            # Skip lines that start with '#' (these are comment lines)
            if line.startswith('#'):     # Split the line into fields
                continue
            fields = line.strip().split('\t')
            if fields[0] not in chromolist:
                chromolist.append(fields[0])
            else: continue
chromodict={}
print(chromolist)
for chromosome in chromolist:
    chromodict[chromosome] = {}
variant ={}
with open('input.vcf', 'r') as f:
    # Iterate through the lines in the file
    for line in f:
        # Skip lines that start with '#' (these are comment lines)
        if line.startswith('#'):     # Split the line into fields
            continue
        fields = line.strip().split('\t')
        sample = fields[9]
        if sample.startswith("0/0"):
            continue
        else:
            # Extract the position and alleles
            chromo = fields[0]
            pos = int(fields[1])-1 
            ref = fields[3]
            alt = fields[4].split(",")
            chromodict[chromo][pos] = (ref, alt)

    print(chromodict)

# Read the reference genome in the Fasta format
file = open('output.fasta', 'w')
with open('reference2.fasta', 'r') as f:
    for record in SeqIO.parse(f, 'fasta'):
        sequence=""
        i=0
        while i < len(record.seq):
            if chromodict[record.id].get(i):
                #print("yes/t")
                #print(variant[i+1])
                sequence+= str(chromodict[record.id][i][1][0])
                #print(len(variant[i][0]))
                i= i+len(chromodict[record.id][i][0])
            else:
                sequence+= str(record.seq)[i]
                i+=1
        file.write(">" + record.id+"\n")
        for i in range(0, len(sequence), 100):
            file.writelines(sequence[i:i+100] + '\n')

        #print(record.seq)
        #print(sequence)
        
        #print(len(sequence),len(record.seq))
file.close()
end_time = time.time()
elapsed_time = end_time - start_time
print(f'Elapsed time: {elapsed_time:.2f} seconds')
