#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#file opened and read
read_seq=open('seq1.txt','r')
seq1=read_seq.read().strip()

#file closed
read_seq.close()

#length of sequence saved and printed
seqlen1=len(seq1)
print("length of the sequence:", seqlen1)


# In[ ]:


#identifying TATAA motif in the sequence
for i in range(seqlen1-2):
    motif=(seq1[i-2]+seq1[i-1]+seq1[i]+seq1[i+1]+seq1[i+2])
    if motif=="TATAA":
        print("TATAA motif present in seq 1")


# In[ ]:


#variable to store position of motif
tata_pos=0

#identifying TATAA motif position in the sequence
for i in range(seqlen1-2):
    motif=(seq1[i-2]+seq1[i-1]+seq1[i]+seq1[i+1]+seq1[i+2])
    if motif=="TATAA":
        tata_pos=i-2
        print("TATAA motif present at position:", tata_pos)
        
    #alternative message if motif is not found
    elif motif!="TATAA":
        print("TATAA motif not present")


# In[ ]:


#storing Inr motif position and sequences in list
inr_list=[]
inr_seq=[]

#counter for number of Inr motifs present
c=0

nwlen= len(seq1[0:30])

#searching for motifs present within first 30 nts
for i in range(nwlen):
    Y=seq1[i]
    R=seq1[i+1]
    inr_mo=(Y+R)
    if Y=='C' or Y=='T':
        if R=='A' or R=='G':
            c=c+1
            inr_list.append(i)
            inr_seq.append(seq1[i:i+2])
print("number of Inr motifs found:",c)

#printing sequence of Inr motifs
for i in range(len(inr_seq)):
    print(inr_seq[i])


# In[ ]:


#file opened and read
read_seq2=open('seq2.txt','r')
seq2=read_seq2.read().strip()

#file closed
read_seq2.close()

#length of sequence
seqlen2=len(seq2)

#variable to store position of TATAbox
tata_pos=0

#storing position of TSS and sequence in lists
tsseq_list=[]
tss_list=[]

#searching for TATA boxes
for i in range(seqlen2-2):
    tatabox=(seq2[i-2]+seq2[i-1]+seq2[i]+seq2[i+1]+seq2[i+2])
    if tatabox=="TATAA":
        tata_pos=i-2
#looking for TSSs within 24-31nt of the TATA box
        for j in range(tata_pos+24,tata_pos+31):
            Y=seq2[j-1]
            R=seq2[j]
            if Y=='C' or Y=='T':
                if R=='A' or R=='G':
                    tss_list.append(j)
                    tsseq_list.append(seq2[j])
                
print("potential TSSs:", tsseq_list, "at positions:", tss_list, "respectively")




# In[ ]:


#variable to store position of TATAbox
tata_pos=0

#identifying TATAA motif position in the sequence
tsseq_list=[]
tss_list=[]

for i in range(seqlen2-2):
    tatabox=(seq2[i-2]+seq2[i-1]+seq2[i]+seq2[i+1]+seq2[i+2])
    if tatabox=="TATAA":
        tata_pos=i-2
        for j in range(tata_pos+24,tata_pos+31):
            Y=seq2[j-1]
            R=seq2[j]
            inr_mo=(Y+R)
            if inr_mo=='CA' or inr_mo=='CG' or inr_mo=='TA' or inr_mo=='TG':
                tss_list.append(j)
                tsseq_list.append(seq2[j])
                
                #printing core promoter region
                print("core promoter regions:", seq2[j-51:j-1], seq2[j], seq2[j+1:j+51])


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




