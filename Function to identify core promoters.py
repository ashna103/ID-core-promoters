#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#function to identify promoters in given sequence file
def idcorep(seqfile):
    
    #open, read and close file
    read_seq=open(seqfile,'r')
    dnaseq=read_seq.read().strip()
    read_seq.close()
    
    #save sequence length and other variables
    seqlen=len(dnaseq)
    tata_pos=0
    tsseq_list=[]
    tss_list=[]
    
    #to prevent repetition of the print command 
    print("core promoter regions:")
    
    #identify tata box
    for i in range(seqlen-2):
        tatabox=(dnaseq[i-2]+dnaseq[i-1]+dnaseq[i]+dnaseq[i+1]+dnaseq[i+2])
        if tatabox=="TATAA":
            tata_pos=i-2
            #identify Inr region, TSS and core promoter region
            for j in range(tata_pos+24,tata_pos+31):
                Y=dnaseq[j-1]
                R=dnaseq[j]
                inr_mo=(Y+R)
                if inr_mo=='CA' or inr_mo=='CG' or inr_mo=='TA' or inr_mo=='TG':
                    tsseq_list.append(dnaseq[j])
                    print(dnaseq[j-51:j-1], dnaseq[j], dnaseq[j+1:j+51])
                    
                    
    return tsseq_list

#result
idcorep('seq2.txt')



    

