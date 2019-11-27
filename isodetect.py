# coding:utf-8
# -*- coding: utf-8 -*-
#!/usr/bin/env python
from sys import argv
from collections import Counter
import sys
import string
import re
import os  
raw=argv[1]
fna=argv[2]
task=argv[3]
def fa2fa(org,out):
    fin=open(org,"r")
    fout=open(out,"w")
    seq=""
    while 1:
        line=fin.readline()
        if not line:
            break
        if line.startswith(">"):
            if seq!="":
                fout.write(seq+"\n")
            seq=""
            fout.write(line)
        else:
            line=line.strip("\n")
            seq+=line
    fout.write(seq+"\n")
    fin.close()
    fout.close()
def getft(org,out):
	inputfa = open(org)
	output = open(out, 'w')
	temp="0"
	ftset=[ 'A', 'AA', 'AAA', 'AAC', 'AAG', 'AAT', 'AC', 'ACA', 'ACC', 'ACG', 'ACT', 'AG', 'AGA', 'AGC', 'AGG', 'AGT', 'AT', 'ATA', 'ATC', 'ATG', 'ATT', 'C', 'CA', 'CAA', 'CAC', 'CAG', 'CAT', 'CC', 'CCA', 'CCC', 'CCG', 'CCT', 'CG', 'CGA', 'CGC', 'CGG', 'CGT', 'CT', 'CTA', 'CTC', 'CTG', 'CTT', 'G', 'GA', 'GAA', 'GAC', 'GAG', 'GAT', 'GC', 'GCA', 'GCC', 'GCG', 'GCT', 'GG', 'GGA', 'GGC', 'GGG', 'GGT', 'GT', 'GTA', 'GTC', 'GTG', 'GTT', 'T', 'TA', 'TAA', 'TAC', 'TAG', 'TAT', 'TC', 'TCA', 'TCC', 'TCG', 'TCT', 'TG', 'TGA', 'TGG', 'TGT', 'TT', 'TTA', 'TTC', 'TTG', 'TTT','TGC']
	# output.write('\t'.join(ftset)+'\t'+'len'+'\n')
	for line in inputfa:
		if not line.startswith(">"):
			line=line.replace('\n','')
			for x in ftset:
				output.write(str(line.count(x))+'\t')
			output.write(str(len(line))+'\n')
	output.flush()
	inputfa.close( )
	output.close()
# def

def getnotagFA(fa,align,outf):
	alignidset=set()
	seq=read_fa(fa)
	out=open(outf,'w+')
	seqset=set(seq.keys())
	for line in open(align):
		rawID= str(line.split()[2])
		alignidset.add(rawID)
	for x in (seqset-alignidset):
		out.write(">"+x+"\n"+seq[x]+"\n")
		passset=set(seq.keys())
def read_fna_fa(fa_file):
        seq={}
        fa_out=open("seq",'w+')
        fa_reader=open(fa_file)
        for line in fa_reader:
                if line.startswith('>'):
                        ID=line.replace('>','').replace('\n','')
                        seq[ID]=''
                else:
                        seq[ID]+=line.replace('\n','')
        fa_reader.close()
        for x in seq.keys():
                fa_out.write(x+"\n"+seq[x]+"\n")
                pass        
        return seq
def fna_ee(fna_file,out):
    fna=read_fna_fa(fna_file)
    out=open(out,'w+')
    print("all annotated RNA seqs",len(fna))
    li={}
    for ID in fna.keys():
        if "join" in ID:
            loci=re.sub("\D", " ", ID[ID.index("join"):]).split()
            f=0
            ll=len(loci)-1
            #temp=[]
            x=len(loci)/2 -1
            for i in range(1,ll):
                if i%2!=0:
                    el=int(loci[i])-int(loci[i-1])+1
                    f=f+el
                    seq=fna[ID]
                    #temp.append(seq[(f - 10):(f + 10)])                  
                    #333  print>> out,">"+ID.split()[0]+"."+str(i/2+1)+"\n"+(seq[(f - 10):(f + 10)])
                    li[seq[(f - 10):(f + 10)]]=(ID.split()[0]+"."+str(i/2+1))
        #li[ID]=temp
    ###########out.flush()
    i=0
    for x in li.keys():
        out.write(">"+str(i)+"\n"+x+"\n")
        i+=1
    out.flush()
    out.close()
    #return li
# read_fa input: .FASTA output struct: seq[ID]
def read_fa(fa_file):
	seq={}
	fa_reader=open(fa_file)
	for line in fa_reader:
		if line.startswith('>'):
			ID=line.replace('>','').replace('\n','').split()[0]
			seq[ID]=''
		else:
			seq[ID]+=line.replace('\n','')
	fa_reader.close()
	return seq
def foralignsort(sortfile,seq):
	res=open(sortfile)
	idx={}
	seqx={}
	# t=()
	for line in res:
		line=line.split()
		rawID= str(line[2])
		exseq= str(line[4])	
		idx.setdefault(rawID,set()).add(exseq)
	for x in idx.keys():
		idx[x]=sorted(idx[x])
		pass
	print("seqs:",len(idx))
	Seqset_IDset={}
	for x in idx.keys():
		mm=str(idx[x])
		Seqset_IDset.setdefault(mm,[]).append(x)
	n=len(Seqset_IDset) 
	print("cluster:",n)
	c1=[]
	c2=set()
	outx=open(task+".1c.fa",'w+')
	for x in Seqset_IDset.keys():
		countid=len(Seqset_IDset[x])
		if (countid)>1:
			name=task+".temp."+str(hash(x))
			out=open(name,'w+')
			for countidx in range(0,countid):
				out.write(">"+Seqset_IDset[x][countidx]+"\n"+seq[Seqset_IDset[x][countidx]]+"\n")
			out.flush()
			out.close()
		else :
			outx.write(">"+Seqset_IDset[x][countid-1]+"\n"+seq[Seqset_IDset[x][countid-1]]+"\n")
		pass
# gtf_ee: get the *********
def gtf_ee(fna_file):
	pass

def rename(task,inf,out):
	fa_reader=open(inf)
	outf=open(out,'w+')
	con=0
	for line in fa_reader:
		if line.startswith('>'):
			outf.write(">"+task+".%d\n"% (con))
			con+=1
		else:
			outf.write(line)
	outf.flush()
	outf.close()
	fa_reader.close()
	pass
# def end
print("get fa")
fa2fa(raw,raw+".fa")
fa2fa(fna,fna+".fa")
#print("get ft")
#getft(raw+".fa",raw+".ft")
#getft(fna+".fa",fna+".ft")
print("align and deal")
# 2.1	get feature sequences
	# fna_ee : fna_file to exon-exon junction fasta file
# 2.2	alignment the feature sequences ->

fna_ee(fna+".fa",task+"_eej.fa")
#fna_ee(fna,task+"_eej.fa")
bowtie_build="bowtie-build %s %s -q" % (raw,task+"_db")
bowtie_align="bowtie -f -a -v 2  %s %s %s" % (task+"_db",task+"_eej.fa",task+"_align")
bowtie_sort="sort -k 3 %s > %s " % ( task+"_align",task+"_align_sort")
os.system("echo xx")
os.system(bowtie_build)
os.system(bowtie_align)
os.system(bowtie_sort)
## ***************
#   3.	cluster and get the file 
seq1=read_fa(raw+".fa")
foralignsort(task+"_align_sort",seq1)
# 	3.1 generate the no-fs
print("get no-fs file")
#getnotagFA(raw+".fa",task+"_align_sort",task+".nofs.fa")
# ***************
# #function:exon-exon junction from fna ||  file_to_fa


# # 1.2	generate the ccs->
# 2.1	get feature sequences
# 2.2	alignment the feature sequences ->
 #  3.	cluster and get the file
	# 3.1.deal the clusters
# ***************
print("cluster")
clu="for i in $(ls %s.temp.*); do sh clu.sh $i; done" % (task)
os.system(clu)
merge="cat uclu.%s.temp.*.fax > %s.cons.fa" % (task,task)
os.system(merge)
rm="rm uclu.%s.temp.*.fax"% (task)
os.system(rm)
rm="rm %s.temp.*"% (task)
os.system(rm)
print("evaluate")
#getft(task+".cons.fa",task+".cons.fa.ft")
# ***************
# out: $task.cons.fasta
# 	3.2.deal the non-FS
# 评估，cluster ，gmap
sortNOfs="nohup usearch9.2.64_i86linux32 -sortbylength %s.%s.fa -fastaout %s.%s.sort.fa &" % (task,"nofs",task,"nofs")
cluNOfs="nohup usearch9.2.64_i86linux32 -cluster_fast %s.%s.sort.fa  -consout %s.%s.clu.fa -id 0.75 &" % (task,"nofs",task,"nofs")
print("evaluate nofs")
rmno="rm %s.nofs.sort.fa &" % (task)
os.system(sortNOfs)
os.system(cluNOfs)
os.system(rmno)
# # out: $task.nofs.clu.fa
sort1cfs="usearch9.2.64_i86linux32 -sortbylength %s.%s.fa -fastaout %s.%s.sort.fa" % (task,"1c",task,"1c")
clu1cs="usearch9.2.64_i86linux32 -cluster_fast %s.%s.sort.fa  -consout %s.%s.clu.fa -id 0.75 &" % (task,"1c",task,"1c")
print("evaluate 1c")
rm1c="rm %s.%s.sort.fa &" % (task,"1c")
os.system(sort1cfs)
os.system(clu1cs)
os.system(rm1c)
# 	3.3.deal the 1-cluster-s
# out : $task.1c.clu.fa 
# 4.	evaluate the isoforms
mergefile="cat %s %s %s >%s.final.fa" %(task+".cons.fa",task+".1c.clu.fa",task+".nofs.clu.fa",task)
os.system(mergefile)
# 5.	stastic the results
rename(task,task+".cons.fa",task+".cons.rename.fa")
rename(task,task+".1c.clu.fa",task+".1c.rename.fa")
rename(task,task+".nofs.clu.fa",task+".nofs.rename.fa")
fa2fa(task+".final.fa",task+".out.fa")
rename(   task,task+".out.fa",task+".finall.fa" )
os.system("rm %s.final.fa %s.out.fa" %(task,task))
os.system("rm %s.output.fa %s.out.fa" %(task,task))
#Rscript train.R $1.ft $2.ft $3.model
#print("train model")
#trainmodel="Rscript train.R %s %s %s" % (fna+".ft",raw+".ft",task+".model")
# os.system(trainmodel)
#print("evaluate the reads")
#eva="Rscript eva.R %s %s %s" % (fna+".ft","human",task+".ft",task+".p")
#os.system(eva)


