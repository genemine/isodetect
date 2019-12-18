# isodetect
A method to detect splice isoforms from long reads

IsoDetect is a method to detect isoforms from TGS data. IsoDetect extracts short feature sequences from annotated isoforms, map feature sequences to long reads, and partition long reads into groups with the shared feature sequences for isoform detection. The major difference of IsoDetect from existing isoform detection methods such as Tofu and Isoseq3 is that it makes use of the feature sequences of annotated isoforms and cluster long reads into groups based on the feature sequences. This characteristic is unique to our method.

# Dependencies
IsoDetect is implemented in Python. Prerequisites of IsoDetect include Bowtie (version 0.12.7), USEARCH (version 9.2.64), GMAP (version 2017-06-20). 

# Run
From the Linux command line, users can use the following code to run IsoDetect:
sh run_isodetect.sh

For any questions, please contact: hongdong@csu.edu.cn
