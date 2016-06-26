import os
import linecache
import sys
import errno
from sys import exit
cutadapt = "~/.local/bin/cutadapt"
tophat2 = "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/tophat/tophat2"
No_of_processors_tophat2 = 8
bam2fastq = "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/bam2fastq-1.1.0/bam2fastq"
bowtie2 = "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/bowtie/bowtie2"
No_of_processors_bowtie2 = 8
samtools = "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/samtools-1.3.1/samtools"
htseq_count = "htseq-count"
picard = "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/picard-tools-2.4.1/picard.jar"
genome_annotation = "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/Genomes/hg19/hg19_genome_annotation/hg19_genes.gtf"
genome_bowtie_index= "/export/cse/NGSDATA/ADLQATAR/FASTQ/2016/Genomes/hg19/fasta/hg19"
adaptor = "GGCCAAGGCG"


def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
def list_files(path):
    # returns a list of names (with extension, without full path) of all files 
    # in folder path
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(name)
    return files
def make_dir_for_each_file(list_files, path):
    dirs = []
    for file_ in list_files:
        dotIndex = file_.index('.')
        dir_ = file_[:dotIndex]
        path_f = path+"//"+dir_
        make_sure_path_exists(path_f)
        dirs.append(path_f)
    return dirs
def remove_adaptor(input_,output_="adaptorTrim.fastq"):
    cmd_re_adap = cutadapt + " -m 16 -b "+adaptor+" -o "+output_+" "+input_
    os.system(cmd_re_adap)
def align_tophat2(adp_trim,Out_dir="tophat_out"):
    cmd_al_tophat = tophat2 + " -p " + str(No_of_processors_tophat2)+" -o "+Out_dir+" --keep-fasta-order --GTF "+genome_annotation+" "+genome_bowtie_index+" "+ adp_trim
    os.system(cmd_al_tophat)
def convert_bam_to_fastq(input_,output_="unmapped.fastq"):
    cmd_conv_to_fastq = bam2fastq + " -o "+output_+" -q "+input_
    os.system(cmd_conv_to_fastq)
def align_bowtie(input_, output_= "unmapped_remap.bam"):
    cmd_al_bow = bowtie2 + " --local --very-sensitive-local -p "+str(No_of_processors_bowtie2) + " --mm -x "+ genome_bowtie_index +" -U " + input_+" | "+ samtools+" view -uhS -F4 - | "+samtools+" sort -o "+output_
    os.system(cmd_al_bow)
def merge_bam_files_picard(input_1, input_2, output_="aligned.bam"):
    cmd_merge_bams = "java -jar "+picard+" MergeSamFiles USE_THREADING=true MSD=true AS=true " +" I=" + input_1+" I="+input_2+" O="+output_
    os.system(cmd_merge_bams)
def convert_bam_sam(input_,output_="aligned.sam"):
    cmd_conv_to_sam = samtools+" view -h -o "+output_+" "+input_
    os.system(cmd_conv_to_sam)
def count_reads(input_,output_="output.sam"):
    cmd_count_reads = htseq_count + " " + input_+" "+ genome_annotation +" -s no -a 10  -o "+output_+" >counts.xls"
    os.system(cmd_count_reads)
def RNA_Pipeline(input_,result_dir=""):
    #remove adaptor sequences
    adap_trim_result = result_dir+"//adaptorTrim.fastq"
    remove_adaptor(input_,adap_trim_result)
    #align with tophat2
    tophat_result = result_dir+"//tophat_output"
    align_tophat2(adap_trim_result, tophat_result)
    #convert unmapped to fastq from bam
    conv_fastq_out = result_dir+"//unmapped.fastq"
    convert_bam_to_fastq(tophat_result+"//unmapped.bam",conv_fastq_out)
    #align the unmapped with bowtie2
    bow_out = result_dir+"//unmapped_remap.bam"
    align_bowtie(conv_fastq_out, bow_out)
    #merge bam files with Picard module MergeSamFiles
    merg_out = result_dir+"//aligned.bam"
    merge_bam_files_picard(tophat_result+"//accepted_hits.bam", bow_out, merg_out)
    #convert aligned to sam from bam
    conv_sam_out = result_dir+"//aligned.sam"
    convert_bam_sam(merg_out, conv_sam_out)
    #Count reads using htseq-count
    count_out = result_dir+"//count.sam"
    count_reads(conv_sam_out, count_out)
           
if __name__ == "__main__":
    if len(sys.argv)<2:
        print "Please add directory to the command line.\nUsage: python RNAPipline.py directory"
    else:
        current_path =os.getcwd()
        dir_in = sys.argv[1]
        try:
            ind_slash = dir_in.index('/')
            dir_in = dir_in[:ind_slash]
        except:
            pass
        error_message = dir_in+" does not exist"
        if os.path.isdir(dir_in):
            path = current_path+"//"+dir_in   
            make_sure_path_exists(path+"_result")
            files = []
            files = list_files(path)
            if len(files)<1:
                print "You do not have any file in "+path+" directory"
            else:
                dirs = make_dir_for_each_file(files, path+"_result")
                for i in range(len(files)):
                    RNA_Pipeline(path+"//"+files[i],dirs[i])
                             
        else:
            print error_message
