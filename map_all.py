#BWAcycle tool
#Read mapping for influenza sequencing on Illumina
#Requires Matplotlib, Weblogo, Biopython, Pillow
import sys
import os
import math
import shutil
import statistics
import operator
from datetime import datetime
import base64
import argparse
import zipfile
import matplotlib.pyplot as ppl 
from Bio import SeqIO
from PIL import Image

#Path to required files
pth = os.path.dirname(os.path.abspath(__file__))
trimmomatic_path = pth+os.sep+'stuff'+os.sep+'trimmomatic'+os.sep
fastqc_path = pth+os.sep+'stuff'+os.sep+'FastQC'
bwa_path = pth+os.sep+'stuff'+os.sep+'bwa'+os.sep
#bwa_path = pth+os.sep+'stuff'+os.sep
samtools_path = pth+os.sep+'stuff'+os.sep
bam_readcount_path = pth+os.sep+'stuff'+os.sep
minimap_path = pth+os.sep+'stuff'+os.sep
centrifuge_path = pth+os.sep+'stuff'+os.sep
picard_path = pth+os.sep+'stuff'+os.sep
freebayes_path = pth+os.sep+'stuff'+os.sep
vardict_path = pth+os.sep+'stuff'+os.sep+'VarDict'+os.sep+'lib'+os.sep
reference_path = pth+os.sep+'Reference'+os.sep



def test():
	print("\nTest\n")

def copy_data(sample_number, plate_number, data_type):
	#Создаем папку по лабораторному номеру образца, копируем риды в нее
	#Create sample folder, copy data to sample folder
	if not os.path.isdir(sample_number):
		os.mkdir(sample_number)
	if data_type == 'illumina':
		dist1 = sample_number+os.sep+sample_number+'_R1.fastq.gz'
		src1 = sample_number+'_S'+plate_number+'_L001_R1_001.fastq.gz'
		dist2 = sample_number+os.sep+sample_number+'_R2.fastq.gz'
		src2 = sample_number+'_S'+plate_number+'_L001_R2_001.fastq.gz'
		shutil.copyfile(src1, dist1)
		shutil.copyfile(src2, dist2)
	if data_type == 'nanopore':
		dist1 = sample_number+os.sep+sample_number+'_unfiltered.fastq.gz'
		src1 = sample_number+'_'+plate_number+'.fastq.gz'
		shutil.copyfile(src1, dist1)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	print("\nCopy data to new folder for sample "+sample_number+"\n")

def trim_data(sample_number, trimmomatic_path):
	#Запускаем Trimmomatic для обрезки входных данных по качеству
	#Use Trimmomatic for quality trimming
	if data_type == 'illumina':
		trimmomatic_files = sample_number+os.sep+sample_number+'_R1.fastq.gz '+sample_number+os.sep+sample_number+'_R2.fastq.gz '+sample_number+os.sep+sample_number+'_fp.fq '+sample_number+os.sep+sample_number+'_fu.fq '+sample_number+os.sep+sample_number+'_rp.fq '+sample_number+os.sep+sample_number+'_ru.fq'
		src = trimmomatic_path+'adapters'+os.sep+'NexteraPE-PE.fa'
		dest = 'NexteraPE-PE.fa'
		shutil.copyfile(src, dest)
		print("\nRun Trimmomatic for quality clipping of data  for sample "+sample_number+". Further output from Trimmomatic:\n")
		os.system('java -jar '+trimmomatic_path+'trimmomatic-0.38.jar PE -threads '+str(os.cpu_count())+' -phred33 '+trimmomatic_files+' ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 HEADCROP:15')
		os.remove('NexteraPE-PE.fa')
	elif data_type == 'nanopore':
		print("\nRun Trimmomatic for quality clipping of data  for sample "+sample_number+". Further output from Trimmomatic:\n")
		trimmomatic_files = sample_number+os.sep+sample_number+'_unfiltered.fastq.gz '+sample_number+os.sep+sample_number+'.fastq'
		os.system('java -jar '+trimmomatic_path+'trimmomatic-0.38.jar SE -threads '+str(os.cpu_count())+' -phred33 '+trimmomatic_files+' LEADING:3 TRAILING:3 MINLEN:30 HEADCROP:20')
	print("\nTrimmomatic finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def conc_data(sample_number):
	#Собираем 4 файла выдачи Trimmomatic в один FASTQ, считаем количество строк в полученном файле
	#Concatenate 4 Trimmomatic output files, count read number
	print("\nConcatenate clipped data to one file for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	with open(sample_number+os.sep+sample_number+'_fp.fq', 'r') as f:
		fp = f.readlines()
	with open(sample_number+os.sep+sample_number+'_fu.fq', 'r') as f:
		fu = f.readlines()
	with open(sample_number+os.sep+sample_number+'_rp.fq', 'r') as f:
		rp = f.readlines()
	with open(sample_number+os.sep+sample_number+'_ru.fq', 'r') as f:
		ru = f.readlines()
	reads = fp + fu + rp + ru
	total_str = len(reads)//4
	print("\nTotal reads obtained: "+str(len(reads))+"\n")
	with open(sample_number+os.sep+sample_number+'.fastq', 'w') as f:
		for item in reads:
			f.write(item)
	print("\nRemove temporary files for sample "+sample_number+"\n")
	os.remove(sample_number+os.sep+sample_number+'_fp.fq')
	os.remove(sample_number+os.sep+sample_number+'_fu.fq')
	os.remove(sample_number+os.sep+sample_number+'_rp.fq')
	os.remove(sample_number+os.sep+sample_number+'_ru.fq')
	return total_str

def fastqc_data(inp_file):
	#Запускаем FASTQC для оценки качества в обработанном FASTQ
	#Run FASTQC for quality report
	print("\nRun FastQC to obtain quality report for file "+inp_file+". Further output from FastQC:\n")
	if os.sep == "/":
		fastqc_par = 'java -Xmx512m -classpath '+fastqc_path+':'+fastqc_path+'/sam-1.103.jar:'+fastqc_path+'/jbzip2-0.9.jar -Djava.awt.headless=true uk.ac.babraham.FastQC.FastQCApplication '+inp_file
	else:
		fastqc_par = 'java -Xmx512m -classpath '+fastqc_path+';'+fastqc_path+'/sam-1.103.jar;'+fastqc_path+'/jbzip2-0.9.jar -Djava.awt.headless=true uk.ac.babraham.FastQC.FastQCApplication '+inp_file
	os.system(fastqc_par)
	print("\nFastQC finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def index_bwa(ref_file):
	#Делаем индекс для работы BWA из ref_file
	#Prepare index for BWA alignment
	print("\nMake index from reference file "+ref_file+". Further output from BWA index:\n")
	bwa_cmd = bwa_path+'bwa index '+ref_file
	os.system(bwa_cmd)
	print("\nBWA index finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def run_bwa(ref_file, sample_number, step):
	#Проверяем существование индекса
	#Check if index already exists
	if not os.path.isfile(ref_file+'.bwt'):
		index_bwa(ref_file)
	#Запускаем выравнивание BWA
	#Run BWA alingment
	print("\nAlign reads by BWA mem for sample "+sample_number+". Further output from BWA mem:\n")
	if data_type == "illumina":
		if (step == '' or step == 'v'):
			bwa_cmd = bwa_path+'bwa mem -aM -R "@RG\\tID:'+sample_number+'\\tSM:'+sample_number+'\\tPL:ILLUMINA\\tLB:library1\\t" -t '+str(os.cpu_count())+' '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+sample_number+'_'+step+'.sam'
		else:
			bwa_cmd = bwa_path+'bwa mem -aM -R "@RG\\tID:'+sample_number+'\\tSM:'+sample_number+'\\tPL:ILLUMINA\\tLB:library1\\t" -t '+str(os.cpu_count())+' -k 13  -B 2 -O [4,4] -E [1,1] -L [3,3] -U 1 '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+sample_number+'_'+step+'.sam'
	elif data_type == "nanopore":
		pl = "NANOPORE"
		if (step == '' or step == 'v'):
			bwa_cmd = bwa_path+'bwa mem -aM -R "@RG\\tID:'+sample_number+'\\tSM:'+sample_number+'\\tPL:NANOPORE\\tLB:library1\\t" -t '+str(os.cpu_count())+' -k 14 -W 20 -r 10 -A 1 -B 1 -O 1 -E 1 -L 0 '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+sample_number+'_'+step+'.sam'
		else:
			bwa_cmd = bwa_path+'bwa mem -aM -R "@RG\\tID:'+sample_number+'\\tSM:'+sample_number+'\\tPL:NANOPORE\\tLB:library1\\t" -t '+str(os.cpu_count())+' -k 13  -B 2 -O [2,2] -L [1,1] -E [1,1] -U 1 '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+sample_number+'_'+step+'.sam'
	os.system(bwa_cmd)
	print("\nBWA mem alignment finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	
def copy_result(sample_number, gene, subtype):
	#Копирование результатов FASTA, SVG и PNG в папку Result
	#Copy resulting files
	print("\nCopy resulting files to Result folder\n")
	if not os.path.isdir("result"+os.sep+subtype):
		os.mkdir("result"+os.sep+subtype)
	if not os.path.isdir("result"+os.sep+subtype+os.sep+gene):
		os.mkdir("result"+os.sep+subtype+os.sep+gene)
	if not os.path.isdir("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"):
		os.mkdir("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta")
	if not os.path.isdir("result"+os.sep+subtype+os.sep+gene+os.sep+"svg"):
		os.mkdir("result"+os.sep+subtype+os.sep+gene+os.sep+"svg")
	if not os.path.isdir("result"+os.sep+subtype+os.sep+gene+os.sep+"png"):
		os.mkdir("result"+os.sep+subtype+os.sep+gene+os.sep+"png")
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+".fasta"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+sample_number+"_"+gene+".fasta"
	shutil.copyfile(src, dest)
	for graph in ("Shannon", "Coverage", "BaseQuality", "MapQuality", "SeqDiversity"):
		src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_"+graph+".svg"
		dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"svg"+os.sep+sample_number+"_"+gene+"_"+graph+".svg"
		shutil.copyfile(src, dest)
	for graph in ("Shannon", "Coverage", "BaseQuality", "MapQuality", "SeqDiversity"):
		src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_"+graph+".png"
		dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_"+graph+".png"
		shutil.copyfile(src, dest)
	with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+sample_number+"_"+gene+".fasta", 'r') as inp:
		fa = inp.readlines()
	ft = []
	if os.path.isfile("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta"):
		with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta", 'r') as total:
			ft = total.readlines()
	ft += fa
	with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta", 'w') as outp:
		for item in ft:
			outp.write(item)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def html_report(report_file, report_arr, sample_number):
	#Генерируем HTML-отчет
	#Generate HTML-report
	fastqc_zip = zipfile.ZipFile(sample_number+os.sep+sample_number+'_fastqc.zip', 'r')
	for graph in ('per_base_quality', 'per_base_sequence_content', 'sequence_length_distribution'):
		source = fastqc_zip.open(sample_number+'_fastqc/Images/'+graph+'.png')
		target = open(sample_number+os.sep+sample_number+'_'+graph+'.png', "wb")
		shutil.copyfileobj(source, target)
		source.close()
		target.close()
	for graph in ('per_base_quality', 'per_base_sequence_content', 'sequence_length_distribution', 'read_distribution'):
		qfile = open(sample_number+os.sep+sample_number+'_'+graph+'.png', 'rb')
		qfiledata = qfile.read()
		qdata = base64.b64encode(qfiledata)
		qstring = qdata.decode()
		report_arr[sample_number][graph] = qstring
		qfile.close()
	for gene in report_arr[sample_number]['gene_list']:
		wbl = Image.open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_logo.eps')
		wbl.save(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Weblogo.png', 'PNG')
		for graph in ('BaseQuality', 'MapQuality', 'SeqDiversity', 'Shannon', 'Coverage', 'Weblogo'):
			qfile = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+graph+'.png', 'rb')
			qfiledata = qfile.read()
			qdata = base64.b64encode(qfiledata)
			qstring = qdata.decode()
			report_arr[sample_number][gene][graph+'_graph'] = qstring
			qfile.close()
	print("\nGenerating HTML report\n")
	html_page = '<!DOCTYPE html>\n<html>\n<head>\n<meta charset="UTF-8">\n<title>BWAcycle analysis report</title></head>\n<body>\n<h1 id =top><b>BWAcycle sequence analysis report</b></h1>\n<br/>\n<table width = "100%" border = 1 cellspacing=0 cellpading=0 bordercolor = "black">\n'
	html_page += '<tr>\n<td>Sample name</td><td>Sample number</td><td>Subtype</td><td>Genes</td><td></td><td></td></tr>\n'
	for sample in report_arr['sample_list']:
		html_page += '<tr>\n<td><a href="#'+sample+'">'+report_arr[sample]['sample_name']+'</a></td><td>'+sample+'</td><td>'+report_arr[sample]['subtype']+'</td>\n'
		for gene in report_arr[sample]['gene_list']:
			if report_arr[sample][gene]['problems'] < 10:
				html_page += '<td bgcolor ="green"><a href="#'+sample+'_'+gene+'">'+gene+'</a></td>'
			else:
				html_page += '<td bgcolor ="yellow"><a href="#'+sample+'_'+gene+'">'+gene+'</a></td>'
		html_page += '</tr>\n'
	html_page += '</table><br/><br/><br/>\n<table>\n'
	for sample in report_arr['sample_list']:
		html_page += '<tr id='+sample+'>\n<td><h2>'+report_arr[sample]['sample_name']+'</h2><br/><br/>Sample name: '+report_arr[sample]['sample_name']+'<br/>Sample number: '+sample+'<br/>Number in samplesheet: '+report_arr[sample]['plate_number']+'<br/>Index 1: '+report_arr[sample]['index1']+'<br/>Index 2: '+report_arr[sample]['index2']+'<br/>Total read number: '+str(report_arr[sample]['read_number'])+'<br/>Segment selection results: '+str(report_arr[sample]['seg_stats'])+'<br/><br/>\n<h3>Assembled genes: </h3>'
		for gene in report_arr[sample]['gene_list']:
			html_page += '<a href="#'+sample+'_'+gene+'">'+gene+'</a><br/>'
		html_page += '<br/><br/><h3>FASTQC report for quality trimmed data</h3><br/><br/><h4 >Per base quality distribution graph</h4><br/><img src="data:image/png;base64,'+report_arr[sample]['per_base_quality']+'"><br/><h4 >Per base sequence content graph</h4><br/><img src="data:image/png;base64,'+report_arr[sample]['per_base_sequence_content']+'"><br/><h4 >Sequence length distribution graph</h4><br/><img src="data:image/png;base64,'+report_arr[sample]['sequence_length_distribution']+'"><br/><h4 >Per gene read distribution graph</h4><br/><img src="data:image/png;base64,'+report_arr[sample]['read_distribution']+'"><br/>\n'
		for gene in report_arr[sample]['gene_list']:
			html_page += '<h3 id='+sample+'_'+gene+'>'+gene+'</h3><br/><br/>Gene name: '+gene+'<br/>Number of mapped reads: '+str(report_arr[sample][gene]['mapped_reads'])+'<br/>Gene length: '+str(report_arr[sample][gene]['length'])+'<br/>Median coverage: '+str(report_arr[sample][gene]['coverage'])+'<br/>Number of potentially problematic positions (coverage <10 or more than 50 percent indels): '+str(report_arr[sample][gene]['problems'])+' '+str(report_arr[sample][gene]['problems_list'].keys())+' <br/>Median base quality (PHRED): '+str(report_arr[sample][gene]['basequality'])+'<br/>Median mapping quality (PHRED): '+str(report_arr[sample][gene]['mapquality'])+'<br/>Median sequence diversity (Nelson et al., doi:10.1016/j.meegid.2014.11.026): '+str(report_arr[sample][gene]['seqdiversity'])+'<br/>Mean Shannon enthropy: '+str(report_arr[sample][gene]['shannon'])+'<br/><br/>'
			html_page += '<h4 >SNP analysis results: </h4><br/>'
			html_page += 'Naive SNP calling results: <br/>'
			if 'naive_SNPS' in report_arr[sample][gene].keys() and len(report_arr[sample][gene]['naive_SNPS']) > 0:
				for el in report_arr[sample][gene]['naive_SNPS'].keys():
					html_page += '&rarr; '+str(el)+': base frequencies (pA, pC, pG, pT): '+str(report_arr[sample][gene]['naive_SNPS'][el])+'<br/>'
			html_page += 'Shannon enthropy SNP calling results: <br/>'
			html_page += '&rArr; Positions with Shannon enthropy more than mean+10stdev: <br/>'
			if 'Shannon10_SNPS' in report_arr[sample][gene].keys() and len(report_arr[sample][gene]['Shannon10_SNPS']) > 0:
				for el in report_arr[sample][gene]['Shannon10_SNPS'].keys():
					html_page += '&rarr; '+str(el)+': Shannon enthropy value: '+str(report_arr[sample][gene]['Shannon10_SNPS'][el])+'<br/>'
			html_page += '&rArr; Positions with Shannon enthropy more than mean+5stdev: <br/>'
			if 'Shannon5_SNPS' in report_arr[sample][gene].keys() and len(report_arr[sample][gene]['Shannon5_SNPS']) > 0:
				for el in report_arr[sample][gene]['Shannon5_SNPS'].keys():
					html_page += '&rarr; '+str(el)+': Shannon enthropy value: '+str(report_arr[sample][gene]['Shannon5_SNPS'][el])+'<br/>'
			html_page += 'FreeBayes SNP calling results: <br/>'
			if report_arr[sample][gene]['length'] != report_arr[sample][gene]['length_rmdup']:
				html_page += 'Warning! Length of sequence has changed after removing of PCR duplicates!<br/>'
			if 'FB_SNPS' in report_arr[sample][gene].keys() and len(report_arr[sample][gene]['FB_SNPS']) > 0:
				for el in report_arr[sample][gene]['FB_SNPS'].keys():
					html_page += '&rarr; '+str(el)+': '+str(report_arr[sample][gene]['FB_SNPS'][el])+'<br/>'
			html_page += '<h4 >Coverage graph for '+gene+' of '+report_arr[sample]['sample_name']+'</h4><br/><img src="data:image/png;base64,'+report_arr[sample][gene]['Coverage_graph']+'"><br/>'
			html_page += '<h4 >Base quality graph for '+gene+' of '+report_arr[sample]['sample_name']+'</h4><br/><img src="data:image/png;base64,'+report_arr[sample][gene]['BaseQuality_graph']+'"><br/>'
			html_page += '<h4 >Mapping quality graph for '+gene+' of '+report_arr[sample]['sample_name']+'</h4><br/><img src="data:image/png;base64,'+report_arr[sample][gene]['MapQuality_graph']+'"><br/>'
			html_page += '<h4 >Sequence diversity graph for '+gene+' of '+report_arr[sample]['sample_name']+'</h4><br/><img src="data:image/png;base64,'+report_arr[sample][gene]['SeqDiversity_graph']+'"><br/>'
			html_page += '<h4 >Shannon enthropy graph for '+gene+' of '+report_arr[sample]['sample_name']+'</h4><br/><img src="data:image/png;base64,'+report_arr[sample][gene]['Shannon_graph']+'"><br/>'
			html_page += '<h4 >Weblogo graph for '+gene+' of '+report_arr[sample]['sample_name']+'</h4><br/><img src="data:image/png;base64,'+report_arr[sample][gene]['Weblogo_graph']+'"><br/>'
			html_page += '<a href="#'+sample+'">Back to sample</a><br/>\n' 
		html_page += '<a href="#top">Back to top</a><br/></td></tr>\n'
	html_page += '</table></body>\n</html>'
	print("\nSaving HTML report\n")
	with open(report_file, 'w') as rpf:
		rpf.write(html_page)

def run_samtools2(sample_number, ref_file, step, report_arr):
	#Преобразуем полученное от BWA выравнивание из SAM в BAM, сортируем и индексируем BAM-файл, получаем статистику по выравниванию, удаляем из BAM невыровненные риды
	#Converting alignment, obtaining statistics data
	print("\nConvert data to BAM-file with Samtools for sample "+sample_number+". Further output from Samtools:\n")
	samtools_cmd1 = samtools_path+'samtools faidx '+ ref_file
	samtools_cmd2 = samtools_path+'samtools view -b  -t '+ref_file+'.fai '+sample_number+os.sep+sample_number+'_'+step+'.sam | '+samtools_path+'samtools sort -l 9 -@ '+str(os.cpu_count())+' -o '+sample_number+os.sep+sample_number+'_'+step+'.bam'
	samtools_cmd3 = samtools_path+'samtools index '+sample_number+os.sep+sample_number+'_'+step+'.bam'
	samtools_cmd4 = samtools_path+'samtools flagstat '+sample_number+os.sep+sample_number+'_'+step+'.bam > '+sample_number+os.sep+sample_number+'_'+step+'_flagstat.txt'
	samtools_cmd5 = samtools_path+'samtools idxstats '+sample_number+os.sep+sample_number+'_'+step+'.bam > '+sample_number+os.sep+sample_number+'_'+step+'_idxstats.txt'
	os.system(samtools_cmd1)
	os.system(samtools_cmd2)
	os.system(samtools_cmd3)   #-F 0x4
	os.system(samtools_cmd4)
	os.system(samtools_cmd5)
	print("\nSamtools finished\n")
	print("\nRemove temporary files for sample "+sample_number+"\n")
	os.remove(sample_number+os.sep+sample_number+'_'+step+'.sam')
	print("\nSamtools finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	
def get_cons2(sample_number, ref_file, step):
	#Получаем консенсусную последовательность и статистику по позициям
	#Get consensus sequence
	print("\nExtracting base counts from BAM-file with Bam-readcount for sample "+sample_number+". Further output from Bam-readcount:\n")
	bamrc_cmd = bam_readcount_path+'bam-readcount -w 0 -f '+ref_file+' '+sample_number+os.sep+sample_number+'_'+step+'.bam > '+sample_number+os.sep+sample_number+'_'+step+'_readcount.txt'
	os.system(bamrc_cmd)
	print("\nBam-readcount finished\n")
	with open(sample_number+os.sep+sample_number+'_'+step+'_readcount.txt', 'r') as f:
		rc_str_arr = f.readlines()
	rc = {}
	rc['genes'] = []
	for el in rc_str_arr:
		rc_arr = el.split()
		tmp_arr = rc_arr[0].split('_')
		if len(tmp_arr) < 3:
			seg = rc_arr[0]
		else:
			seg = tmp_arr[len(tmp_arr)-2] +'_'+ tmp_arr[len(tmp_arr)-1]
		if seg not in rc['genes']:
			rc['genes'].append(seg)
			rc[seg] = {}
		rc[seg][int(rc_arr[1])] = {"A":rc_arr[5].split(':'), "C":rc_arr[6].split(':'), "G":rc_arr[7].split(':'), "T":rc_arr[8].split(':'), "INDEL":[], "COV":int(rc_arr[3]), "REF":rc_arr[2]}
		if len(rc_arr) > 10:
			ind_arr = []
			for elt in rc_arr[10:]:
				ind_arr.append(elt.split(':'))
			rc[seg][int(rc_arr[1])]["INDEL"] = ind_arr
	for gene in rc['genes']:
		print("\nCalculating consensus for sample "+sample_number+" and gene "+gene+"\n")
		con = ""
		dels = [0,0]
		probl_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_problems'+step+'.txt', 'w')
		ks = list(rc[gene].keys())
		ks.sort()
		for i in ks:
			if rc[gene][i]["COV"] < 10:
				print("\nCoverage in gene "+gene+" is lower than 10 in position "+str(i)+". Please review basecount data.\n")
				probl_file.write(str(i)+"\tCoverage in gene "+gene+" is lower than 10 in position "+str(i)+". Please review basecount data.\n")
			if dels[0] > 0 and dels[1] > rc[gene][i]["COV"]:
				dels[0]-=1
			else:
				pos = rc[gene][i]["REF"]
				if int(rc[gene][i]["A"][1]) > int(rc[gene][i]["C"][1]) and int(rc[gene][i]["A"][1]) > int(rc[gene][i]["G"][1]) and int(rc[gene][i]["A"][1]) > int(rc[gene][i]["T"][1]):
					pos = "A"
				if int(rc[gene][i]["C"][1]) > int(rc[gene][i]["A"][1]) and int(rc[gene][i]["C"][1]) > int(rc[gene][i]["G"][1]) and int(rc[gene][i]["C"][1]) > int(rc[gene][i]["T"][1]):
					pos = "C"
				if int(rc[gene][i]["G"][1]) > int(rc[gene][i]["C"][1]) and int(rc[gene][i]["G"][1]) > int(rc[gene][i]["A"][1]) and int(rc[gene][i]["G"][1]) > int(rc[gene][i]["T"][1]):
					pos = "G"
				if int(rc[gene][i]["T"][1]) > int(rc[gene][i]["C"][1]) and int(rc[gene][i]["T"][1]) > int(rc[gene][i]["G"][1]) and int(rc[gene][i]["T"][1]) > int(rc[gene][i]["A"][1]):
					pos = "T"
				if len(rc[gene][i]["INDEL"]) > 0:
					indel_cov = [0, int(rc[gene][i]["INDEL"][0][1])]
					for j in range(len(rc[gene][i]["INDEL"])):
						if int(rc[gene][i]["INDEL"][j][1]) > indel_cov[1]:
							indel_cov = [j, int(rc[gene][i]["INDEL"][j][1])]
					if rc[gene][i]["INDEL"][indel_cov[0]][0][0] == "+" and float(indel_cov[1])/float(rc[gene][i]["COV"]) > 0.5:
						pos += rc[gene][i]["INDEL"][indel_cov[0]][0][1:]
						print("\nFound more than 50 percent of insertions in gene "+gene+" in position "+str(i)+". Please review basecount data.\n")
						probl_file.write(str(i)+"\tInsertions present in more than 50 percent of reads in position "+str(i)+". Please review basecount data.\n")
					if rc[gene][i]["INDEL"][indel_cov[0]][0][0] == "-" and float(indel_cov[1])/float(rc[gene][i]["COV"]) > 0.5:
						dels = [len(rc[gene][i]["INDEL"][indel_cov[0]][0][1:]), int(rc[gene][i]["INDEL"][indel_cov[0]][1])]
						print("\nFound more than 50 percent of deletions in gene "+gene+" in position "+str(i)+". Please review basecount data.\n")
						probl_file.write(str(i)+"\tDeletions present in more than 50 percent of reads in position "+str(i)+". Please review basecount data.\n")
				con += pos
		ou = []
		if step != '':
			if os.path.isfile(sample_number+os.sep+sample_number+'_'+step+'.fasta'):
				with open(sample_number+os.sep+sample_number+'_'+step+'.fasta', 'r') as out_file:
					ou = out_file.readlines()
			with open(sample_number+os.sep+sample_number+'_'+step+'.fasta', 'w') as out_file:
				if len(ou) > 0:
					for st in ou:
						out_file.write(st)
				out_file.write('>'+sample_name+'_'+gene+'\n'+con+'\n')	
		else:
			if os.path.isfile(sample_number+os.sep+sample_number+'.fasta'):
				with open(sample_number+os.sep+sample_number+'.fasta', 'r') as out_file:
					ou = out_file.readlines()
			with open(sample_number+os.sep+sample_number+'.fasta', 'w') as out_file:
				if len(ou) > 0:
					for st in ou:
						out_file.write(st)
				out_file.write('>'+sample_name+'_'+gene+'\n'+con+'\n')
			with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.fasta', 'w') as out_file:
				out_file.write('>'+sample_name+'_'+gene+'\n'+con+'\n')
		rc[gene]["CON"] = con
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return rc

def get_stats(rc_arr, sample_number, sample_name, report_arr):
	#Собираем статистику, строим графики
	#Obtain statistics, plot graphs
	data_names = ['Unmapped']
	data_values = [report_arr[sample_number]['unmapped_reads']]
	for gene in rc_arr['genes']:
		data_names.append(gene)
		data_values.append(report_arr[sample_number][gene]['mapped_reads'])
	ppl.pie(data_values, labels = data_names, autopct='%1.2f%%', radius = 1.1, explode = [0.15] + [0 for _ in range(len(data_names) - 1)] )
	ppl.title('Reads distribution for sample: '+sample_name+'\nTotal: '+str(report_arr[sample_number]['read_number'])+' reads')
	ppl.axis('equal')
	ppl.savefig(sample_number+os.sep+sample_number+'_read_distribution.svg')
	ppl.savefig(sample_number+os.sep+sample_number+'_read_distribution.png')
	ppl.clf()
	ppl.cla()
	with open(sample_number+os.sep+sample_number+'_rmdup_depth.txt', 'r') as dpth:
		dpth_arr = dpth.readlines()
	dp = {}
	for el in dpth_arr:
		el_sp = el.split("\t")
		if el_sp[0] not in dp.keys():
			dp[el_sp[0]] = [[], []]
		dp[el_sp[0]][0].append(int(el_sp[1]))
		dp[el_sp[0]][1].append(int(el_sp[2]))
	#Получение консенсусной последовательности, таблиц и графиков из массива частот встречаемости
	#Parsing nucleotide frequency array
	ans = {}
	for gene in rc_arr['genes']:
		cov_arr = rc_arr[gene]
		table = []
		transfac = []
		print("\nCalculating statistics for sample "+sample_number+" and gene "+gene+"\n")
		table.append("Position\tA\tC\tG\tT\tpA\tpC\tpG\tpT\tCoverage\tConsensus\tBaseQuality\tMapQuality\tDiversity\tShannon\tA_SNP\tC_SNP\tG_SNP\tT_SNP\tShannon5_SNPS\tShannon10_SNPS")
		transfac.append("ID\tID\nBF\tBF\nP0\tA\tC\tG\tT")
		cons = cov_arr.pop("CON")
		print("\nLength of gene "+gene+" of sample "+sample_number+" after final run: "+str(len(cons))+"\n")
		Shan = []
		Cov = []
		Div = []
		bq_arr = []
		mq_arr = []
		sln = len(cov_arr)
		prob_count = 0
		ks = list(cov_arr.keys())
		ks.sort()
		for j in ks:
			Acount = int(cov_arr[j]["A"][1])
			Gcount = int(cov_arr[j]["G"][1])
			Ccount = int(cov_arr[j]["C"][1])
			Tcount = int(cov_arr[j]["T"][1])
			CovJ = cov_arr[j]["COV"]
			cn = {"A":Acount, "C":Ccount, "G":Gcount, "T":Tcount}
			pos = max(cn, key=cn.get)
			bq_arr.append(float(cov_arr[j][pos][3]))
			mq_arr.append(float(cov_arr[j][pos][2]))
	#Расчет энтропии Шеннона
			PA = 0.0
			PG = 0.0
			PC = 0.0
			PT = 0.0
			if float(CovJ) != 0.0:
				PA = float(Acount)/float(CovJ)
				PG = float(Gcount)/float(CovJ)
				PC = float(Ccount)/float(CovJ)
				PT = float(Tcount)/float(CovJ)
			AENT = 0.0
			GENT = 0.0
			CENT = 0.0
			TENT = 0.0
			if PA != 0.0:
				AENT = -1*PA*math.log2(PA)
			if PG != 0.0:
				GENT = -1*PG*math.log2(PG)
			if PC != 0.0:
				CENT = -1*PC*math.log2(PC)
			if PT != 0.0:
				TENT = -1*PT*math.log2(PT)
			ShanJ = AENT + GENT + CENT + TENT
			DivJ = 0
			if CovJ > 10:
				DivJ = (Acount*Ccount + Acount*Gcount + Acount*Tcount + Ccount*Gcount + Ccount*Tcount + Gcount*Tcount)/((CovJ*CovJ-CovJ)/2)
	#Вывод результатов в таблицу
			table.append(str(j)+"\t"+str(Acount)+"\t"+str(Ccount)+"\t"+str(Gcount)+"\t"+str(Tcount)+"\t"+str(PA)+"\t"+str(PC)+"\t"+str(PG)+"\t"+str(PT)+"\t"+str(CovJ)+"\t"+pos+"\t"+cov_arr[j][pos][3]+"\t"+cov_arr[j][pos][2]+"\t"+str(DivJ)+"\t"+str(ShanJ))
			transfac.append(str(j)+"\t"+str(Acount)+"\t"+str(Ccount)+"\t"+str(Gcount)+"\t"+str(Tcount))
	#Поиск SNP
			p_arr = [PA, PC, PG, PT]
			p_arr.sort(reverse=True)
			if data_type == "illumina":
				snp_sens = 0.05
			elif data_type == "nanopore":
				snp_sens = 0.1
			if p_arr[1] >= snp_sens and CovJ > 10:
				table[len(table)-1] += "\t"+str(PA)+"\t"+str(PC)+"\t"+str(PG)+"\t"+str(PT)
				if 'naive_SNPS' not in report_arr[sample_number][gene].keys():
					report_arr[sample_number][gene]['naive_SNPS'] = {}
				report_arr[sample_number][gene]['naive_SNPS'][j] = p_arr
			else:
				table[len(table)-1] +="\t\t\t\t"
			Shan.append(ShanJ)
			Cov.append(CovJ)
			Div.append(DivJ)
		print("\nCalculation finished\n")
		with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.transfac', 'w') as transfac_file:
			for el in transfac:
				transfac_file.write(el+"\n")
		os.system("weblogo --alphabet ACGT -F eps -U probability --resolution 600 -s large --scale-width no -n 250 --color-scheme classic -D transfac < "+sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+".transfac > "+sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_logo.eps") 		
	#Рисование графика энтропии Шеннона
		print("\nDrawing Shannon graph\n")
		mSh = round(statistics.mean(Shan), 6)
		vSh = statistics.stdev(Shan)
		v2Sh = vSh*2
		v3Sh = vSh*3
		v5Sh = vSh*5
		v10Sh = vSh*10
		mvector = []
		vvector = []
		v2vector = []
		v3vector = []
		v5vector = []
		v10vector = []
		qua_lim_arr = []
		qua_min_arr = []
		cov_lim_arr = []
		for i in range(0, len(ks)):
			mvector.append(mSh)
			vvector.append(mSh+vSh)
			v2vector.append(mSh+v2Sh)
			v3vector.append(mSh+v3Sh)
			v5vector.append(mSh+v5Sh)
			v10vector.append(mSh+v10Sh)
			qua_lim_arr.append(30)
			qua_min_arr.append(20)
			cov_lim_arr.append(100)
			if Shan[i] > mSh+v5Sh:
				if 'Shannon5_SNPS' not in report_arr[sample_number][gene].keys():
					report_arr[sample_number][gene]['Shannon5_SNPS'] = {}
				report_arr[sample_number][gene]['Shannon5_SNPS'][ks[i]] = Shan[i]
				table[i+1] += "\t"+str(Shan[i])
			else:
				table[i+1] += "\t"
			if Shan[i] > mSh+v10Sh:
				if 'Shannon10_SNPS' not in report_arr[sample_number][gene].keys():
					report_arr[sample_number][gene]['Shannon10_SNPS'] = {}
				report_arr[sample_number][gene]['Shannon10_SNPS'][ks[i]] = Shan[i]
				table[i+1] += "\t"+str(Shan[i])
			else:
				table[i+1] += "\t"
		ppl.plot(ks, Shan, 'b', ks, mvector, 'g--', ks, v3vector, 'r--', ks, v5vector, 'y--', ks, v10vector, 'k--')
		ppl.title(sample_name+" "+gene+",\nmean Shannon enthropy = "+str(mSh)+",\nplus 3(red), 5(yellow), 10(black) stdev")
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Shannon.svg')
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Shannon.png')
		ppl.clf()
		ppl.cla()
		print("\nShannon graph ready\n")
	#Рисование графика покрытия
		print("\nDrawing coverage graph\n")
		m = round(statistics.median(Cov), 2)
		ppl.plot(ks, Cov, 'b', ks, cov_lim_arr, 'r--', dp[sample_name+'_'+gene][0], dp[sample_name+'_'+gene][1], 'g')
		if len(ks) != len(dp[sample_name+'_'+gene][0]):
			print("\nWarning! Segment length has changed after removing duplicates!")
		ppl.title(sample_name+" "+gene+",\n coverage before (blue) and after (green) removing of PCR duplicates,\nmedian coverage before rmdup = "+str(m))
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Coverage.svg')
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Coverage.png')
		print("\nCoverage graph ready\n")
		ppl.clf()
		ppl.cla()
	#Рисование графика качества прочтения
		print("\nDrawing read quality graph\n")
		bqm = round(statistics.median(bq_arr), 2)
		ppl.plot(ks, bq_arr, 'b', ks, qua_lim_arr, 'g--', ks, qua_min_arr, 'y--')
		ppl.title(sample_name+" "+gene+",\nmedian base quality = "+str(bqm))
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_BaseQuality.svg')
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_BaseQuality.png')
		print("\nQuality graph ready\n")
		ppl.clf()
		ppl.cla()
	#Рисование графика качества выравнивания
		print("\nDrawing mapping quality graph\n")
		mqm = round(statistics.median(mq_arr), 2)
		ppl.plot(ks, mq_arr, 'b', ks, qua_lim_arr, 'r--')
		ppl.title(sample_name+" "+gene+",\nmedian mapping quality = "+str(mqm))
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_MapQuality.svg')
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_MapQuality.png')
		print("\nQuality graph ready\n")
		ppl.clf()
		ppl.cla()
	#Рисование графика неоднородности
		print("\nDrawing sequence diversity graph\n")
		sdm = round(statistics.mean(Div), 6)
		vd = statistics.stdev(Div)
		sdm_vector = []
		sdmv_vector = []
		sdmv2_vector = []
		sdmv3_vector = []
		for i in range(0, len(ks)):
			sdm_vector.append(sdm)
			sdmv_vector.append(sdm + vd)
			sdmv2_vector.append(sdm + vd*2)
			sdmv3_vector.append(sdm + vd *3)
		ppl.plot(ks, Div, 'b', ks, sdm_vector, 'g--', ks, sdmv_vector, 'r--', ks, sdmv2_vector, 'y--',ks, sdmv3_vector, 'k--')
		ppl.title(sample_name+" "+gene+",\nmean sequence diversity = "+str(sdm)+",\nplus 3(red), 5(yellow), 10(black) stdev")
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_SeqDiversity.svg')
		ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_SeqDiversity.png')
		print("\nSequence diversity ready\n")
		ppl.clf()
		ppl.cla()
		if "FB_SNPS" in report_arr[sample_number][gene].keys() and len(report_arr[sample_number][gene]["FB_SNPS"]) > 0:
			for sn in report_arr[sample_number][gene]["FB_SNPS"].keys():
				table[int(sn)] += "\t"+report_arr[sample_number][gene]["FB_SNPS"][sn]
		with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_problems.txt', 'r') as probl_file:
			probs = probl_file.readlines()
			prob_count = len(probs)
			pl = {}
			for st in probs:
				sts = st.split("\t")
				pl[sts[0]] = sts[1]
		with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.table', 'w') as table_file:
			for el in table:
				table_file.write(el+"\n")
		ans[gene] = {}
		ans[gene]['length'] = len(cons)
		ans[gene]['coverage'] = m
		ans[gene]['problems'] = prob_count
		ans[gene]['basequality'] = bqm
		ans[gene]['mapquality'] = mqm
		ans[gene]['seqdiversity'] = sdm
		report_arr[sample_number][gene]['length'] = len(cons)
		report_arr[sample_number][gene]['coverage'] = m
		report_arr[sample_number][gene]['problems'] = prob_count
		report_arr[sample_number][gene]['problems_list'] = pl
		report_arr[sample_number][gene]['basequality'] = bqm
		report_arr[sample_number][gene]['mapquality'] = mqm
		report_arr[sample_number][gene]['seqdiversity'] = sdm
		report_arr[sample_number][gene]['shannon'] = mSh
		report_arr[sample_number][gene]['length_rmdup'] = len(dp[sample_name+'_'+gene][0])
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return ans

def run_minimap2(ref_file, sample_number, step):
	#Выравнивание ридов с minimap2
	#Run minimap2 read alignment
	print("\nAlign reads by Minimap2 for sample "+sample_number+". Further output from Minimap2:\n")
	if data_type == 'illumina':
		minimap_cmd = minimap_path+'minimap2 -x sr -R "@RG\\tID:'+sample_number+'\\tSM:'+sample_number+'\\tPL:ILLUMINA\\tLB:library1\\t" -t '+str(os.cpu_count()-1)+' --end-bonus=1 -2 -a -o '+sample_number+os.sep+sample_number+'_'+step+'.sam '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq'
	elif data_type == 'nanopore':
		minimap_cmd = minimap_path+'minimap2 -x map-ont -R "@RG\\tID:'+sample_number+'\\tSM:'+sample_number+'\\tPL:NANOPORE\\tLB:library1\\t" -t '+str(os.cpu_count()-1)+' --end-bonus=1 -2 -a -o '+sample_number+os.sep+sample_number+'_'+step+'.sam '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq'
	os.system(minimap_cmd)
	print("\nMinimap2 alignment finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def run_centifuge(sample_number, mode, total, report_arr):
	#Определение сегментов для выравнивания
	print("\nUsing Centrifuge for selection of segments")
	centrifuge_cmd = centrifuge_path+'centrifuge-class.exe -x '+reference_path+mode+' -p '+str(os.cpu_count())+' -U '+sample_number+os.sep+sample_number+'.fastq --report-file '+sample_number+os.sep+sample_number+'_report.txt --min-hitlen 15 > '+os.devnull
	os.system(centrifuge_cmd)
	print("Centrifuge classification finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	with open(sample_number+os.sep+sample_number+'_report.txt', 'r') as rep:
		rep_arr = rep.readlines()
	seg_list = {'1':'root', '2':'Influenza_A', '3':'A_PB2', '4':'A_PB1', '5':'A_PA', '6':'A_H1', '7':'A_H2', '8':'A_H3', '9':'A_H4', '10':'A_H5', '11':'A_H6', '12':'A_H7', '13':'A_H8', '14':'A_H9', '15':'A_H10', '16':'A_H11', '17':'A_H12', '18':'A_H13', '19':'A_H14', '20':'A_H15', '21':'A_H16', '22':'A_H17', '23':'A_H18', '24':'A_NP', '25':'A_N1', '26':'A_N2', '27':'A_N3', '28':'A_N4', '29':'A_N5', '30':'A_N6', '31':'A_N7', '32':'A_N8', '33':'A_N9', '34':'A_N10', '35':'A_N11', '36':'A_M', '37':'A_NS', '38':'Influenza_B', '39':'B_PB2', '40':'B_PB1', '41':'B_PA', '42':'B_BHAvic', '43': 'B_BHAyam', '44': 'B_NP', '45':'B_BNAvic', '46':'B_BNAyam', '47':'B_M', '48':'B_NS'}
	surface = [['A_H1', 'A_H2', 'A_H3', 'A_H4', 'A_H5', 'A_H6', 'A_H7', 'A_H8', 'A_H9', 'A_H10', 'A_H11', 'A_H12', 'A_H13', 'A_H14', 'A_H15', 'A_H16', 'A_H17','A_H18', 'B_BHAvic', 'B_BHAyam'], ['A_N1', 'A_N2', 'A_N3', 'A_N4', 'A_N5', 'A_N6', 'A_N7', 'A_N8', 'A_N9', 'A_N10', 'A_N11', 'B_BNAvic', 'B_BNAyam']]
	surf_list = [{},{}]
	gene_list = []
	report_arr[sample_number]['centrifuge'] = {}
	if data_type == 'illumina':
		sens = 500
	elif data_type =='nanopore':
		sens = 25
	for el in rep_arr:
		els = el.split()
		if els[0] != "name":
			seg_name = seg_list[els[0]]
			report_arr[sample_number][seg_name] = {}
			report_arr[sample_number][seg_name]['mapped_reads'] = int(els[4])
			report_arr[sample_number]['centrifuge'][seg_name] = int(els[4])
			if int(els[4]) > sens:
				if seg_name in surface[0]:
					surf_list[0][seg_name] = int(els[4])
				elif seg_name in surface[1]:
					surf_list[1][seg_name] = int(els[4])
				else:
					gene_list.append(seg_name)
	if len(surf_list[0]) > 0:
		surf_list0_sorted = sorted(surf_list[0].items(), key=operator.itemgetter(1), reverse=True)
		gene_list.append(surf_list0_sorted[0][0])
	if len(surf_list[1]) > 0:
		surf_list1_sorted = sorted(surf_list[1].items(), key=operator.itemgetter(1), reverse=True)
		gene_list.append(surf_list1_sorted[0][0])
	refs = []
	for seq_record in SeqIO.parse(reference_path+mode+'.fasta', 'fasta'):
		refs.append(seq_record)
	with open(sample_number+os.sep+sample_number+'_start_ref.fasta', 'w') as f:
		for gene in gene_list:
			for ref in refs:
				if gene == ref.id:
					f.write('>'+ref.id+'\n'+str(ref.seq)+'\n')
	report_arr[sample_number]['gene_list'] = gene_list
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return gene_list

def call_variants(sample_number):
	#Поиск SNP c использованием FreeBayes
	print("Realign data for SNP calling with FreeBayes\n")
	if align_mode == "bwa-mem":
		run_bwa(sample_number+os.sep+sample_number+'.fasta', sample_number, 'v')
	elif align_mode == "minimap2":
		run_minimap2(sample_number+os.sep+sample_number+'.fasta', sample_number, 'v')
	samtools_cmd1 = samtools_path+'samtools faidx '+ sample_number+os.sep+sample_number+'.fasta'
	samtools_cmd2 = samtools_path+'samtools dict '+ sample_number+os.sep+sample_number+'.fasta > '+ sample_number+os.sep+sample_number+'.dict'
	if data_type == 'illumina':
		samtools_cmd3 = samtools_path+'samtools view -b -F 0x4 -t '+sample_number+os.sep+sample_number+'.fasta.fai '+sample_number+os.sep+sample_number+'_v.sam | '+samtools_path+'samtools sort -n -l 9 -@ '+str(os.cpu_count())+' | '+samtools_path+'samtools fixmate -rm - - | '+samtools_path+'samtools sort -l 9 -@ '+str(os.cpu_count())+' | '+samtools_path+'samtools calmd -rbA - '+sample_number+os.sep+sample_number+'.fasta > '+sample_number+os.sep+sample_number+'_var.bam'
	elif data_type == 'nanopore':
		samtools_cmd3 = samtools_path+'samtools view -b -F 0x4 -t '+sample_number+os.sep+sample_number+'.fasta.fai '+sample_number+os.sep+sample_number+'_v.sam | '+samtools_path+'samtools sort -l 9 -@ '+str(os.cpu_count())+' > '+sample_number+os.sep+sample_number+'_var.bam'
	samtools_cmd4 = samtools_path+'samtools index '+sample_number+os.sep+sample_number+'_var.bam'
	picard_cmd = 'java -Xmx512m -jar '+picard_path+'picard.jar MarkDuplicates I='+sample_number+os.sep+sample_number+'_var.bam'+' O='+sample_number+os.sep+sample_number+'_rmdup.bam'+' M='+sample_number+os.sep+sample_number+'_markdup.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true'
	samtools_cmd5 = samtools_path+'samtools index '+sample_number+os.sep+sample_number+'_rmdup.bam'
	samtools_cmd6 = samtools_path+'samtools depth '+sample_number+os.sep+sample_number+'_rmdup.bam > '+sample_number+os.sep+sample_number+'_rmdup_depth.txt'
	freebayes_cmd1 = freebayes_path+'freebayes -f '+sample_number+os.sep+sample_number+'.fasta -p 1 '+sample_number+os.sep+sample_number+'_rmdup.bam > '+sample_number+os.sep+sample_number+'_raw_freebayes.vcf'
	vcflib_cmd1 = freebayes_path+'vcfallelicprimitives -kg '+sample_number+os.sep+sample_number+'_raw_freebayes.vcf | '+freebayes_path+'vcffilter -f "QUAL > 20" > '+sample_number+os.sep+sample_number+'.vcf'
	os.system(samtools_cmd1)
	os.system(samtools_cmd2)
	os.system(samtools_cmd3)
	os.system(samtools_cmd4)
	print("Removing PCR duplicates with picard. Further output from Picard:\n")
	os.system(picard_cmd)
	os.system(samtools_cmd5)
	os.system(samtools_cmd6)
	wst = get_cons2(sample_number, sample_number+os.sep+sample_number+'.fasta', 'rmdup')
	print("Run SNP search using FreeBayes\n")
	os.system(freebayes_cmd1)
	os.system(vcflib_cmd1)
	os.remove(sample_number+os.sep+sample_number+'_v.sam')
	os.remove(sample_number+os.sep+sample_number+'_var.bam')
	os.remove(sample_number+os.sep+sample_number+'_var.bam.bai')
	with open(sample_number+os.sep+sample_number+'.vcf', 'r') as vcf_file:
		snps_arr = vcf_file.readlines()
	for snp_str in snps_arr:
		if snp_str[0] != "#":
			snp_str_sp = snp_str.split("\t")
			gn_arr = snp_str_sp[0].split("_")
			gn = gn_arr[len(gn_arr)-2]+"_"+gn_arr[len(gn_arr)-1]
			if "FB_SNPS" not in report_arr[sample_number][gn].keys():
				report_arr[sample_number][gn]["FB_SNPS"] = {}
				report_arr[sample_number][gn]["FB_SNPS"][snp_str_sp[1]] = snp_str_sp[3] + "->"+snp_str_sp[4]+"\t"+snp_str_sp[9]+" ("+snp_str_sp[8]+")"
			else:
				report_arr[sample_number][gn]["FB_SNPS"][snp_str_sp[1]] = snp_str_sp[3] + "->"+snp_str_sp[4]+"\t"+snp_str_sp[9]+" ("+snp_str_sp[8]+")"

def parse_idxstats(sample_number, step):
	#Разбираем выдачу Samtools idxstats
	#Parse Samtools idxstats output
	print("\nParsing Samtools idxstats output\n")
	surface = [['A_H1pdm', 'A_H1av', 'A_H1hs', 'A_H2', 'A_H3', 'A_H4', 'A_H5', 'A_H6', 'A_H7', 'A_H8', 'A_H9', 'A_H10', 'A_H11', 'A_H12', 'A_H13', 'A_H14', 'A_H15', 'A_H16', 'A_H17','A_H18', 'B_BHAvic', 'B_BHAyam'], ['A_N1', 'A_N2', 'A_N3', 'A_N4', 'A_N5', 'A_N6', 'A_N7', 'A_N8', 'A_N9', 'A_N10', 'A_N11', 'B_BNAvic', 'B_BNAyam']]
	surf_list = [{},{}]
	gene_list = []
	if data_type == 'illumina':
		sens = 700
	elif data_type =='nanopore':
		sens = 25
	with open(sample_number+os.sep+sample_number+'_'+step+'_idxstats.txt','r') as idx:
		for line in idx.readlines():
			line_split = line.split("\t")
			if line_split[0] != "*" and int(line_split[2]) != 0:
				tmp_arr = line_split[0].split('_')
				if len(tmp_arr) < 3:
					seg = line_split[0]
				else:
					seg = tmp_arr[len(tmp_arr)-2] +'_'+ tmp_arr[len(tmp_arr)-1]
				if seg not in report_arr[sample_number].keys():
					report_arr[sample_number][seg] = {}
				report_arr[sample_number][seg]['mapped_reads'] = int(line_split[2])
				if 'seg_stats' not in report_arr[sample_number].keys():
					report_arr[sample_number]['seg_stats'] = {}
				report_arr[sample_number]['seg_stats'][seg] = int(line_split[2])
				if int(line_split[2]) > sens:
					if seg in surface[0]:
						surf_list[0][seg] = int(line_split[2])
					elif seg in surface[1]:
						surf_list[1][seg] = int(line_split[2])
					else:
						gene_list.append(seg)
			elif line_split[0] == "*":
				report_arr[sample_number]['unmapped_reads'] = int(line_split[3])
	for arr in surf_list:
		if len(arr) > 0:
			arr_sorted = sorted(arr.items(), key=operator.itemgetter(1), reverse=True)
			gene_list.append(arr_sorted[0][0])
	refs = []
	if step == "s":
		print("\nPrepare reference set\n")
		for seq_record in SeqIO.parse(reference_path+mode+'.fasta', 'fasta'):
			refs.append(seq_record)
		with open(sample_number+os.sep+sample_number+'_start_ref.fasta', 'w') as f:
			ref_segs = ''
			for gene in gene_list:
				ref_segs = ref_segs + gene + ' '
				for ref in refs:
					if gene == ref.id:
						f.write('>'+ref.id+'\n'+str(ref.seq)+'\n')
		report_arr[sample_number]['gene_list'] = gene_list
		samtools_cmd = samtools_path+'samtools view -b  -o '+sample_number+os.sep+sample_number+'_'+step+'_s_tmp.bam '+sample_number+os.sep+sample_number+'_s.bam '+ref_segs
		os.system(samtools_cmd)
		src = sample_number+os.sep+sample_number+'_'+step+'_s_tmp.bam'
		dest = sample_number+os.sep+sample_number+'_s.bam'
		shutil.move(src, dest)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return gene_list
	
	

	
print("Start\n")
#Parse input data
descriptionText = "Illumina read mapping for Influenza A and B"
parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-list", dest="list", required=True, help="Path to sample list file")
parser.add_argument("-mode", dest="mode", required=False, help="Choose segment selection mode: all", default="all")
parser.add_argument("-align_mode", dest="align_mode", required=False, help="Choose alignment algorithm: bwa-mem, minimap2", default="bwa-mem")
parser.add_argument("-data_type", dest="data_type", required=False, help="Choose input data type: illumina, nanopore", default="illumina")
args = parser.parse_args()
mode = args.mode
#mode = 'all'
align_mode = args.align_mode
#align_mode = 'bwa-mem'
data_type = args.data_type
#Open sample list
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
with open(args.list, 'r') as sl:
	sample_list = sl.readlines()
if not os.path.isdir('result'):
	os.mkdir("result")

#Open report file
mapping_table = open('mapping_table.txt', 'w')
mapping_table.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'\n')
mapping_table.close()
report_arr = {}
report_arr['mode'] = mode
report_arr['align_mode'] = align_mode
report_arr['data_type'] = data_type
report_arr['list'] = args.list
report_arr['sample_list'] = []
for sample in sample_list:
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	sample_splitted = sample.split()
	sample_number = sample_splitted[0]
	sample_name = sample_splitted[1]
	plate_number = sample_splitted[2]
	if len(sample_splitted) == 4:
		index1 = sample_splitted[3]
		index2 = ''
	if len(sample_splitted) == 5:
		index1 = sample_splitted[3]
		index2 = sample_splitted[4]
	report_arr[sample_number] = {}
	report_arr[sample_number]['sample_number'] = sample_number
	report_arr[sample_number]['sample_name'] = sample_name
	report_arr[sample_number]['plate_number'] = plate_number
	report_arr[sample_number]['index1'] = index1
	report_arr[sample_number]['index2'] = index2
	report_arr[sample_number]['start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	copy_data(sample_number, plate_number, data_type)
	sample_stats = open(sample_number+os.sep+sample_number+"_stats.txt", "w")
	sample_stats.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	if not os.path.isfile(sample_number+os.sep+sample_number+".fastq") and data_type == 'illumina':
		fastqc_data(sample_number+os.sep+sample_number+'_R1.fastq.gz')
		fastqc_data(sample_number+os.sep+sample_number+'_R2.fastq.gz')
		trim_data(sample_number, trimmomatic_path)
		ttl = conc_data(sample_number)
		fastqc_data(sample_number+os.sep+sample_number+'.fastq')
	elif not os.path.isfile(sample_number+os.sep+sample_number+".fastq") and data_type == 'nanopore':
		fastqc_data(sample_number+os.sep+sample_number+'_unfiltered.fastq.gz')
		trim_data(sample_number, trimmomatic_path)
		fastqc_data(sample_number+os.sep+sample_number+'.fastq')
		with open(sample_number+os.sep+sample_number+'.fastq', 'r') as fq:
			fq_arr = fq.readlines()
			ttl = len(fq_arr)//4
	else:
		with open(sample_number+os.sep+sample_number+'.fastq', 'r') as fq:
			fq_arr = fq.readlines()
			ttl = len(fq_arr)//4
	report_arr[sample_number]['read_number'] = ttl	
	mapping_arr = {}
	print("First run for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	#gene_list = run_centifuge(sample_number, mode, ttl, report_arr)
	if not os.path.isfile(sample_number+os.sep+sample_number+'_ext_s.sam'):
		if align_mode == "bwa-mem":
			run_bwa(reference_path+mode+'.fasta', sample_number, 's')
		elif align_mode == "minimap2":
			run_minimap2(reference_path+mode+'.fasta', sample_number, 's')
	else:
		dest = sample_number+os.sep+sample_number+'_s.sam'
		src = sample_number+os.sep+sample_number+'_ext_s.sam'
		shutil.copyfile(src, dest)
	run_samtools2(sample_number, reference_path+mode+'.fasta', 's', report_arr)
	gene_list = parse_idxstats(sample_number, 's')
	for gene in gene_list:
		print(gene+"\n")
		if not os.path.isdir(sample_number+os.sep+gene):
			os.mkdir(sample_number+os.sep+gene)
	surf1 = gene_list[len(gene_list)-2].split("_")
	surf2 = gene_list[len(gene_list)-1].split("_")
	subtype = surf1[1]+surf2[1]
	report_arr[sample_number]['subtype'] = subtype
	if subtype in ["H3N2", "H1pdmN1", "BHAvicBNAvic", "BHAyamBNAyam"]:
		src = reference_path+subtype+'.fasta'
		dest = sample_number+os.sep+sample_number+'_start_ref.fasta'
		shutil.copyfile(src,dest)
	cnt0 = get_cons2(sample_number, sample_number+os.sep+sample_number+'_start_ref.fasta', 's')
	print("Subtype for sample: "+sample_number+" - "+subtype+"\n")
	if not os.path.isfile(sample_number+os.sep+sample_number+'_ext_1.sam'):
		if align_mode == "bwa-mem":
			run_bwa(sample_number+os.sep+sample_number+'_start_ref.fasta', sample_number, '1')
		elif align_mode == "minimap2":
			run_minimap2(sample_number+os.sep+sample_number+'_start_ref.fasta', sample_number, '1')
	else:
		dest = sample_number+os.sep+sample_number+'_1.sam'
		src = sample_number+os.sep+sample_number+'_ext_1.sam'
		shutil.copyfile(src, dest)
	run_samtools2(sample_number, sample_number+os.sep+sample_number+'_start_ref.fasta', '1', report_arr)
	gl = parse_idxstats(sample_number, '1')
	cnt = get_cons2(sample_number, sample_number+os.sep+sample_number+'_start_ref.fasta', '1')
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	for gene in gene_list:
		if len(cnt[gene])-1 != len(cnt[gene]["CON"]):
			print("Warning! Possible indels found in gene "+gene+"! ("+str(len(cnt[gene])-1)+" != "+str(len(cnt[gene]["CON"]))+")")
		print("Length of gene "+gene+" of sample "+sample_number+" after first run: "+str(len(cnt[gene])-1)+"\n")
	print("Second run for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	if align_mode == "bwa-mem":
		run_bwa(sample_number+os.sep+sample_number+'_1.fasta', sample_number, '')
	elif align_mode == "minimap2":
		run_minimap2(sample_number+os.sep+sample_number+'_1.fasta', sample_number, '')
	run_samtools2(sample_number, sample_number+os.sep+sample_number+'_1.fasta', '', report_arr)
	gl = parse_idxstats(sample_number, '')
	cnt1 = get_cons2(sample_number, sample_number+os.sep+sample_number+'_1.fasta', '')
	call_variants(sample_number)
	cover_ar = get_stats(cnt1, sample_number, sample_name, report_arr)
	for gene in gene_list: 
		mapping_arr[gene] = report_arr[sample_number][gene]['mapped_reads']
	#print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		copy_result(sample_number, gene, subtype)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")


	print("Write mapping data to report file for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	report_arr['sample_list'].append(sample_number)
	mapping_table = open('mapping_table.txt', 'a')
	mapping_table.write(sample_name+'_'+sample_number+'\n')
	sample_stats.write(sample_name+'_'+sample_number+'\n')
	mapping_table.write(subtype+'\n')
	sample_stats.write(subtype+'\n')
	mapping_table.write(str(ttl)+'\n')
	sample_stats.write(str(ttl)+'\n')
	mapping_arr_sorted = sorted(mapping_arr.items(), key=operator.itemgetter(1), reverse=True)
	for el in mapping_arr_sorted:
		mapping_table.write("\t"+el[0]+"\t"+str(el[1])+"\t")
		sample_stats.write("\t"+el[0]+"\t"+str(el[1])+"\t")
		if cover_ar.get(el[0]):
			tm = cover_ar.get(el[0])
			mapping_table.write(str(tm['length'])+"\t"+str(tm['coverage'])+"\t"+str(tm['problems'])+"\t"+str(tm['basequality'])+"\t"+str(tm['mapquality'])+"\t"+str(tm['seqdiversity'])+"\n")
			sample_stats.write(str(tm['length'])+"\t"+str(tm['coverage'])+"\t"+str(tm['problems'])+"\t"+str(tm['basequality'])+"\t"+str(tm['mapquality'])+"\t"+str(tm['seqdiversity'])+"\n")
		else:
			mapping_table.write("< than 1 percent of total reads or not selected for analysis\n")
			sample_stats.write("< than 1 percent of total reads or not selected for analysis\n")
	sample_stats.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	mapping_table.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	sample_stats.close()
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	report_arr[sample_number]['finish_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	mapping_table.close()
	html_report('report.html', report_arr, sample_number)
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
print("Finish")
