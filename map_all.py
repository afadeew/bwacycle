#BWAcycle tool
#Read mapping for influenza sequencing on Illumina
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


trimmomatic_path = '..'+os.sep+'stuff'+os.sep+'trimmomatic'+os.sep
fastqc_path = '..'+os.sep+'stuff'+os.sep+'FastQC'
bwa_path = '..'+os.sep+'stuff'+os.sep+'bwa'+os.sep
samtools_path = '..'+os.sep+'stuff'+os.sep
bam_readcount_path = '..'+os.sep+'stuff'+os.sep
reference_path = '..'+os.sep+'Reference'+os.sep


def test():
	print("Test")

def copy_data(sample_number, plate_number):
	#Создаем папку по лабораторному номеру образца, копируем риды в нее
	if not os.path.isdir(sample_number):
		os.mkdir(sample_number)
	dist1 = sample_number+os.sep+sample_number+'_R1.fastq.gz'
	src1 = sample_number+'_S'+plate_number+'_L001_R1_001.fastq.gz'
	dist2 = sample_number+os.sep+sample_number+'_R2.fastq.gz'
	src2 = sample_number+'_S'+plate_number+'_L001_R2_001.fastq.gz'
	shutil.copyfile(src1, dist1)
	shutil.copyfile(src2, dist2)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	print("Copy data to new folder for sample "+sample_number+"\n")

def trim_data(sample_number, trimmomatic_path):
	#Запускаем Trimmomatic для обрезки входных данных по качеству
	trimmomatic_files = sample_number+os.sep+sample_number+'_R1.fastq.gz '+sample_number+os.sep+sample_number+'_R2.fastq.gz '+sample_number+os.sep+sample_number+'_fp.fq '+sample_number+os.sep+sample_number+'_fu.fq '+sample_number+os.sep+sample_number+'_rp.fq '+sample_number+os.sep+sample_number+'_ru.fq'
	print("Run Trimmomatic for quality clipping of data  for sample "+sample_number+". Further output from Trimmomatic:\n")
	os.system('java -jar '+trimmomatic_path+'trimmomatic-0.38.jar PE -threads '+str(os.cpu_count())+' -phred33 '+trimmomatic_files+' ILLUMINACLIP:'+trimmomatic_path+'adapters'+os.sep+'NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 HEADCROP:15')
	print("Trimmomatic finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def conc_data(sample_number):
	#Собираем 4 файла выдачи Trimmomatic в один FASTQ, считаем количество строк в полученном файле
	print("Concatenate clipped data to one file for sample "+sample_number+"\n")
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
	total_str = len(reads)
	print("Total reads obtained: "+str(len(reads)//4))
	with open(sample_number+os.sep+sample_number+'.fastq', 'w') as f:
		for item in reads:
			f.write(item)
	print("Remove temporary files for sample "+sample_number+"\n")
	os.remove(sample_number+os.sep+sample_number+'_fp.fq')
	os.remove(sample_number+os.sep+sample_number+'_fu.fq')
	os.remove(sample_number+os.sep+sample_number+'_rp.fq')
	os.remove(sample_number+os.sep+sample_number+'_ru.fq')
	return total_str

def fastqc_data(sample_number, fastqc_path):
	#Запускаем FASTQC для оценки качества в обработанном FASTQ
	print("Run FastQC to obtain quality report for sample "+sample_number+". Further output from FastQC:\n")
	if os.sep == "/":
		fastqc_par = 'java -Xmx512m -classpath '+fastqc_path+':'+fastqc_path+'/sam-1.103.jar:'+fastqc_path+'/jbzip2-0.9.jar -Djava.awt.headless=true uk.ac.babraham.FastQC.FastQCApplication '+sample_number+os.sep+sample_number+'.fastq'
	else:
		fastqc_par = 'java -Xmx512m -classpath '+fastqc_path+';'+fastqc_path+'/sam-1.103.jar;'+fastqc_path+'/jbzip2-0.9.jar -Djava.awt.headless=true uk.ac.babraham.FastQC.FastQCApplication '+sample_number+os.sep+sample_number+'.fastq'
	os.system(fastqc_par)
	print("FastQC finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def index_bwa(ref_file):
	#Делаем индекс для работы BWA из ref_file
	print("Make index from reference file "+ref_file+". Further output from BWA index:\n")
	bwa_cmd = bwa_path+'bwa index '+ref_file
	os.system(bwa_cmd)
	print("BWA index finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def run_bwa(ref_file, sample_number, gene, step):
	#Проверяем существование индекса
	if not os.path.isfile(ref_file+'.bwt'):
		index_bwa(ref_file)
	#Запускаем выравнивание BWA
	print("Align reads by BWA mem for sample "+sample_number+" and gene "+gene+". Further output from BWA mem:\n")
	if step != '':
		bwa_cmd = bwa_path+'bwa mem -aM -R "@RG\tID:'+sample_number+'\tSM:'+gene+'\tLB:library1\t" -t '+str(os.cpu_count())+' -k 13 -B 1 '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam'
	else:
		bwa_cmd = bwa_path+'bwa mem -aM -R "@RG\tID:'+sample_number+'\tSM:'+gene+'\tLB:library1\t" -t '+str(os.cpu_count())+' '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam'
	os.system(bwa_cmd)
	print("BWA mem alignment finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def parse_flagstat(sample_number,gene, step):
	#Разбираем файл со статистикой выравнивания
	print("Parsing flagstat data for sample "+sample_number+" and gene "+gene+"\n")
	with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_flagstat.txt') as f:
		flags = f.readlines()
	total_str = flags[0].split()
	total = int(total_str[0])
	print("Total reads: "+str(total)+"\n")
	#mapped_str = flags[2].split() #samtools 0.1.19
	mapped_str = flags[4].split()
	mapped = []
	mapped.append(int(mapped_str[0]))
	mapped.append(total)
	print("Mapped reads: "+str(mapped[0])+" of "+str(mapped[1])+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return mapped
	
def copy_result(sample_number, gene, HA, NA):
	#Копирование результатов FASTA, SVG и PNG в папку Result
	print("Copy resulting files to Result folder\n")
	subtype = HA + NA
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
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Shannon.svg"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"svg"+os.sep+sample_number+"_"+gene+"_Shannon.svg"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Coverage.svg"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"svg"+os.sep+sample_number+"_"+gene+"_Coverage.svg"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_BaseQuality.svg"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"svg"+os.sep+sample_number+"_"+gene+"_BaseQuality.svg"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_MapQuality.svg"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"svg"+os.sep+sample_number+"_"+gene+"_MapQuality.svg"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Shannon.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_Shannon.png"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Coverage.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_Coverage.png"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_BaseQuality.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_BaseQuality.png"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_MapQuality.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_MapQuality.png"
	shutil.copyfile(src, dest)
	with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+sample_number+"_"+gene+".fasta", 'r') as inp:
		fa = inp.readlines()
	ft = []
	if os.path.isfile("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+"total.fasta"):
		with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta", 'r') as total:
			ft = total.readlines()
	ft+=fa
	with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta", 'w') as outp:
		for item in ft:
			outp.write(item)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")

def html_report(report_file, report_arr, sample_number):
	fastqc_zip = zipfile.ZipFile(sample_number+os.sep+sample_number+'_fastqc.zip', 'r')
	source = fastqc_zip.open(sample_number+'_fastqc/Images/per_base_quality.png')
	target = open(sample_number+os.sep+sample_number+'_per_base_quality.png', "wb")
	shutil.copyfileobj(source, target)
	source.close()
	target.close()
	source = fastqc_zip.open(sample_number+'_fastqc/Images/per_base_sequence_content.png')
	target = open(sample_number+os.sep+sample_number+'_per_base_sequence_content.png', "wb")
	shutil.copyfileobj(source, target)
	source.close()
	target.close()
	source = fastqc_zip.open(sample_number+'_fastqc/Images/sequence_length_distribution.png')
	target = open(sample_number+os.sep+sample_number+'_sequence_length_distribution.png', "wb")
	shutil.copyfileobj(source, target)
	source.close()
	target.close()
	bqfile = open(sample_number+os.sep+sample_number+'_per_base_quality.png', 'rb')
	bqfiledata = bqfile.read()
	bqdata = base64.b64encode(bqfiledata)
	bqstring = bqdata.decode()
	report_arr['read_quality'] = bqstring
	bqfile.close()
	scfile = open(sample_number+os.sep+sample_number+'_per_base_sequence_content.png', 'rb')
	scfiledata = scfile.read()
	scdata = base64.b64encode(scfiledata)
	scstring = scdata.decode()
	report_arr['sequence_content'] = scstring
	scfile.close()
	ldfile = open(sample_number+os.sep+sample_number+'_sequence_length_distribution.png', 'rb')
	ldfiledata = ldfile.read()
	lddata = base64.b64encode(ldfiledata)
	ldstring = lddata.decode()
	report_arr['length_distribution'] = ldstring
	ldfile.close()
	html_page = "<html>\n<head>\n</head>\n<body>\n<p>\n<b>Sequence analysis report</b>\n</p>\n<br/>\n<table>\n<tr>\n"
	for i in range(1,13):
		html_page += "<td>\n<b>"+str(i)+"</b>\n</td>\n"
	html_page += "</tr>\n"
	html_page += "</table>\n"
	html_page += "</body>\n</html>"
	report_file.write(html_page)

def run_samtools2(sample_number,gene,step):
	#Преобразуем полученное от BWA выравнивание из SAM в BAM, сортируем и индексируем BAM-файл, получаем статистику по выравниванию, удаляем из BAM невыровненные риды
	print("Convert data to BAM-file, obtain flagstat data by Samtools for sample "+sample_number+" and gene "+gene+". Further output from Samtools:\n")
	samtools_cmd1 = samtools_path+'samtools view -bSu@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam' 
	samtools_cmd2 = samtools_path+'samtools sort -l 0 -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam -o '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	samtools_cmd3 = samtools_path+'samtools index -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	samtools_cmd4 = samtools_path+'samtools flagstat -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_flagstat.txt'
	os.system(samtools_cmd1)
	os.system(samtools_cmd2)
	os.system(samtools_cmd3)   #-F 0x4
	os.system(samtools_cmd4)
	print("Samtools finished\n")
	print("Remove temporary files for sample "+sample_number+" and gene "+gene+"\n")
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam')
	print("Compress BAM-files by Samtools for sample "+sample_number+" and gene "+gene+". Further output from Samtools:\n")
	samtools_cmd5 = samtools_path+'samtools view -bu -F 0x4 -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam'
	samtools_cmd6 = samtools_path+'samtools sort -@ '+str(os.cpu_count())+' -l 9 '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam -o '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	samtools_cmd7 = samtools_path+'samtools index -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	os.system(samtools_cmd5)
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam')
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam.bai')
	os.system(samtools_cmd6)
	os.system(samtools_cmd7)
	print("Samtools finished\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	print("Remove temporary files for sample "+sample_number+" and gene "+gene+"\n")
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam')
	
def get_cons2(sample_number, gene, ref_file, step):
	#Получаем консенсусную последовательность и статистику по позициям
	print("Extracting base counts from BAM-file with Bam-readcount for sample "+sample_number+" and gene "+gene+". Further output from Bam-readcount:\n")
	bamrc_cmd = bam_readcount_path+'bam-readcount -w 0 -f '+ref_file+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_readcount.txt'
	os.system(bamrc_cmd)
	print("Bam-readcount finished\n")
	with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_readcount.txt', 'r') as f:
		rc_str_arr = f.readlines()
	rc = {}
	for el in rc_str_arr:
		rc_arr = el.split()
		rc[int(rc_arr[1])] = {"A":rc_arr[5].split(':'), "C":rc_arr[6].split(':'), "G":rc_arr[7].split(':'), "T":rc_arr[8].split(':'), "INDEL":[], "COV":int(rc_arr[3]), "REF":rc_arr[2]}
		if len(rc_arr) > 10:
			ind_arr = []
			for elt in rc_arr[10:]:
				ind_arr.append(elt.split(':'))
			rc[int(rc_arr[1])]["INDEL"] = ind_arr
	print("Calculating consensus for sample "+sample_number+" and gene "+gene+"\n")
	con = ""
	dels = [0,0]
	probl_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_problems'+step+'.txt', 'w')
	ks = list(rc.keys())
	ks.sort()
	for i in ks:
		if rc[i]["COV"] < 10:
			print("Coverage is lower than 10 in position "+str(i+1)+". Please review basecount data.\n")
			probl_file.write("Coverage is lower than 10 in position "+str(i+1)+". Please review basecount data.\n")
		if dels[0] > 0 and dels[1] > rc[i]["COV"]:
			dels[0]-=1
		else:
			pos = rc[i]["REF"]
			if int(rc[i]["A"][1]) > int(rc[i]["C"][1]) and int(rc[i]["A"][1]) > int(rc[i]["G"][1]) and int(rc[i]["A"][1]) > int(rc[i]["T"][1]):
				pos = "A"
			if int(rc[i]["C"][1]) > int(rc[i]["A"][1]) and int(rc[i]["C"][1]) > int(rc[i]["G"][1]) and int(rc[i]["C"][1]) > int(rc[i]["T"][1]):
				pos = "C"
			if int(rc[i]["G"][1]) > int(rc[i]["C"][1]) and int(rc[i]["G"][1]) > int(rc[i]["A"][1]) and int(rc[i]["G"][1]) > int(rc[i]["T"][1]):
				pos = "G"
			if int(rc[i]["T"][1]) > int(rc[i]["C"][1]) and int(rc[i]["T"][1]) > int(rc[i]["G"][1]) and int(rc[i]["T"][1]) > int(rc[i]["A"][1]):
				pos = "T"
			if len(rc[i]["INDEL"]) > 0:
				indel_cov = [0, int(rc[i]["INDEL"][0][1])]
				for j in range(len(rc[i]["INDEL"])):
					if int(rc[i]["INDEL"][j][1]) > indel_cov[1]:
						indel_cov = [j, int(rc[i]["INDEL"][j][1])]
				if rc[i]["INDEL"][indel_cov[0]][0][0] == "+" and float(indel_cov[1])/float(rc[i]["COV"]) > 0.5:
					pos += rc[i]["INDEL"][indel_cov[0]][0][1:]
					print("Found more than 50 percent of insertions in position "+str(i)+". Please review basecount data.\n")
					probl_file.write("Found more than 50 percent of insertions in position "+str(i)+". Please review basecount data.\n")
				if rc[i]["INDEL"][indel_cov[0]][0][0] == "-" and float(indel_cov[1])/float(rc[i]["COV"]) > 0.5:
					dels = [len(rc[i]["INDEL"][indel_cov[0]][0][1:]), int(rc[i]["INDEL"][indel_cov[0]][1])]
					print("Found more than 50 percent of deletions in position "+str(i)+". Please review basecount data.\n")
					probl_file.write("Found more than 50 percent of deletions in position "+str(i)+". Please review basecount data.\n")
			con += pos
	if step != '':
		out_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.fasta', 'w')
		out_file.write('>'+sample_name+'_'+gene+'_intermediate_cons'+'\n'+con+'\n')
	else:
		out_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.fasta', 'w')
		out_file.write('>'+sample_name+'_'+gene+'\n'+con+'\n')
	out_file.close()
	probl_file.close()
	rc["CON"] = con
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return rc

def get_stats(cov_arr, sample_number, gene, sample_name):
	#Получение консенсусной последовательности, таблиц и графиков из массива частот встречаемости
	print("Calculating statistics for sample "+sample_number+" and gene "+gene+"\n")
	t = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.table', 'w')
	t.write("Position\tA\tC\tG\tT\tpA\tpC\tpG\tpT\tShannon\tCoverage\tConsensus\tA_SNP\tC_SNP\tG_SNP\tT_SNP\n")
	cons = cov_arr.pop("CON")
	Shan = []
	Cov = []
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
#Вывод результатов в таблицу
		t.write(str(j)+"\t"+str(Acount)+"\t"+str(Ccount)+"\t"+str(Gcount)+"\t"+str(Tcount)+"\t"+str(PA)+"\t"+str(PC)+"\t"+str(PG)+"\t"+str(PT)+"\t"+str(ShanJ)+"\t"+str(CovJ)+"\t"+pos)
#Поиск SNP
		if (PA > 0.05 and PA < 0.95):
			t.write("\t"+str(PA))
		else:
			t.write("\t")
		if (PC > 0.05 and PC < 0.95):
			t.write("\t"+str(PC))
		else:
			t.write("\t")
		if (PG > 0.05 and PG < 0.95):
			t.write("\t"+str(PG))
		else:
			t.write("\t")
		if (PT > 0.05 and PT < 0.95):
			t.write("\t"+str(PT))
		else:
			t.write("\t")
		t.write("\n")
		Shan.append(ShanJ)
		Cov.append(CovJ)
	t.close()
	print("Calculation finished\n")
#Рисование графика энтропии Шеннона
	print("Drawing Shannon graph\n")
	mSh = statistics.mean(Shan)
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
	cov_lim_arr = []
	for i in range(0, len(ks)):
		mvector.append(mSh)
		vvector.append(mSh+vSh)
		v2vector.append(mSh+v2Sh)
		v3vector.append(mSh+v3Sh)
		v5vector.append(mSh+v5Sh)
		v10vector.append(mSh+v10Sh)
		qua_lim_arr.append(30)
		cov_lim_arr.append(100)
	ppl.plot(ks, Shan, 'b', ks, mvector, 'g--', ks, v3vector, 'r--', ks, v5vector, 'y--', ks, v10vector, 'k--')
	ppl.title(sample_name+" "+gene+" Shannon enthropy, plus 3, 5, 10 stdev")
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Shannon.svg')
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Shannon.png')
	ppl.clf()
	ppl.cla()
	print("Shannon graph ready\n")
#Рисование графика покрытия
	print("Drawing coverage graph\n")
	m = round(statistics.mean(Cov), 2)
	ppl.plot(ks, Cov, 'b', ks, cov_lim_arr, 'r--')
	ppl.title(sample_name+" "+gene+", mean coverage = "+str(m))
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Coverage.svg')
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Coverage.png')
	print("Coverage graph ready\n")
	ppl.clf()
	ppl.cla()
#Рисование графика качества прочтения
	print("Drawing read quality graph\n")
	bqm = round(statistics.mean(bq_arr), 2)
	ppl.plot(ks, bq_arr, 'b', ks, qua_lim_arr, 'r--')
	ppl.title(sample_name+" "+gene+", mean base quality = "+str(bqm))
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_BaseQuality.svg')
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_BaseQuality.png')
	print("Quality graph ready\n")
	ppl.clf()
	ppl.cla()
#Рисование графика качества выравнивания
	print("Drawing mapping quality graph\n")
	mqm = round(statistics.mean(mq_arr), 2)
	ppl.plot(ks, mq_arr, 'b', ks, qua_lim_arr, 'r--')
	ppl.title(sample_name+" "+gene+", mean mapping quality = "+str(mqm))
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_MapQuality.svg')
	ppl.savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_MapQuality.png')
	print("Quality graph ready\n")
	ppl.clf()
	ppl.cla()
	with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_problems.txt', 'r') as probl_file:
		probs = probl_file.readlines()
		prob_count = len(probs)
	ans = []
	ans.append(len(cons))
	ans.append(m)
	ans.append(prob_count)
	ans.append(bqm)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	return ans



print("Start\n")
#Parse input data
descriptionText = "Illumina read mapping for influenza"
parser = argparse.ArgumentParser(description = descriptionText,formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-list", dest="list", required="true", help="Path to sample list file")
parser.add_argument("-mode", dest="mode", required="true", help="Path to alignment file")
args = parser.parse_args()
#Open sample list
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
with open(args.list, 'r') as sl:
	sample_list = sl.readlines()
if not os.path.isdir('result'):
	os.mkdir("result")
#Choose list of genes
if args.mode == "A":
	gene_list = ['H1','H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'PB2', 'PB1', 'PA', 'NP', 'M', 'NS']
elif args.mode == "A_seasonal":
	gene_list = ['H1', 'H3', 'N1', 'N2', 'PB2', 'PB1', 'PA', 'NP', 'M', 'NS']
elif args.mode == "H1":
	gene_list = ['H1', 'N1', 'PB2', 'PB1', 'PA', 'NP', 'M', 'NS']
elif args.mode == "H3":
	gene_list = ['H3', 'N2', 'PB2', 'PB1', 'PA', 'NP', 'M', 'NS']
elif args.mode == "B":
	gene_list = ['BHAvic', 'BNAvic', 'BHAyam', 'BNAyam', 'BPB2', 'BPB1', 'BPA', 'BNP', 'BM', 'BNS']
elif args.mode == "test":
	gene_list = ['H1', 'H3', 'N1', 'N2', 'BHAvic', 'BNAvic', 'BHAyam', 'BNAyam']
else:
	gene_list = ['H1','H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'PB2', 'PB1', 'PA', 'NP', 'M', 'NS', 'BHAvic', 'BNAvic', 'BHAyam', 'BNAyam', 'BPB2', 'BPB1', 'BPA', 'BNP', 'BM', 'BNS']
#Open report file
mapping_table = open('mapping_table.txt', 'w')
mapping_table.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'\n')
mapping_table.close()
report_arr = {}
for sample in sample_list:
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	sample_splitted = sample.split()
	sample_number = sample_splitted[0]
	sample_name = sample_splitted[1]
	plate_number = sample_splitted[2]
	report_arr[sample_number] = {}
	report_arr[sample_number]['number'] = sample_number
	report_arr[sample_number]['name'] = sample_name
	report_arr[sample_number]['plate'] = plate_number
	report_arr[sample_number]['start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	copy_data(sample_number, plate_number)
	sample_stats = open(sample_number+os.sep+sample_number+"_stats.txt", "w")
	sample_stats.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	trim_data(sample_number, trimmomatic_path)
	ttl = conc_data(sample_number)
	report_arr[sample_number]['read_number'] = ttl
	fastqc_data(sample_number, fastqc_path)	
	second_run = []
	mapping_arr = {}
	print("First run for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	for gene in gene_list:
		print(gene+"\n")
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		if not os.path.isdir(sample_number+os.sep+gene):
			os.mkdir(sample_number+os.sep+gene)
		run_bwa(reference_path+gene+'.fasta', sample_number, gene, '1')
		run_samtools2(sample_number, gene, '1')
		mapped = parse_flagstat(sample_number, gene, '1')
		mapping_arr[gene] = mapped[0]
		if mapped[0]/mapped[1] > 0.01 and gene not in ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'BHAyam', 'BHAvic', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'BNAyam', 'BNAvic']:
			second_run.append(gene)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	if args.mode == "A":
		HA_list = ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16']
		NA_list = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9']
		HAM = 'H1'
		NAM = 'N1'
	elif args.mode == "A_seasonal":
		HA_list = ['H1', 'H3']
		NA_list = ['N1', 'N2']
		HAM = 'H1'
		NAM = 'N1'
	elif args.mode == "H1":
		HA_list = ['H1']
		NA_list = ['N1']
		HAM = 'H1'
		NAM = 'N1'
	elif args.mode == "H3":
		HA_list = ['H3']
		NA_list = ['N2']
		HAM = 'H3'
		NAM = 'N2'
	elif args.mode == "B":
		HA_list = ['BHAyam', 'BHAvic']
		NA_list = ['BNAyam', 'BNAvic']
		HAM = 'BHAyam'
		NAM = 'BNAyam'
	elif args.mode == "test":
		HA_list = ['H1', 'H3', 'BHAyam', 'BHAvic']
		NA_list = ['N1', 'N2', 'BNAyam', 'BNAvic']
		HAM = 'H1'
		NAM = 'N1'
	else:
		HA_list = ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'BHAyam', 'BHAvic']
		NA_list = ['N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'BNAyam', 'BNAvic']
		HAM = 'H1'
		NAM = 'N1'
	print(str(mapping_arr)+"\n") #!!!!!!!!!!!!!!!!!!!!!!!!!!!
	max_HA = mapping_arr[HAM]
	for HA in HA_list:
		if mapping_arr[HA] > max_HA:
			HAM = HA
			max_HA = mapping_arr[HAM]
	max_NA = mapping_arr[NAM]
	for NA in NA_list:
		if mapping_arr[NA] > max_NA:
			NAM = NA
			max_NA = mapping_arr[NAM]
	subtype = HAM + NAM
	print("Subtype for sample: "+sample_number+" - "+subtype+"\n")
	cover_ar = {}
	print("Second run for sample "+sample_number+"\n")
	second_run.append(HAM)
	second_run.append(NAM)
	for gene in second_run:
		print(gene+"\n")
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		cnt = get_cons2(sample_number, gene, reference_path+gene+'.fasta', '1')
		print("Length of gene "+gene+" of sample "+sample_number+" after first run: "+str(len(cnt)-1)+"\n")
		run_bwa(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_1.fasta', sample_number, gene, '')
		run_samtools2(sample_number, gene, '')
		mapped1 = parse_flagstat(sample_number, gene, '')
		mapping_arr[gene] = mapped1[0]
		cnt1 = get_cons2(sample_number, gene, sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_1.fasta', '')
		print("Length of gene "+gene+" of sample "+sample_number+" after second run: "+str(len(cnt1)-1)+"\n")
		cv = get_stats(cnt1, sample_number, gene, sample_name)
		cover_ar[gene] = cv
		#print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		copy_result(sample_number, gene, HAM, NAM)
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")


	print("Write mapping data to report file for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	mapping_table = open('mapping_table.txt', 'a')
	mapping_table.write(sample_name+'_'+sample_number+'\n')
	sample_stats.write(sample_name+'_'+sample_number+'\n')
	mapping_table.write(subtype+'\n')
	sample_stats.write(subtype+'\n')
	mapping_table.write(str(mapped[1])+'\n')
	sample_stats.write(str(mapped[1])+'\n')
	mapping_arr_sorted = sorted(mapping_arr.items(), key=operator.itemgetter(1), reverse=True)
	for el in mapping_arr_sorted:
		mapping_table.write("\t"+el[0]+"\t"+str(el[1])+"\t")
		sample_stats.write("\t"+el[0]+"\t"+str(el[1])+"\t")
		if cover_ar.get(el[0]):
			tm = cover_ar.get(el[0])
			mapping_table.write(str(tm[0])+"\t"+str(tm[1])+"\t"+str(tm[2])+"\t"+str(tm[3])+"\n")
			sample_stats.write(str(tm[0])+"\t"+str(tm[1])+"\t"+str(tm[2])+"\t"+str(tm[3])+"\n")
		else:
			mapping_table.write("< than 1 percent of total reads or not selected for analysis\n")
			sample_stats.write("< than 1 percent of total reads or not selected for analysis\n")
	sample_stats.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	mapping_table.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	sample_stats.close()
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	report_arr[sample_number]['finish_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	mapping_table.close()
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
print("Finish")
