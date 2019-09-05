#BWAcycle tool
#Read mapping for influenza sequencing on Illumina
import sys
import os
import math
import shutil
from pylab import *
import statistics
import operator
from datetime import datetime
import base64
import argparse
import zipfile

trimmomatic_path = '..'+os.sep+'stuff'+os.sep+'trimmomatic'+os.sep
fastqc_path = '..'+os.sep+'stuff'+os.sep+'FastQC'
bwa_path = '..'+os.sep+'stuff'+os.sep+'bwa'+os.sep
samtools_path = '..'+os.sep+'stuff'+os.sep+'samtools'+os.sep
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
	print("Copy data to new folder for sample "+sample_number+"\n")

def trim_data(sample_number, trimmomatic_path):
	#Запускаем Trimmomatic для обрезки входных данных по качеству
	trimmomatic_files = sample_number+os.sep+sample_number+'_R1.fastq.gz '+sample_number+os.sep+sample_number+'_R2.fastq.gz '+sample_number+os.sep+sample_number+'_fp.fq '+sample_number+os.sep+sample_number+'_fu.fq '+sample_number+os.sep+sample_number+'_rp.fq '+sample_number+os.sep+sample_number+'_ru.fq'
	print("Run Trimmomatic for quality clipping of data  for sample "+sample_number+". Further output from Trimmomatic:\n")
	os.system('java -jar '+trimmomatic_path+'trimmomatic-0.38.jar PE -threads '+str(os.cpu_count())+' -phred33 '+trimmomatic_files+' ILLUMINACLIP:'+trimmomatic_path+'adapters'+os.sep+'NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30 HEADCROP:15')
	print("Trimmomatic finished\n")

def conc_data(sample_number):
	#Собираем 4 файла выдачи Trimmomatic в один FASTQ, считаем количество строк в полученном файле
	print("Concatenate clipped data to one file for sample "+sample_number+"\n")
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
	print(str(len(reads)//4))
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
	if not os.path.isfile('/usr/bin/fastqc'):
		fastqc_par = 'java -Xmx512m -classpath '+fastqc_path+';'+fastqc_path+'/sam-1.103.jar;'+fastqc_path+'/jbzip2-0.9.jar -Djava.awt.headless=true uk.ac.babraham.FastQC.FastQCApplication '+sample_number+os.sep+sample_number+'.fastq'
	else:
		fastqc_par = 'fastqc '+sample_number+os.sep+sample_number+'.fastq'
	os.system(fastqc_par)
	print("FastQC finished\n")

def index_bwa(ref_file):
	#Делаем индекс для работы BWA из ref_file
	print("Make index from reference file "+ref_file+". Further output from BWA index:\n")
	bwa_cmd = bwa_path+'bwa index '+ref_file
	os.system(bwa_cmd)
	print("BWA index finished\n")

def run_bwa(ref_file, sample_number, gene, step):
	#Запускаем выравнивание BWA
	print("Align reads by BWA mem  for sample "+sample_number+" and gene "+gene+". Further output from BWA mem:\n")
	if step != '':
		bwa_cmd = bwa_path+'bwa mem -t '+str(os.cpu_count())+' -k 13 -B 1 '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam'
	else:
		bwa_cmd = bwa_path+'bwa mem -t '+str(os.cpu_count())+' '+ref_file+' '+sample_number+os.sep+sample_number+'.fastq > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam'
	os.system(bwa_cmd)
	print("BWA mem alignment finished\n")

def run_samtools(sample_number,gene,step):
	#Преобразуем полученное от BWA выравнивание из SAM в BAM, сортируем и индексируем BAM-файл, получаем статистику по выравниванию, удаляем из BAM невыровненные риды
	print("Convert data to BAM-file, obtain flagstat data by Samtools for sample "+sample_number+" and gene "+gene+". Further output from Samtools:\n")
	samtools_cmd1 = samtools_path+'samtools view -b -S -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam'
	samtools_cmd2 = samtools_path+'samtools sort -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step
	samtools_cmd3 = samtools_path+'samtools index '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	samtools_cmd4 = samtools_path+'samtools flagstat '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_flagstat.txt'
	os.system(samtools_cmd1)
	os.system(samtools_cmd2)
	os.system(samtools_cmd3)   #-F 0x4
	os.system(samtools_cmd4)
	print("Samtools finished\n")
	print("Remove temporary files for sample "+sample_number+" and gene "+gene+"\n")
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam')
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam')
	src = sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	dest = sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_uncompressed.bam'
	shutil.copyfile(src, dest)
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam')
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam.bai')
	print("Compress BAM-files by Samtools for sample "+sample_number+" and gene "+gene+". Further output from Samtools:\n")
	samtools_cmd5 = samtools_path+'samtools view -b -F 0x4 -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_uncompressed.bam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_compressed.bam'
	samtools_cmd6 = samtools_path+'samtools sort -@ '+str(os.cpu_count())+' -l 9 '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_compressed.bam '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step
	samtools_cmd7 = samtools_path+'samtools index '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	os.system(samtools_cmd5)
	os.system(samtools_cmd6)
	os.system(samtools_cmd7)
	print("Samtools finished\n")
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_uncompressed.bam')
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_compressed.bam')
	print("Remove temporary files for sample "+sample_number+" and gene "+gene+"\n")

def parse_flagstat(sample_number,gene, step):
	#Разбираем файл со статистикой выравнивания
	print("Parsing flagstat data for sample "+sample_number+" and gene "+gene+"\n")
	with open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_flagstat.txt') as f:
		flags = f.readlines()
	total_str = flags[0].split()
	total = int(total_str[0])
	print("Total reads: "+str(total)+"\n")
	mapped_str = flags[2].split()
	mapped = []
	mapped.append(int(mapped_str[0]))
	mapped.append(total)
	print("Mapped reads: "+str(mapped[0])+" of "+str(mapped[1])+"\n")
	return mapped
	
def make_pileup(sample_number, gene, ref_file, step):
	#Получаем PILEUP-файл, содержащий статистику по позициям
	print("Making pileup-file by Samtools for sample "+sample_number+" and gene "+gene+". Further output from Samtools:\n")
	samtools_cmd = samtools_path+'samtools mpileup -f'+ref_file+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.txt'
	os.system(samtools_cmd)
	print("Samtools finished\n")

def get_cons(cov_arr, ref_file, sample_number, gene, sample_name, step):
	#Получение консенсусной последовательности из массива частот встречаемости
	print("Calculating consensus for sample "+sample_number+" and gene "+gene+"\n")
	with open(ref_file, 'r') as r:
		refs = r.readlines()
	ref = refs[1]
	sln = len(cov_arr)
	cons = ""
	offset = 0
	probl_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_problems'+step+'.txt', 'w')
	for i in range(len(cov_arr)):
		Acount = cov_arr[i][0]
		Gcount = cov_arr[i][2]
		Ccount = cov_arr[i][1]
		Tcount = cov_arr[i][3]
		gapcount = cov_arr[i][4]
		inscount = cov_arr[i][5]
		cov = cov_arr[i][6]
		if cov < 10:
			print("Coverage is lower than 10 in position "+str(i+1)+". Please review pileup data.\n")
			probl_file.write("Coverage is lower than 10 in position "+str(i+1)+". Please review pileup data.\n")
		if gapcount/cov > 0.1:
			print("Found more than 10 percent of deletions in position "+str(i+1)+". Please review pileup data.\n")
			probl_file.write("Found more than 10 percent of deletions in position "+str(i+1)+". Please review pileup data.\n")
		if inscount/cov > 0.1:
			print("Found more than 10 percent of insertions in position "+str(i+1)+". Please review pileup data.\n")
			probl_file.write("Found more than 10 percent of insertions in position "+str(i+1)+". Please review pileup data.\n")
		if step != '':
			if inscount != cov and i+offset < len(ref):
				con = ref[i+offset]
			elif inscount == cov and i+offset < len(ref):
				con = "N"
				offset = offset + 1
			else:
				con = "N"
			if gapcount == cov:
				offset = offset - 1
				con = ""
		else:
			con = "N"
		if Acount > Ccount and Acount > Gcount and Acount > Tcount:
			con = "A"
		if Ccount > Acount and Ccount > Gcount and Ccount > Tcount:
			con = "C"
		if Gcount > Ccount and Gcount > Acount and Gcount > Tcount:
			con = "G"
		if Tcount > Ccount and Tcount > Gcount and Tcount > Acount:
			con = "T"
		cons = cons + con
	probl_file.close()
	if step != '':
		out_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.fasta', 'w')
		out_file.write('>'+sample_name+'_'+gene+'_intermediate_cons'+'\n'+cons+'\n')
	else:
		out_file = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.fasta', 'w')
		out_file.write('>'+sample_name+'_'+gene+'\n'+cons+'\n')
	out_file.close()

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
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Quality.svg"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"svg"+os.sep+sample_number+"_"+gene+"_Quality.svg"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Shannon.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_Shannon.png"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Coverage.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_Coverage.png"
	shutil.copyfile(src, dest)
	src = sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_Quality.png"
	dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"png"+os.sep+sample_number+"_"+gene+"_Quality.png"
	shutil.copyfile(src, dest)
	if not os.path.isfile("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+"total.fasta"):
		src = "result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+sample_number+"_"+gene+".fasta"
		dest = "result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta"
		shutil.copyfile(src, dest)
	else:
		with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+sample_number+"_"+gene+".fasta", 'r') as inp:
			fa = inp.readlines()
		with open("result"+os.sep+subtype+os.sep+gene+os.sep+"fasta"+os.sep+"total.fasta", 'a') as outp:
			for item in fa:
				outp.write(item)

def shannon_consensus(cov_arr, sample_number, gene, sample_name):
	#Получение консенсусной последовательности, таблиц и графиков из массива частот встречаемости
	print("Calculating final consensus for sample "+sample_number+" and gene "+gene+"\n")
	t = open(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'.table', 'w')
	t.write("Position\tA\tC\tG\tT\tpA\tpC\tpG\tpT\tShannon\tCoverage\tConsensus\tA_SNP\tC_SNP\tG_SNP\tT_SNP\n")
	cons = ''
	m = 0
	Shan = []
	Cov = []
	qua_arr = []
	A_arr = []
	C_arr = []
	G_arr = []
	T_arr = []
	PA_arr = []
	PC_arr = []
	PG_arr = []
	PT_arr = []
	gap_arr = []
	ins_arr = []
	sln = len(cov_arr)
	prob_count = 0
	for j in range(sln):
		Acount = cov_arr[j][0]
		Gcount = cov_arr[j][2]
		Ccount = cov_arr[j][1]
		Tcount = cov_arr[j][3]
		gapcount = cov_arr[j][4]
		inscount = cov_arr[j][5]
		CovJ = cov_arr[j][6]
		A_arr.append(Acount)
		C_arr.append(Ccount)
		G_arr.append(Gcount)
		T_arr.append(Tcount)
		gap_arr.append(gapcount)
		ins_arr.append(inscount)
		if len(cov_arr[j]) > 7:
			qua = cov_arr[j][7]
		else:
			qua = 0
		qua_arr.append(qua)
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
		PA_arr.append(PA*100)
		PC_arr.append(PC*100)
		PG_arr.append(PG*100)
		PT_arr.append(PT*100)
		AENT = 0.0
		GENT = 0.0
		CENT = 0.0
		TENT = 0.0
		if PA != 0.0:
			AENT = -1*PA*log2(PA)
		if PG != 0.0:
			GENT = -1*PG*log2(PG)
		if PC != 0.0:
			CENT = -1*PC*log2(PC)
		if PT != 0.0:
			TENT = -1*PT*log2(PT)
		ShanJ = AENT + GENT + CENT + TENT
#Определение консенсуса
		Con = "N"
		if ((PA > PG) and (PA > PC) and (PA > PT)):
			Con = "A"
		if ((PG > PA) and (PG > PC) and (PG > PT)):
			Con = "G"
		if ((PC > PA) and (PC > PG) and (PC > PT)):
			Con = "C"
		if ((PT > PA) and (PT > PC) and (PT > PG)):
			Con = "T"
		if ((PT == 0) and (PA == 0) and (PC == 0) and (PG == 0)):
			Con = "-"
		cons = cons + Con
		t.write(str(j+1)+"\t"+str(Acount)+"\t"+str(Ccount)+"\t"+str(Gcount)+"\t"+str(Tcount)+"\t"+str(PA)+"\t"+str(PC)+"\t"+str(PG)+"\t"+str(PT)+"\t"+str(ShanJ)+"\t"+str(CovJ)+"\t"+Con)
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
	for i in range(0, sln):
		mvector.append(mSh)
		vvector.append(mSh+vSh)
		v2vector.append(mSh+v2Sh)
		v3vector.append(mSh+v3Sh)
		v5vector.append(mSh+v5Sh)
		v10vector.append(mSh+v10Sh)
		qua_lim_arr.append(30)
		cov_lim_arr.append(100)
	plot(range(sln), Shan, 'b', range(sln), mvector, 'g--', range(sln), v3vector, 'r--', range(sln), v5vector, 'y--', range(sln), v10vector, 'k--')
	title(sample_name+" "+gene+" Shannon enthropy, plus 3, 5, 10 stdev")
	savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Shannon.svg')
	savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Shannon.png')
	clf()
	cla()
	print("Shannon graph ready\n")
#Рисование графика покрытия
	print("Drawing coverage graph\n")
	m = round(statistics.mean(Cov), 2)
	plot(range(0, sln), Cov, 'b', range(0, sln), cov_lim_arr, 'r--')
	title(sample_name+" "+gene+", mean coverage = "+str(m))
	savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Coverage.svg')
	savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Coverage.png')
	print("Coverage graph ready\n")
	clf()
	cla()
#Рисование графика качества прочтения
	print("Drawing quality graph\n")
	qm = round(statistics.mean(qua_arr), 2)
	plot(range(0, sln), qua_arr, 'b', range(0, sln), qua_lim_arr, 'r--')
	title(sample_name+" "+gene+", mean quality = "+str(qm))
	savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Quality.svg')
	savefig(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_Quality.png')
	print("Quality graph ready\n")
	clf()
	cla()
	t.close()
	ans = []
	ans.append(len(cons))
	ans.append(m)
	ans.append(prob_count)
	ans.append(qm)
	return ans

def parse_pileup2(pileup_file, step):
	#Получение массива частот встречаемости нуклеотидов в позиции из файла PILEUP
	print("Parsing pileup data\n")
	with open(pileup_file, 'r') as p:
		p_arr = p.readlines()
	cnt = []
	i = 0
	while i < len(p_arr):
		el = p_arr[i].split()
		if len(el) < 5:
			i +=1
			continue
		else:
			ref = el[2]
			cov = int(el[3])
			seq_str = el[4]
			qua_str = el[5]
			j = 0
			Acount = 0
			Ccount = 0
			Gcount = 0
			Tcount = 0
			dotcount = 0
			gapcount = 0
			inscount = 0
			dels = []
			inss = []
			while j < len(seq_str):
				#Не рассматриваем другие символы
				if seq_str[j] not in ["A", "a", "C", "c", "G", "g", "T", "t", ".", ",", "-", "+"]:
					j +=1
					continue
				#Считаем встречаемость каждого из нуклеотидов
				if seq_str[j] == "A" or seq_str[j] == "a":
					Acount +=1
					j +=1
					continue
				if seq_str[j] == "C" or seq_str[j] == "c":
					Ccount +=1
					j +=1
					continue
				if seq_str[j] == "G" or seq_str[j] == "g":
					Gcount +=1
					j +=1
					continue
				if seq_str[j] == "T" or seq_str[j] == "t":
					Tcount +=1
					j +=1
					continue
				if seq_str[j] == "." or seq_str[j] == ",":
					dotcount +=1
					j +=1
					continue
				#Считаем делеции
				if seq_str[j] == "-":
					gapcount +=1
					j +=1
					num_str = ""
					while seq_str[j] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
						num_str = num_str + seq_str[j]
						j +=1
					if len(num_str) > 0:
						num = int(num_str)
					else:
						num = 0
					dels.append(num)
					j +=num
					continue
				#Считаем вставки
				if seq_str[j] == "+":
					inscount +=1
					j +=1
					num_str = ""
					while seq_str[j] in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
						num_str = num_str + seq_str[j]
						j +=1
					if len(num_str) > 0:
						num = int(num_str)
					else:
						num = 0
					for k in range(num):
						if len(inss) < k+1:
							el_inss = {}
							el_inss["A"] = 0
							el_inss["C"] = 0
							el_inss["G"] = 0
							el_inss["T"] = 0
							inss.append(el_inss)
						if seq_str[j+k] == "A" or seq_str[j+k] == "a":
							inss[k]["A"] +=1
						if seq_str[j+k] == "C" or seq_str[j+k] == "c":
							inss[k]["C"] +=1
						if seq_str[j+k] == "G" or seq_str[j+k] == "g":
							inss[k]["G"] +=1
						if seq_str[j+k] == "T" or seq_str[j+k] == "t":
							inss[k]["T"] +=1
					j +=num
					continue
			if ref == "A" or ref == "a":
				Acount +=dotcount
			if ref == "C" or ref == "c":
				Ccount +=dotcount
			if ref == "G" or ref == "g":
				Gcount +=dotcount
			if ref == "T" or ref == "t":
				Tcount +=dotcount
			#Формируем итоговый массив
			el_cnt = []
			el_cnt.append(Acount)
			el_cnt.append(Ccount)
			el_cnt.append(Gcount)
			el_cnt.append(Tcount)
			el_cnt.append(gapcount)
			el_cnt.append(inscount)
			el_cnt.append(cov)
			qua_arr = []
			for s in qua_str:
				quas = ord(s)
				qua_arr.append(quas)
			qua = statistics.mean(qua_arr) - 33
			el_cnt.append(qua)
			cnt.append(el_cnt)
			i +=1
			#Учитываем делеции
			if gapcount > float(cov) * 0.5:
				dl = round(statistics.mean(dels))
				if step != '':
					i +=dl
					continue
				else:
					el_cnt = [0, 0, 0, 0, gapcount, 0, gapcount, 0]
					cnt.append(el_cnt)
			#Учитываем вставки
			if inscount > float(cov) * 0.5:
				in_cov_arr = []
				for v in inss:
					in_cov = v["A"] + v["C"] + v["G"] + v["T"]
					in_cov_arr.append(in_cov)
				in_cov = round(statistics.mean(in_cov_arr))
				for ins_pos in inss:
					if ins_pos["A"] + ins_pos["C"] + ins_pos["G"] + ins_pos["T"] > in_cov * 0.75:
						el_cnt = []
						el_cnt.append(ins_pos["A"])
						el_cnt.append(ins_pos["C"])
						el_cnt.append(ins_pos["G"])
						el_cnt.append(ins_pos["T"])
						el_cnt.append(0)
						el_cnt.append(ins_pos["A"] + ins_pos["C"] + ins_pos["G"] + ins_pos["T"])
						el_cnt.append(ins_pos["A"] + ins_pos["C"] + ins_pos["G"] + ins_pos["T"])
						el_cnt.append(0)
						cnt.append(el_cnt)
	return cnt

def html_report(report_file, report_arr):
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
	samtools_cmd1 = samtools_path+'samtools view -bSu -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.sam > '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam'
	samtools_cmd2 = samtools_path+'samtools sort -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step
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
	samtools_cmd6 = samtools_path+'samtools sort -@ '+str(os.cpu_count())+' -l 9 '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step
	samtools_cmd7 = samtools_path+'samtools index -@ '+str(os.cpu_count())+' '+sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam'
	os.system(samtools_cmd5)
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam')
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'.bam.bai')
	os.system(samtools_cmd6)
	os.system(samtools_cmd7)
	print("Samtools finished\n")
	os.remove(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_'+step+'_unsorted.bam')
	print("Remove temporary files for sample "+sample_number+" and gene "+gene+"\n")


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
	report = {}
	report['number'] = sample_number
	report['name'] = sample_name
	report['plate'] = plate_number
	report['start_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	copy_data(sample_number, plate_number)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	sample_stats = open(sample_number+os.sep+sample_number+"_stats.txt", "w")
	sample_stats.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	trim_data(sample_number, trimmomatic_path)
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	ttl = conc_data(sample_number)
	report['read_number'] = ttl
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	fastqc_data(sample_number, fastqc_path)
	fastqc_zip = zipfile.ZipFile(sample_number+os.sep+sample_number+'_fastqc.zip', 'r')
	source = fastqc_zip.open(sample_number+'_fastqc/Images/per_base_quality.png')
	target = open(sample_number+os.sep+sample_number+'_per_base_quality.png', "wb")
	shutil.copyfileobj(source, target)
	target.close()
	qfile = open(sample_number+os.sep+sample_number+'_per_base_quality.png', 'rb')
	qfiledata = qfile.read()
	qdata = base64.b64encode(qfiledata)
	qstring = qdata.decode()
	report['read_quality'] = qstring
	qfile.close()	
	second_run = []
	mapping_arr = {}
	print("First run for sample "+sample_number+"\n")
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	for gene in gene_list:
		print(gene+"\n")
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		if not os.path.isdir(sample_number+os.sep+gene):
			os.mkdir(sample_number+os.sep+gene)
		index_bwa(reference_path+gene+'.fasta')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		run_bwa(reference_path+gene+'.fasta', sample_number, gene, '1')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		run_samtools2(sample_number, gene, '1')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
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
		make_pileup(sample_number, gene, reference_path+gene+'.fasta', '1')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		cnt = parse_pileup2(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_1.txt', '1')
		#cnt_file = open(sample_number+os.sep+gene+os.sep+sample_number+"_"+gene+"_cnt.txt", "w") #!!!!!!!!!!!!!!!!!!
		print("Length of gene "+gene+" of sample "+sample_number+" after first run: "+str(len(cnt))+"\n") #!!!!!!!!!!!!!!!!!!
		#for el in cnt: #!!!!!!!!!!!!!!!!!!
		#	cnt_file.write(str(el)+"\n") #!!!!!!!!!!!!!!!!!!
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		get_cons(cnt, reference_path+gene+'.fasta', sample_number, gene, sample_name, '1')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		index_bwa(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_1.fasta')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		run_bwa(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_1.fasta', sample_number, gene, '')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		run_samtools(sample_number, gene, '')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		mapped1 = parse_flagstat(sample_number, gene, '')
		mapping_arr[gene] = mapped1[0]
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		make_pileup(sample_number, gene, reference_path+gene+'.fasta', '')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		cnt1 = parse_pileup2(sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_.txt', '')
		print("Length of gene "+gene+" of sample "+sample_number+" after second run: "+str(len(cnt1))+"\n") #!!!!!!!!!!!!!!!!!!
		#cnt_file.write("\n\n\n\n\n\n\n\n") #!!!!!!!!!!!!!!!!!!
		#for el in cnt1: #!!!!!!!!!!!!!!!!!!
		#	cnt_file.write(str(el)+"\n") #!!!!!!!!!!!!!!!!!!
		#cnt_file.close() #!!!!!!!!!!!!!!!!!!
		get_cons(cnt1, sample_number+os.sep+gene+os.sep+sample_number+'_'+gene+'_1.fasta', sample_number, gene, sample_name, '')
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
		cv = shannon_consensus(cnt1, sample_number, gene, sample_name)
		cover_ar[gene] = cv
		print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
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
			mapping_table.write("< than 1 percent of total reads\n")
			sample_stats.write("< than 1 percent of total reads\n")
	sample_stats.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	mapping_table.write(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	sample_stats.close()
	print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
	report['finish_time'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	mapping_table.close()
print(datetime.now().strftime("%Y-%m-%d %H:%M:%S")+"\n")
print("Finish")
