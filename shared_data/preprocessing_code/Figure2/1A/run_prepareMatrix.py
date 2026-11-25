#-*- coding:utf-8 -*-
import re,sys,os
from collections import OrderedDict

cancerbar = ["01","02","03","04","05","06","07","08","09"]
patientInfo = {}
sexcov = {}

statusIndex, daysIndex, deathIndex = 0, 0, 0

genderIndex = 0

for line in open(sys.argv[2],"r"):
    if re.search("submitter_id",line):
        info = line.rstrip('\n').split('\t')
        for index,element in enumerate(info):
            if re.search("days_to_last_follow_up",element):
                daysIndex = index + 1
            if re.search("vital_status",element):
                statusIndex = index + 1
            if re.search("days_to_death",element):
                deathIndex = index + 1  
            if re.search("gender",element):
                genderIndex = index + 1
        continue
    # 开始整理病人的数据
    info = line.rstrip("\n").split("\t")
    #print(info[2].split('"')[1])
    patient = info[2].split('"')[1]
    if info[statusIndex] == "NA":
        continue 
    status = info[statusIndex].split('"')[1]
    #print(status)
    daystodeath = info[deathIndex]  ##有明确的死亡时间
    flowupordeath = info[daysIndex]  ## 随访时间、如果死亡时间有, 则时间设置为死亡时间
    days =  0
    if status == "Dead":
        status = 1
    elif status == "Alive":
        status = 0
    else:
        continue
    if daystodeath != "NA":
        days = daystodeath
    else:
        days = flowupordeath
    if days == "NA":  ##时间不明确的舍弃
        continue
    if float(days)<0:  ##有的天数有负号，舍去
        continue
    patientInfo[patient] = "\t".join([str(status),str(days)])
    sexcov[patient] = info[genderIndex] if genderIndex!=0 else "unknown"

print(patientInfo)
print(sexcov)

##### 整理数据
ribp = []

for line in open("../coreAMPs.txt","r"):
    if re.search("^genename",line):
        continue
    ribp.append(line.split('\t')[0])

# 填充基因信息
out = open("geneId.txt","a")
#key	ProbeSet	Symbol	Name	EntrezID
#IMP4	ENSG00000136718	IMP4	IMP U3 Small Nucleolar Ribonucleoprotein 4	92856

ensgToname = OrderedDict()
entona = {}

for line in open("hg19.genecode.V27.gtf", "r"):
    if re.search("^#",line):
        continue
    info = line.split('\t')
    if info[2]=="gene":
        ensg, name = "",""
        for index,element in enumerate(info[8].split(";")):
            if re.search("gene_id",element):
                ensg = element.split(".")[0].split('"')[1]
            if re.search("gene_name",element):
                name = element.split('"')[1]
        entona[name] = ensg
        ensgToname[ensg] = name

entoID = {}
for line in open("hgnc-symbol-check.tsv", "r"):
    if re.search("Approved symbol",line) or re.search("Previous symbol",line):
        info = line.split("\t")
        print(info)
        name = info[0]
        fullname = info[3]
        hgnc = info[4].replace('"', "").split(":")[1]
        entoID[name] = name+'\t'+fullname+'\t'+hgnc

out.write("key"+"\t"+"ProbeSet"+'\t'+"Symbol"+'\t'+"Name"+"\t"+"EntrezID"+'\n')
for key in ribp:
    out.write(key+'\t'+entona[key]+'\t'+entoID[key]+'\n')
out.close()

## 填充表达谱
patientIndex = OrderedDict()
expressionD = OrderedDict()
out = open("exprTrain.txt","a")
for line in open(sys.argv[1],"r"):
    if re.search("TCGA",line):
        for index, element in enumerate(line.rstrip('\n').split('\t')):
            if element == '""':
                continue
            #print(element)
            if element.split("-")[3][0] + element.split("-")[3][1] in cancerbar: ##这个病人是一个tumor
                if "-".join(element.split('"')[1].split("-")[0:3]) in patientInfo:   ##只保存you信息的人
                    patientIndex["-".join(element.split('"')[1].split("-")[0:3])] = index + 1 
        #print(patientIndex)
        continue
    info = line.rstrip('\n').split('\t')
    ensg = info[0].split("\"")[1].split(".")[0]
    if ensg in ensgToname and ensgToname[ensg] in ribp:   #如果是
        expressionD[ensgToname[ensg]] = []
        for key in patientIndex:
            expressionD[ensgToname[ensg]].append(info[patientIndex[key]])

headers = ["geneSymbol"]
for key in patientIndex:
    headers.append(key)
out.write('\t'.join(headers)+'\n')

done = []
for key in ensgToname:
    if ensgToname[key] in ribp and ensgToname[key] not in done:
        done.append(ensgToname[key])
        out.write(ensgToname[key]+'\t'+"\t".join(expressionD[ensgToname[key]])+'\n')
out.close()

# 填充协变量信息
out = open("covTrain.txt","a")
out.write("sex"+'\n')
for key in patientIndex:
    out.write(sexcov[key]+'\t'+'\n')
out.close()

# 填充生存信息
out = open("status.txt","a")
out.write("status"+'\t'+"tme"+'\n')
for key in patientIndex:
    out.write(patientInfo[key]+'\n')
out.close()