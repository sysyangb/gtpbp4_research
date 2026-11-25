#-*- coding:utf-8 -*-
import re,sys,os

out = open("GTPBP4.msData.txt","a")

out.write("Accession"+'\t'+"geneName"+'\t'+"IP"+'\t'+"control"+'\t'+"foldchange"+'\n')


accesion = []
accesionToname = {}
ipDicT = {}

for line in open("IP.txt",'r'):
    if re.search("^Checked",line):
        continue
    info = line.split("\t")
    # 去除角蛋白
    if re.search("^Keratin, ",info[3]):
        continue
    if info[2] not in accesion:
        accesion.append(info[2])
    ipDicT[info[2]] = int(info[5])
    accesionToname[info[2]] = info[19]

controlDict = {}
for line in open("Control.txt",'r'):
    if re.search("^Checked",line):
        continue
    info = line.split("\t")
    # 去除角蛋白
    if re.search("^Keratin, ",info[3]):
        continue
    if info[2] not in accesion:
        accesion.append(info[2])
    controlDict[info[2]] = int(info[5])
    accesionToname[info[2]] = info[18]

for key in accesion:
    temp = []
    temp.append(key)
    temp.append(accesionToname[key])
    numIP, numControl = "NA","NA"
    if key in ipDicT:
        numIP = ipDicT[key]
    temp.append(str(numIP))
    if key in controlDict:
        numControl = controlDict[key]
    temp.append(str(numControl))
    if numControl == "NA":
        temp.append("100")
    elif numIP == "NA":
        temp.append("NA")
    else:
        temp.append(str(numIP/numControl))
    out.write("\t".join(temp)+'\n')


out.close()

