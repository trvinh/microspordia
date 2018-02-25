# -*- coding: utf-8 -*-
import sys
import getopt
import glob
import time
from bs4 import BeautifulSoup
import urllib2
import itertools
# import unicodedata	# for converting unicode to ascii

### fetch xml (KGML format) files from KEGG to get connections between KOs
### 2018.01.11
### http://www.kegg.jp/kegg-bin/download?entry=ko02010&format=kgml

def appendToDict(key,value,aDict):
	if not key in aDict.keys():
		aDict[key] = []
		aDict[key].append(value)
	else:
		if not value in aDict[key]:
			aDict[key].append(value)

def main(argv):
	inFile = ''
	try:
		opts, args = getopt.getopt(argv,"i:h",["inFile","help"])
	except getopt.GetoptError:
		print('getKeggInteraction.py -i pathways.list')
		sys.exit(2)

	for opt,arg in opts:
		if opt in ('-h','--help'):
			print('getKeggInteraction.py -i pathways.list')
			sys.exit()
		elif opt in ('-i','--inFile'):
			inFile = arg

	logFile = open("log.cxn","w")
	with open(inFile) as fp:
		for line in fp:
			pathID = line.rstrip().split("\t")[0]

			url = "http://rest.kegg.jp/get/ko"+pathID+"/kgml"
			print(url)

			try:
				##### open url if existing and create output file
				xmlFile = urllib2.urlopen(url)

				##### read file into beatifulsoup object
				xmlIn = BeautifulSoup(xmlFile,"xml")

				##### PARSING XML FILE

				### get list of all KOs
				id2ko = {}	# id2ko{(nodeID:koID)}

				for entry in xmlIn.findAll("entry"):
					if entry.get("type") == "ortholog":
						koIDs = entry.get("name")
						entryID = entry.get("id")
						id2ko[entryID] = koIDs

				### get links between KOs
				flag = 0
				pairKO = {}	# pairKO{(Source\tTarget:Type)}
				for relation in xmlIn.findAll("relation"):
					entry1 = relation.get("entry1")
					entry2 = relation.get("entry2")
					relType = relation.get("type")
					if(entry1 in id2ko and entry2 in id2ko):
						for ko1 in id2ko[entry1].split(" "):
							for ko2 in id2ko[entry2].split(" "):
								if(ko1.replace("ko:","") != ko2.replace("ko:","")):
									# print(ko1.replace("ko:","")+"\t"+ko2.replace("ko:","")+"\t"+relType)
									pairTMP = [ko1.replace("ko:","").encode('ascii','ignore'),ko2.replace("ko:","").encode('ascii','ignore')]
									pairTMP.sort()
									pairKO["\t".join(pairTMP)] = relType
									# output.write(ko1.replace("ko:","")+"\t"+ko2.replace("ko:","")+"\t"+relType+"\n")
									flag = 1
					# else:
					# 	print(entry1+" or "+entry2+" not in list")
					# time.sleep(2)

				for group in xmlIn.findAll("entry"):
					if group.get("type") == "group":
						allCom = []
						for component in group.findAll("component"):
							componentID = component.get("id").encode("ascii")
							allCom.append(componentID)

						for pair in list(itertools.combinations(allCom, 2)):
							if(pair[0] in id2ko and pair[1] in id2ko):
								# time.sleep(4)
								for ko1 in id2ko[pair[0]].split(" "):
									for ko2 in id2ko[pair[1]].split(" "):
										# print(ko1.replace("ko:","")+"\t"+ko2.replace("ko:","")+"\tgroup")
										if(ko1.replace("ko:","") != ko2.replace("ko:","")):
											pairTMP = [ko1.replace("ko:","").encode('ascii','ignore'),ko2.replace("ko:","").encode('ascii','ignore')]
											pairTMP.sort()
											pairKO["\t".join(pairTMP)] = relType
											# output.write(ko1.replace("ko:","")+"\t"+ko2.replace("ko:","")+"\tgroup"+"\n")
											flag = 1

				if flag == 0:
					logFile.write(pathID+"\tno-relation\t"+line)
				else:
					output = open(pathID+".cxn","w")
					for pairID in pairKO:
						output.write(pairID+"\t"+pairKO[pairID]+"\n")
					output.close()
			except urllib2.HTTPError, e:
			    logFile.write(url+"\tnot-exist\t"+line)
			except urllib2.URLError, e:
			    logFile.write(url+"\tnot-exist\t"+line)

	print("Finished! Check log.cxn for missing information!")
	logFile.close()
if __name__ == "__main__":
	if len(sys.argv[1:])==0:
		print('orthoxmlParser.py -i input')
		sys.exit(2)
	else:
		main(sys.argv[1:])
