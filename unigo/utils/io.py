def loadUniprotIDsFromCliFiles(expressedProtIDFile, deltaProtIDFile):
	expUniprotID   = []
	deltaUniprotID = []
	with open(expressedProtIDFile) as f:
		for l in f:
			expUniprotID.append(l.rstrip())
	with open(deltaProtIDFile) as f:
		for l in f:
			deltaUniprotID.append(l.rstrip())

	if not check_proteins_subset(expUniprotID, deltaUniprotID):
		raise Exception("Differentially expressed proteins are not completely included in total proteins")

	return expUniprotID, deltaUniprotID

def check_proteins_subset(major_list:list[str], sublist:list[str]) -> bool:
	return set(sublist).issubset(major_list)
