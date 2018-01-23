#Protein sequence statistics
#www.pavannichani.com

library(Peptides)

#Determines the amino acid composition of the sequence
aaComp(cox1[1])

#Returns the aliphatic index of protein sequence 
aIndex(cox1)

#Predicts the net charge of the protein
charge(seq, pH = 7, pKscale = "Lehninger")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Bjellqvist")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="EMBOSS")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Murray")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Sillero")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Solomon")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Stryer")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Dawson")
#charge(seq="FLPVLAGLTPSIVPKLVCLLTKKC",pH=7, pKscale="Rodwell")
charge(cox1)

#Overall hydrophobicity
hydrophobicity(cox1[1])