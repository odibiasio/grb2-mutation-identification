# Directory to FASTA files for csv transfer
Homo_Sapiens = "/Users/livdibiasio/Desktop/Grb2_data/HomoSapiens_Grb2.fasta"
Rat = "/Users/livdibiasio/Desktop/Grb2_data/Rat_Grb2.fasta"
Mouse = "/Users/livdibiasio/Desktop/Grb2_data/Mouse_Grb2.fasta"
Chicken = "/Users/livdibiasio/Desktop/Grb2_data/Chicken_Grb2.fasta"
African_Clawed_Frog = "/Users/livdibiasio/Desktop/Grb2_data/AfricanClawedFrog_Grb2.fasta"
Zebrafish = "/Users/livdibiasio/Desktop/Grb2_data/Zebrafish_Grb2.fasta"

# Checking if absolute file paths work
print()
print(" { CHECK }   ")
print()
import os 
if os.path.exists(Homo_Sapiens):
	print("human path exists")
if os.path.exists(Rat):
	print("rat path exists")
if os.path.exists(Mouse):
	print("mouse path exists")
if os.path.exists(Chicken):
	print("chicken path exists")
if os.path.exists(African_Clawed_Frog):
	print("african clawed frog path exists")
if os.path.exists (Zebrafish):
	print("zebrafish path exists")
	print()
	print(" { SEQUENCES }   ")
	print()

# Printing out FASTA files as parsed sequences
	# Homo Sapiens
	from Bio import SeqIO
	for record in SeqIO.parse(Homo_Sapiens,"fasta"):
		# printing out amino acid sequence
			print("Homo Sapiens")
			print("Protein ID: ", record.id)
			print("Sequence: ", record.seq) 
		# identify and print the SH2 and SH3 domains
			Homo_Sapiens_List = str(record.seq)
			print()
			Homo_Sapiens_domains = {
				"SH2": Homo_Sapiens_List[59:151],
				"SH3_N": Homo_Sapiens_List[0:57],
				"SH3_C": Homo_Sapiens_List[155:214]
			}
			print("SH2 Domain: ", Homo_Sapiens_domains["SH2"])
			print ("SH3 Domain: ", Homo_Sapiens_domains["SH3_N"] + Homo_Sapiens_domains["SH3_C"])
			print("-"*50)
	# Rat
	from Bio import SeqIO
	for record in SeqIO.parse(Rat,"fasta"):
		# printing out amino acid sequence
			print("Rat")
			print("Protein ID: ", record.id)
			print("Sequence: ", record.seq)
		# identify and print the SH2 and SH3 domains
			Rat_List = str(record.seq)
			print()
			Rat_domains = {
				"SH2": Rat_List[59:151],
				"SH3_N": Rat_List[0:57],
				"SH3_C": Rat_List[155:214]
			}
			print("SH2 Domain: ", Rat_domains["SH2"])
			print ("SH3 Domain: ", Rat_domains["SH3_N"] + Rat_domains["SH3_C"])
			print("-"*50)

	# Mouse
	from Bio import SeqIO
	for record in SeqIO.parse(Mouse,"fasta"):
		# printing out amino acid sequence
			print("Mouse")
			print("Protein ID: ", record.id)
			print("Sequence: ", record.seq)
		# identify and print the SH2 and SH3 domains
			Mouse_List = str(record.seq)
			print()
			Mouse_domains = {
				"SH2": Mouse_List[59:151],
				"SH3_N": Mouse_List[0:57],
				"SH3_C": Mouse_List[155:214]
			}
			print("SH2 Domain: ", Mouse_domains["SH2"])
			print ("SH3 Domain: ", Mouse_domains["SH3_N"] + Mouse_domains["SH3_C"])
			print("-"*50)

	# Chicken
	from Bio import SeqIO
	for record in SeqIO.parse(Chicken,"fasta"):
		# printing out amino acid sequence
			print("Chicken")
			print("Protein ID: ", record.id)
			print("Sequence: ", record.seq)
		# identify and print the SH2 and SH3 domains
			Chicken_List = str(record.seq)
			print()
			Chicken_domains = {
				"SH2": Chicken_List[59:151],
				"SH3_N": Chicken_List[0:57],
				"SH3_C": Chicken_List[155:214]
			}
			print("SH2 Domain: ", Chicken_domains["SH2"])
			print ("SH3 Domain: ", Chicken_domains["SH3_N"] + Chicken_domains["SH3_C"])
			print("-"*50)

	# African Clawed Frog
	from Bio import SeqIO
	for record in SeqIO.parse(African_Clawed_Frog,"fasta"):
		# printing out amino acid sequence
			print("African Clawed Frog")
			print("Protein ID: ", record.id)
			print("Sequence: ", record.seq)
		# identify and print the SH2 and SH3 domains
			African_Clawed_Frog_List = str(record.seq)
			print()
			African_Clawed_Frog_Domains= {
				"SH2": African_Clawed_Frog_List[59:151],
				"SH3_N": African_Clawed_Frog_List[0:57],
				"SH3_C": African_Clawed_Frog_List[167:226]
			}
			print("SH2 Domain: ", African_Clawed_Frog_Domains["SH2"])
			print ("SH3 Domain: ", African_Clawed_Frog_Domains["SH3_N"] + African_Clawed_Frog_Domains["SH3_C"])
			print("-"*50)

	# Zebrafish
	from Bio import SeqIO
	for record in SeqIO.parse(Zebrafish,"fasta"):
		# printing out amino acid sequence
			print("Zebrafish")
			print("Protein ID: ", record.id)
			print("Sequence: ", record.seq)
		# identify and print the SH2 and SH3 domains
			Zebrafish_List = str(record.seq)
			print()
			Zebrafish_Domains = {
				"SH2": Zebrafish_List[59:151],
				"SH3_N": Zebrafish_List[0:57],
				"SH3_C": Zebrafish_List[155:214]
			}
			print("SH2 Domain: ", Zebrafish_Domains["SH2"])
			print ("SH3 Domain: ", Zebrafish_Domains["SH3_N"] + Zebrafish_Domains["SH3_C"])
			print("-"*50)

# Sequence Anaylsis 
seqs = {
    "Homo_Sapiens": "MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIEMKPHPWFFGKIPRAKAEEMLSKQRHDGAFLIRESESAPGDFSLSVKFGNDVQHFKVLRDGAGKYFLWVVKFNSLNELVDYHRSTSVSRNQQIFLRDIEQVPQQPTYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYVTPVNRNV",
    "Rat": "MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIEMKPHPWFFGKIPRAKAEEMLSKQRHDGAFLIRESESAPGDFSLSVKFGNDVQHFKVLRDGAGKYFLWVVKFNSLNELVDYHRSTSVSRNQQIFLRDIEQVPQQPTYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYVTPVNRNV",
    "Mouse": "MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKDGFIPKNYIEMKPHPWFFGKIPRAKAEEMLSKQRHDGAFLIRESESAPGDFSLSVKFGNDVQHFKVLRDGAGKYFLWVVKFNSLNELVDYHRSTSVSRNQQIFLRDIEQMPQQPTYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYVTPVNRNV",
    "Chicken": "MEAIAKYDFKATADDELSFKRGDILKVLNEECDQNWYKAELNGKGGFIPKNYIEMKPHPWFFGKIPRAKAEEMLGKQRHDGAFLIRESESAPGDFSLSVKFGNDVQQFKVLRDGAGKYLLWVVKFNSLNELVDYHRSTSVSRNQQIFLRDIEQVPQQPTYVQALFDFDPQEEGELGFRRGDFIQVLDNSDPNWWKGACHGQTGMFPRNYVTPVNRNI",
    "Frog": "MEAIAKYDFKATADDELSFKRGDVLKVLNEECDQNWYKAELNGKDGFIPKNYIEMKAHPWFFGKIPRAKAEEMLGKQRHDGAFLIRESESAPGDFSLSVKFGNDVQHFKVLRDGAGKYFLWVVKFNSLNELVDYHRSTSVSRNQQIFLRDIEQVPQVHGGDRATSLPQQPTYVQALFDFDPQEDGELGFRRGDFIQVVDNSDPNWWKGTCLSQTGMFPRNYVTPVNRNM",
    "Zebrafish": "MEAIAKYDFKATADDELSFKRGEILKVLNEECDQNWYKAELNGKEGFIPKNYIEMKPHPWFYGKIPRAKAEEMLNKQRHDGAFLIRESESAPGDFSLSVKFGNDVQHFKVLRDGAGKYFLWVVKFNSLNSLVDYHRSTSVSRNQPIFLRDIEQVPQNSTYVQALFDFDPQEDGELGFRRGDFIQVLDNSDPNWWKGACHGQTGMFPRNYVTPVNQNM",
}

#END RESULT: rank all species vs Homo_Sapiens using a point system

print()
print(" { RANKING }   ")
print()

def identity_points(human_seq: str, other_seq: str):
    # compare only the overlapping region so different lengths won't crash
    min_len = min(len(human_seq), len(other_seq))
    human_seq = human_seq[:min_len]
    other_seq = other_seq[:min_len]

    points = 0
    for a, b in zip(human_seq, other_seq):
        if a == b:
            points += 1

    percent = (points / min_len) * 100
    return points, percent, min_len


# 1) Get the human reference sequence
human_seq = seqs["Homo_Sapiens"]

# 2) Score everyone vs human
scores = {}
for species, seq in seqs.items():
    if species == "Homo_Sapiens":
        continue  # skip self comparison
    pts, pct, similarity = identity_points(human_seq, seq)
    scores[species] = (pts, pct, similarity)

# 3) Rank highest percent identity to lowest
ranked = sorted(scores.items(), key=lambda x: x[1][1], reverse=True)

# 4) Print results
print("\nRanked vs Homo_Sapiens (Whole Protein)")
print("-"*10)
print("Species        Points    Possible   Percent")
for species, (pts, pct, similarity) in ranked:
    print(f"{species:12s} {pts:6d}   {similarity:10d}   {pct:6.2f}%")

#DOMAIN-BASED RANKING (SH2 and SH3)

# Domain boundaries 
SH2 = slice(59, 151)
SH3_N = slice(0, 57)
SH3_C = slice(155, 214)

def get_sh3(seq: str) -> str:
    return seq[SH3_N] + seq[SH3_C]


def rank_domain(domain_name: str, human_domain: str, domain_func):
    domain_scores = {}

    for species, seq in seqs.items():
        if species == "Homo_Sapiens":
            continue

        pts, pct, similarity = identity_points(human_domain, domain_func(seq))
        domain_scores[species] = (pts, pct, similarity)

    ranked_domain = sorted(domain_scores.items(), key=lambda x: x[1][1], reverse=True)

    print(f"\nRanked vs Homo_Sapiens ({domain_name})")
    print("-"*10)
    print("Species        Points     Possible   Percent")
    for species, (pts, pct, similarity) in ranked_domain:
        print(f"{species:12s} {pts:6d}   {similarity:10d}   {pct:6.2f}%")
    print()


#Run SH2 ranking
rank_domain(
    "SH2 Domain",
    human_seq[SH2],
    lambda s: s[SH2]
)

# Run SH3 ranking (N + C combined)
rank_domain(
    "SH3 Domains (N + C)",
    get_sh3(human_seq),
    get_sh3
)

