#!/usr/bin/env python3

# ---------------------------------------
# Universal Plasmid Maker
# ---------------------------------------

# ---------- DEFAULT PLASMID BACKBONE ----------
# Broad host range (IncQ-like)
ORI_V = "TTGACATGTTGACATGTTGACATG"
REP_A = "ATGGCTGCTGCTGCTGCTGCT"
REP_B = "ATGCCGCCGCCGCCGCCGCC"
REP_C = "ATGAAAGGGAAAGGGAAAGGG"

BACKBONE = ORI_V + REP_A + REP_B + REP_C

# ---------- RESTRICTION SITES ----------
RESTRICTION_SITES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG"
}

# ---------- ANTIBIOTIC GENES ----------
ANTIBIOTIC_GENES = {
    "Ampicillin": "ATGAGTATTCAACATTTCCGTG",
    "Kanamycin": "ATGAGCCATATTCAACGGGAA",
    "Chloramphenicol": "ATGGAGAAAAAAATCACTGG"
}

# ---------- READ FASTA ----------
def read_fasta(filename):
    seq = ""
    with open(filename) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq

# ---------- READ DESIGN FILE ----------
def read_design(filename):
    mcs = []
    antibiotics = []

    with open(filename) as f:
        for line in f:
            name, value = line.strip().split(",")
            value = value.strip()
            if "Cloning" in name:
                mcs.append(value)
            else:
                antibiotics.append(value)
    return mcs, antibiotics

# ---------- BUILD MCS ----------
def build_mcs(enzymes):
    seq = ""
    for enzyme in enzymes:
        seq += RESTRICTION_SITES.get(enzyme, "")
    return seq

# ---------- BUILD ANTIBIOTIC CASSETTE ----------
def build_antibiotic(antibiotics):
    seq = ""
    for ab in antibiotics:
        seq += ANTIBIOTIC_GENES.get(ab, "")
    return seq

# ---------- MAIN ----------
def main():
    insert_seq = read_fasta("Input.fa")
    mcs_enzymes, antibiotics = read_design("Design.txt")

    mcs_seq = build_mcs(mcs_enzymes)
    antibiotic_seq = build_antibiotic(antibiotics)

    plasmid_seq = BACKBONE + mcs_seq + insert_seq + antibiotic_seq

    with open("Output.fa", "w") as f:
        f.write(">Universal_Plasmid\n")
        for i in range(0, len(plasmid_seq), 60):
            f.write(plasmid_seq[i:i+60] + "\n")

    print("Plasmid construction completed.")
    print("Output written to Output.fa")

if __name__ == "__main__":
    main()
