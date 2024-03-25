import sys

file_dir = sys.argv[1]

with open(file_dir, "r") as input_vcf:
    for line in input_vcf:
        if line.startswith("#"):
            continue

        line = line.strip().split("\t")

        if "," in line[4]:
            print(line)
