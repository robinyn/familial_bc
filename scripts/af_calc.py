import sys

sample_list_path = sys.argv[1]
variant_list_path = sys.argv[2]
per_sample_file_path = sys.argv[3]

sample_dict = dict()

with open(sample_list_path, "r") as sample_list:
    for line in sample_list:
        sample_dict[line.strip()]=dict()

with open(variant_list_path, "r") as variant_list:
    for line in variant_list:

        variant = line.strip()

        chr = variant.split("-")[0]
        pos = variant.split("-")[1]
        ref = variant.split("-")[2]
        alt = variant.split("-")[3]

        if len(ref) != 1:
            continue

        for alt_allele in alt.split(","):
            if len(alt_allele) != 1:
                continue

            variant_recoded = "{}-{}-{}-{}".format(chr, pos, ref, alt_allele)

            for sample in sample_dict.keys():
                sample_dict[sample][variant_recoded]=0

with open("test.txt", "w") as output:
    output.write(sample_dict)