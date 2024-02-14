import sys

HEBCS_samples_path = sys.argv[1]
per_sample_path = sys.argv[2]

variant_dict = dict()
HEBCS_samples = []

sample_count = 110752

header = "variant\thom_ref\thet\thom_alt\tfound\tnot_found\tAF_ref\tAF_alt"

with open(HEBCS_samples_path, "r") as HEBCS_samples_file:
    HEBCS_samples_file.readline()

    for index, line in enumerate(HEBCS_samples_file):
        line = line.strip().split("\t")

        sampleID = line[0]

        HEBCS_samples.append(sampleID)

with open(per_sample_path, "r") as sample_file:
    sample_file.readline()

    for line in sample_file:
        line = line.strip().split("\t")

        sampleID = line[0]
        variant = line[1]

        variant = variant.split("-")

        chrom = variant[0]
        pos = variant[1]
        ref_allele = variant[2]
        alt_allele = variant[3]

        genotype = line[2]

        homozygous_ref = "{}{}".format(ref_allele, ref_allele)
        heterozygous = "{}{}".format(ref_allele, alt_allele)
        homozygous_alt = "{}{}".format(alt_allele, alt_allele)

        if sampleID in HEBCS_samples:
            continue

        for allele in alt_allele.split(","):
            variant = "{}-{}-{}-{}".format(chrom, pos, ref_allele, allele)

            if variant not in variant_dict.keys():
                variant_dict[variant]={"hom_ref":0, "het":0, "hom_alt":0, "found":0, "not_found":0, "AF_ref":0, "AF_alt":0}

            variant_dict[variant]["found"]+=1

            if genotype==homozygous_ref:
                variant_dict[variant]["hom_ref"]+=1
            elif genotype==heterozygous:
                variant_dict[variant]["het"]+=1
            elif genotype==homozygous_alt:
                variant_dict[variant]["hom_alt"]+=1

for variant in variant_dict.keys():
    print(variant)

    variant_dict[variant]["not_found"] = sample_count - variant_dict[variant]["found"]
    variant_dict[variant]["AF_ref"]=((2*(variant_dict[variant]["hom_ref"] + variant_dict[variant]["not_found"]))+((variant_dict[variant]["het"])/2))/sample_count
    variant_dict[variant]["AF_alt"]=((2*(variant_dict[variant]["hom_alt"] + variant_dict[variant]["not_found"]))+((variant_dict[variant]["het"])/2))/sample_count

with open("allele_frequencies.tsv", "w") as output_file:
    output_file.write(header)

    for variant in variant_dict.keys():
        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(variant, \
                                                                variant_dict[variant]["hom_ref"], \
                                                                variant_dict[variant]["het"], \
                                                                variant_dict[variant]["hom_alt"], \
                                                                variant_dict[variant]["found"], \
                                                                variant_dict[variant]["not_found"], \
                                                                variant_dict[variant]["AF_ref"], \
                                                                variant_dict[variant]["AF_alt"]))





