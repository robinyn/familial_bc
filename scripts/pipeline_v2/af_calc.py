import sys

HEBCS_samples_path = sys.argv[1]
per_sample_path = sys.argv[2]
output_file_path = sys.argv[3]

variant_dict = dict()
HEBCS_samples = []
tri_allelic = []

sample_count = 3403

header = "variant\thom_ref\thet\thom_alt\tfound\tnot_found\tAF_ref\tAF_alt\n"

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

        if sampleID in HEBCS_samples:
            continue

        if "," in alt_allele:
            tri_allelic.append("\t".join(line))
            continue

        variant = "{}-{}-{}-{}".format(chrom, pos, ref_allele, alt_allele)

        print("Sample: {}\nVariant: {}".format(sampleID, variant))

        if variant not in variant_dict.keys():
            variant_dict[variant]={"hom_ref":0, "het":0, "hom_alt":0, "found":0, "not_found":0, "AF_ref":0, "AF_alt":0}

        variant_dict[variant]["found"]+=1

        homozygous_ref = "{}{}".format(ref_allele, ref_allele)
        heterozygous_1 = "{}{}".format(ref_allele, alt_allele)
        heterozygous_2 = "{}{}".format(alt_allele, ref_allele)
        homozygous_alt = "{}{}".format(alt_allele, alt_allele)

        if genotype==homozygous_ref:
            variant_dict[variant]["hom_ref"]+=1
        elif genotype==heterozygous_1 or genotype==heterozygous_2:
            variant_dict[variant]["het"]+=1
        elif genotype==homozygous_alt:
            variant_dict[variant]["hom_alt"]+=1
        else:
            print("Invalid genotype: {}".format(genotype))
            exit()

for variant in variant_dict.keys():
    print(variant)

    variant_dict[variant]["not_found"] = sample_count - variant_dict[variant]["found"]
    variant_dict[variant]["AF_ref"]=((2*(variant_dict[variant]["hom_ref"] + variant_dict[variant]["not_found"]))+(variant_dict[variant]["het"]))/(sample_count*2)
    variant_dict[variant]["AF_alt"]=((2*(variant_dict[variant]["hom_alt"]))+(variant_dict[variant]["het"]))/(sample_count*2)

with open(output_file_path, "w") as output_file:
    output_file.write(header)

    for variant in variant_dict.keys():
        output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(variant, \
                                                                variant_dict[variant]["hom_ref"], \
                                                                variant_dict[variant]["het"], \
                                                                variant_dict[variant]["hom_alt"], \
                                                                variant_dict[variant]["found"], \
                                                                variant_dict[variant]["not_found"], \
                                                                variant_dict[variant]["AF_ref"], \
                                                                variant_dict[variant]["AF_alt"]))

with open("tri_allelic.tsv", "w") as output_file:
    output_file.write("\n".join(tri_allelic))


