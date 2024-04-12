import sys

# wget -qO- https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz | zcat | buffer -m 100M -s 100k |  python3 a.py list_of_variants.tsv > tbl.out

variants = set()
populations = ["afr",
               "amr",
               "asj",
               "eas_jpn",
               "eas_kor",
               "eas_oea",
               "fin",
               "nfe_seu",
               "nfe_onf",
               "nfe_est",
               "nfe_swe",
               "nfe_bgr",
               "nfe_nwe",
               "sas",
               "oth"]

controls_maxAN = 54704 * 2

header = ["variant\tref\talt",
          "QUAL",
          "maf",
          "AC",
          "call_rate"]
#header += ["<dp10",
#           "<dp20",
#           "<gq20",
#           "<gq30",
#           "GQ_bins",
#           "DP_bins"] + populations

print("\t".join(header))

variants_file = sys.argv[1]
with open(variants_file, "r") as fd:
    for line in fd:
        variants.add(line.strip())

for line in sys.stdin:
    if line.startswith("#"):
        continue
    tokens = line.rstrip().split("\t")
    if tokens[6] != "PASS":
        continue
    variant = "chr" + tokens[0] + ":" + tokens[1] + "\t" + tokens[3] + "\t" + tokens[4]
    if variant not in variants:
        continue
    qual = tokens[5]
    info = dict()
    for record in tokens[7].split(";"):
        if "=" not in record:
            continue
        key, value = record.split("=")
        info[key] = value

    popinfo = []
    for pop in populations:
        ac = int(info["controls_AC_" + pop])
        an = int(info["controls_AN_" + pop])
        homalt = int(info["controls_nhomalt_" + pop])
        het = ac - 2 * homalt
        homref = an // 2 - homalt - het
        popinfo.append((homref, het, homalt))

    sum_an = sum(map(lambda x: sum(x), popinfo)) * 2
    if sum_an == 0:
        continue
    sums = list(map(lambda i: sum(map(lambda x: x[i], popinfo)), range(0, 3)))
    minor = min(sums[0], sums[2])

    mac = minor * 2 + sums[1]
    maf = mac / sum_an
    if minor == 0 and sums[1] ==  0:
        continue
    cr = sum_an / controls_maxAN

    dp_hist = list(map(int, info["dp_hist_all_bin_freq"].split("|")))
    gq_hist = list(map(int, info["gq_hist_all_bin_freq"].split("|")))
    dp_sum = sum(dp_hist)
    gq_sum = sum(gq_hist)

    dp10 = str(sum(dp_hist[0:2]) / dp_sum)
    gq20 = str(sum(gq_hist[0:4]) / gq_sum)
    dp20 = str(sum(dp_hist[0:4]) / dp_sum)
    gq30 = str(sum(gq_hist[0:6]) / gq_sum)

    row = [variant,
           qual,
           str(maf),
           str(mac),
           str(cr)]
    #row += [dp10,
    #       dp20,
    #       gq30,
    #       info["gq_hist_all_bin_freq"],
    #       info["dp_hist_all_bin_freq"]] + \
    #       list(map(lambda x: str(x[0]) + ";" + str(x[1]) + ";" + str(x[2]), popinfo))
    print("\t".join(row))
