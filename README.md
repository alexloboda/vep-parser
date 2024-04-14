# Supporting tools for dnascore rare-variant analysis.

The following scripts can be used to parse gnomad vcf files to extract minor allele frequencies in control population and 
VEP VCF file to annotate list of variants with gene symbols and predicted variant consequences.

## gnomAD VCF file parser

To run the parser use the following command:

    wget -qO- https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.    sites.vcf.bgz | zcat | python3 gnomad.py list_of_variants.tsv > tbl.out

The list of variants (one per line) should be passed in the following format, for example:

    chr1:55039768\tC\tT
    chr1:55039805\tC\tT


## VEP VCF file parser

To run the parser you need to open annotation file and look for the VEP annotation format in the INFO section of
VCF file. The following command will parse the file vep_annotation.vcf with corresponding VEP annotation strcture:

    python3 vep_parser.py annotation.vcf output.tsv "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|REFSEQ_MATCH|SOURCE|GIVEN_REF|USED_REF|BAM_EDIT"
