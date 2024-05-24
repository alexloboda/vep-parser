import argparse
import re
from enum import Enum
import os
import contextlib

class VepParser:
    vcf_format = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO"

    def __init__(self, format_line):
        self.fields = format_line.split("|")
        self.vcf_fields = self.vcf_format.split("|")

    def parse_file(self, file):
        with open(file) as fd:
            for line in fd:
                if line.startswith("#"):
                    continue
                self.parse_line(line.rstrip())

    def extract_vep(self, info):
        csq = None
        fields = info.split(";")
        for field in fields:
            if field.startswith("CSQ="):
                csq = field[4:]
                break
        if csq is None:
            print("Could not find CSQ field in the INFO column")
            exit(1)
        transcripts = csq.split(",")
        result = []
        for t in transcripts:
            tokens = t.split("|")
            if len(tokens) != len(self.fields):
                print("Number of fields in CSQ does not match the number of fields in the argument")
                print("Info field " + info)
                exit(1)
            result.append(dict(zip(self.fields, tokens)))
        return result

    def parse_line(self, line):
        vcf_fields = line.split("\t")
        variant = vcf_fields[0:7]
        variant_dict = dict((self.vcf_fields[i], variant[i]) for i in range(len(variant)))
        info = vcf_fields[7]
        if info == ".":
            return
        vep = self.extract_vep(info)
        for annotation in vep:
            gene = annotation['SYMBOL']
            if gene == '':
                continue
            if annotation['CANONICAL'] != 'YES':
                continue
            biotype = annotation['BIOTYPE']
            if biotype != 'protein_coding':
                continue
            cs = annotation['Consequence'].split("&")
            if len(cs) == 0:
                continue
            for cons in cs:
                annotation['Consequence'] = cons
                self.process_annotation(annotation, variant_dict)

    def process_annotation(self, annotation: dict, variant: dict):
        raise NotImplementedError()


class UniqueValues(VepParser):
    def __init__(self, vep_format):
        super().__init__(vep_format)
        self.uniq_vals = dict((key, set()) for key in self.fields)

    def process_annotation(self, annotation, _):
        for key in annotation.keys():
            self.uniq_vals[key].add(annotation[key])

    def close(self):
        for k in self.uniq_vals['Consequence']:
            print(k)


class Annotator(VepParser):
    def __init__(self, output, vep_format):
        self.lines = set()
        super().__init__(vep_format)
        self.fd = open(output, "w")
        self.fd.write("chr\tpos\tref\talt\treason\tgene\n")

    def process_annotation(self, annotation, variant):
        # vcf_format = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO"
        line = variant['CHROM'] + "\t" + variant['POS'] + "\t" + variant['REF'] + "\t" + variant['ALT'] +\
            "\t" + annotation['Consequence'] + "\t" + annotation['SYMBOL']
        if line not in self.lines:
            self.lines.add(line)
            self.fd.write(line + "\n")

    def close(self):
        self.fd.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('vep')
    args = parser.parse_args()
    with contextlib.closing(Annotator(args.output, args.vep)) as ann:
        ann.parse_file(args.input)


if __name__ == '__main__':
    main()
