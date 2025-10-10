import argparse
import re
import sys
from enum import Enum
import os
import contextlib

class VepParserError(Exception):
    """Base exception for VEP parser errors"""
    pass


class VepFormatError(VepParserError):
    """Exception raised when VEP format is incorrect"""
    pass


class VepParseError(VepParserError):
    """Exception raised during line parsing"""
    def __init__(self, message, line_number=None, line_content=None):
        self.line_number = line_number
        self.line_content = line_content
        if line_number is not None:
            message = f"Line {line_number}: {message}"
        if line_content is not None:
            message += f"\n  Content: {line_content[:100]}"
        super().__init__(message)


class VepParser:
    vcf_format = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO"

    def __init__(self, format_line = None):
        self.init = False
        self.header = ""
        self.line_number = 0
        if format_line is not None:
            print(format_line)
            self.fields = format_line.split("|")
            self.initialize()
        self.vcf_fields = self.vcf_format.split("|")

    def initialize(self):
        """Virtual method to be overridden by subclasses if needed."""
        self.init = True

    def parse_header(self, header):
        expr = re.compile("##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: ([^\"]+)\">")
        match = expr.search(header)
        if match is None:
            raise VepFormatError("Could not find the CSQ field in the VCF header. Make sure the file is annotated with VEP.")
        self.fields = match.group(1).split("|")

    def parse_file(self, fd):
        for line in fd:
            self.line_number += 1
            if line.startswith("#"):
                self.header += line
                continue
            else:
                if not self.init:
                    self.parse_header(self.header)
                    self.initialize()
                try:
                    self.parse_line(line.rstrip())
                except VepParserError:
                    raise
                except Exception as e:
                    raise VepParseError(
                        f"Unexpected error during parsing: {str(e)}",
                        line_number=self.line_number,
                        line_content=line.rstrip()
                    ) from e

    def extract_vep(self, info):
        csq = None
        fields = info.split(";")
        for field in fields:
            if field.startswith("CSQ="):
                csq = field[4:]
                break
        if csq is None:
            raise VepParseError(
                "Could not find CSQ field in the INFO column",
                line_number=self.line_number
            )
        if not csq or csq == "":
            # Empty CSQ field - no annotations
            return []
        transcripts = csq.split(",")
        result = []
        for t in transcripts:
            tokens = t.split("|")
            if len(tokens) != len(self.fields):
                raise VepParseError(
                    f"Number of fields in CSQ ({len(tokens)}) does not match expected format ({len(self.fields)})",
                    line_number=self.line_number,
                    line_content=info
                )
            result.append(dict(zip(self.fields, tokens)))
        return result

    def parse_line(self, line):
        vcf_fields = line.split("\t")
        if len(vcf_fields) < 8:
            raise VepParseError(
                f"Invalid VCF format: expected at least 8 columns, got {len(vcf_fields)}",
                line_number=self.line_number,
                line_content=line
            )
        
        variant = vcf_fields[0:7]
        variant_dict = dict((self.vcf_fields[i], variant[i]) for i in range(len(variant)))
        info = vcf_fields[7]
        if info == ".":
            return
        vep = self.extract_vep(info)
        for annotation in vep:
            try:
                gene = annotation.get('SYMBOL', '')
                if gene == '':
                    continue
                if annotation.get('CANONICAL', '') != 'YES':
                    continue
                biotype = annotation.get('BIOTYPE', '')
                if biotype != 'protein_coding':
                    continue
                try:
                    annotation['am_pathogenicity'] = float(annotation['am_pathogenicity'])
                except (ValueError, KeyError):
                    annotation['am_pathogenicity'] = None
                cs = annotation.get('Consequence', '').split("&")
                if len(cs) == 0:
                    continue
                for cons in cs:
                    annotation['Consequence'] = cons
                    self.process_annotation(annotation, variant_dict)
            except KeyError as e:
                raise VepParseError(
                    f"Missing expected field in VEP annotation: {str(e)}",
                    line_number=self.line_number,
                    line_content=line
                ) from e

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
    def __init__(self, output, vep_format, am):
        self.am = am
        self.lines = set()
        super().__init__(vep_format)

        header = "chr\tpos\tref\talt\treason\tgene"
        if am:
            header += "\tam_pathogenicity\tam_class"
        self.fd = open(output, "w")
        self.fd.write(header + "\n")

    def initialize(self):
        super().initialize()
        if self.am and 'am_pathogenicity' not in self.fields:
            raise VepFormatError(
                "am_pathogenicity field not found in the VEP format. "
                "Make sure VEP was run with AlphaMissense plugin or remove --am flag."
            )

    def process_annotation(self, annotation, variant):
        # vcf_format = "CHROM|POS|ID|REF|ALT|QUAL|FILTER|INFO"
        line = variant['CHROM'] + "\t" + variant['POS'] + "\t" + variant['REF'] + "\t" + variant['ALT'] +\
            "\t" + annotation['Consequence'] + "\t" + annotation['SYMBOL']
        if self.am:
            path = annotation['am_pathogenicity']
            am_class = ""
            if path is not None:
                am_class = "ambiguous"
                if path < 0.34:
                    am_class = "likely_benign"
                if path > 0.564:
                    am_class = "likely_pathogenic"
            else:
                am_class = "NA"
                path = "NA"
            line += "\t" + str(path)
            line += "\t" + am_class
        if line not in self.lines:
            self.lines.add(line)
            self.fd.write(line + "\n")

    def close(self):
        self.fd.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--am', '--alpha-missense', action = 'store_true')
    parser.add_argument('--vep-format', dest = 'vep')
    parser.add_argument('--output', default='vep.tsv')
    parser.add_argument('input')
    args = parser.parse_args()

    try:
        with contextlib.closing(Annotator(args.output, args.vep, args.am)) as ann:
            if args.input == '-':
                ann.parse_file(sys.stdin)
            else:
                with open(args.input) as fd:
                    ann.parse_file(fd)
    except VepFormatError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except VepParseError as e:
        print(f"Parse error: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    except IOError as e:
        print(f"Error: IO error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
