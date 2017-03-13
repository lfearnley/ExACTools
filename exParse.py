import argparse
import gzip
import re
import sys
import os
from decimal import *


__author__ = 'lfearnley'
# The initial version of this code was based on the ExAC parsing example provided by Konrad Karczewski (konradjk);
# attribution for re-used code is made in-line.


ac_flags = ["AC_AFR", "AC_AMR", "AC_SAS", "AC_EAS", "AC_FIN", "AC_NFE"]
an_flags = ["AN_AFR", "AN_AMR", "AN_SAS", "AN_EAS", "AN_FIN", "AN_NFE"]


def main(args):
    files = list()
    # handle files in a vcf directory...
    if args.vcfdir:
        for filename in os.listdir(args.vcfdir):
            if filename.__contains__("vcf"):
                files.append(os.path.join(args.vcfdir, filename))
    elif args.vcf:
        files.append(args.vcf)
    outputfile = open("variantdetails.txt", 'w')
    outputfile.write(
        "HGVSC\tENSGID\tCHROM\tPOS\tRSID\tAllele\tSIFT\tPolyPhen\tConsequence\tAFR\tAMR\tSAS\tEAS\tFIN\tNFE"),
    # gnomAD files add Ashkenazi Jewish populations; as a temporary fix, append the ASJ pop code to the lists.
    if args.gnomad:
        ac_flags.append("AC_ASJ")
        an_flags.append("AN_ASJ")
        outputfile.write("\tASJ"),
    outputfile.write("\n"),
    for currentfile in files:
        if currentfile.endswith('.gz'):
            f = gzip.open(currentfile)
        else:
            f = open(currentfile)
        processfile(f, outputfile, args.allvars, args.noncanonical)
        f.close()
    outputfile.close()


def processfile(f, outputfile, allvars, noncanonical):
    # Begin parser sourced from @konradjk:
    vep_field_names = None
    header = None
    for line in f:
        line = line.strip()
        # Reading header lines to get VEP and individual arrays
        if line.startswith('#'):
            line = line.lstrip('#')
            if line.find('ID=CSQ') > -1:
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
            if line.startswith('CHROM'):
                header = line.split()
                header = dict(zip(header, range(len(header))))
            continue
        if vep_field_names is None:
            print >> sys.stderr, "VCF file does not have a VEP header line. Exiting."
            sys.exit(1)
        if header is None:
            print >> sys.stderr, "VCF file does not have a header line (CHROM POS etc.). Exiting."
            sys.exit(1)
        # Pull out annotation info from INFO and ALT fields
        fields = line.split('\t')
        info_field = dict(
            [(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[header['INFO']])])
        # Only reading lines with an annotation after this point
        if 'CSQ' not in info_field:
            continue
        annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',') if
                       len(vep_field_names) == len(x.split('|'))]
        # End parser from @konradjk
        allele_split = fields[header['ALT']].split(',')
        coordstring = fields[header['CHROM']] + '\t' + fields[header['POS']] + '\t' + fields[header['ID']]
        processannotation(annotations, outputfile, coordstring, allele_split, info_field, allvars, noncanonical)


def processannotation(annotations, outputfile, coordstring, allele_split, info_field, allvars, noncanonical):
    for annotation in annotations:
        # User may choose to process only canonical transcripts:
        if not noncanonical and not annotation['CANONICAL']:
            continue
        # User may choose to process only HC LoF calls:
        if not allvars and not annotation['LoF'] == 'HC':
            continue
        # Write variant information:
        outputfile.write(annotation['HGVSc']),
        outputfile.write("\t"),
        outputfile.write(annotation['Gene']),
        outputfile.write('\t'),
        outputfile.write(coordstring)
        outputfile.write('\t')
        outputfile.write(annotation['Allele']),
        outputfile.write('\t'),
        outputfile.write(annotation['SIFT']),
        outputfile.write('\t'),
        outputfile.write(annotation['PolyPhen']),
        outputfile.write('\t'),
        outputfile.write(annotation['Consequence']),
        # Get allele index (where multiple alleles present):
        allele_idx = 0
        if len(allele_split) > 1:
            if not annotation['Allele'] in allele_split:
                for idx, val in enumerate(allele_split):
                    if annotation['Allele'] in val:
                        allele_idx = idx
            else:
                allele_idx = allele_split.index(annotation['Allele'])
        # Iterate over all populations, output the frequencies.
        for i in range(0, len(ac_flags)):
            # gnomAD r2 v0.1 has no SAS genomes; write a zero if the flag is missing.
            if ac_flags[i] not in info_field:
                outputfile.write('\t')
                outputfile.write("0")
            else:
                AF = 0
                AC_SPLIT = info_field[ac_flags[i]].split(',')
                AN_SPLIT = info_field[an_flags[i]].split(',')
                num = Decimal(AC_SPLIT[allele_idx])
                if len(AN_SPLIT) > 1:
                    div = Decimal(AN_SPLIT[allele_idx])
                else:
                    div = Decimal(AN_SPLIT[0])
                if div == 0:
                    AF = 0
                elif num > 0:
                    AF = num / div
                outputfile.write('\t'),
                outputfile.write(str(AF)),
        outputfile.write('\n'),


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # Modified after @konradjk:
    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file (from VEP+LoF); may be gzipped', required=False)
    # Added directory-based input (for gnomAD genomes, which are presented per-chromosome:
    parser.add_argument('--vcfdir', help='Input directory containing VCF files (from VEP+LoF); VCFs may be gzipped',
                        required=False)
    # Add flags; distinguish between gnomAD/ExAC VCFs, allow selection of all variants (not just LoF/canonical).
    parser.add_argument('--gnomad', dest='gnomad', action='store_true', help='Use this flag to process gnomAD data')
    parser.add_argument('--allvars', dest='allvars', action='store_true',
                        help='Use this flag to process all variants, not only HC LoFs.')
    parser.add_argument('--noncanonical', dest='noncanonical', action='store_true',
                        help='Use this flag to process noncanonical transcripts.')
    parser.set_defaults(allvars=False)
    parser.set_defaults(noncanonical=False)
    parser.set_defaults(gnomad=False)
    args = parser.parse_args()
    main(args)
