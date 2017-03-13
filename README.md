# ExACTools
Tools for extracting information from ExAC VCF files when you need 'one variant per line' type output; some code sourced from a parser by Konrad Karczewski (@konradjk). Re-used code is attributed inline; any errors are my own.

Notes
-----

The ExAC/gnomAD VCFs store annotations in two places; at the top level of the VCF (parsed in the usual way), and in the CSQ field (LOFTEE annotations). The code supplied should give a general idea as to how to access these (look at the annotations and info_fields objects in particular). As of gnomAD v0.1r2 a significant number of fields are left empty, including what looks like allocated space for ExAC frequencies and data.

A brief summary of command line options follows:

Input files (must specify at least one of these flags):
-------------------------------------------------------

```--vcf <path>``` : a single VCF file (e.g. the exomes VCF); may be gzipped.

```--vcfdir <path>``` : a directory containing multiple VCFs (e.g. the genomes directory where one chromosome = one VCF); files may be gzipped, code (should) also work when there is a mixture of gzipped and uncompressed text.

Flags:
------

```--gnomad``` : Specifies that the file is from gnomAD. This adds an additional population (Ashkenazi Jewish, ASJ) over ExAC. Running this code on gnomAD files without the flag should still work, but the ASJ population will NOT be processed.

```--allvars``` : Default behaviour is to process only high confidence LOF calls from LOFTEE. This removes this check (subject to --noncanonical)

```--noncanonical``` : Default behaviour is to process only canonical transcripts; this will remove this check and process all variants (subject to --allvars)

Default values are false for all of these flags; the code will (without any flags) process HC LOF calls in canonical transcripts in ExAC populations.
