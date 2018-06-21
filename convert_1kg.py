"""
Example script for importing 1000 Genomes data.
"""
import argparse
import subprocess
import multiprocessing
import os
import shutil
import sys

import numpy as np
import msprime
import tsinfer
import attr
import cyvcf2
import tqdm

@attr.s()
class Site(object):
    position = attr.ib(None)
    alleles = attr.ib(None)
    genotypes = attr.ib(None)
    metadata = attr.ib({})

def filter_duplicates(vcf):
    """
    Returns the list of variants from this specified VCF with duplicate sites filtered
    out. If any site appears more than once, throw all variants away.
    """
    # TODO this has not been tested properly.
    row = next(vcf, None)
    bad_pos = -1
    for next_row in vcf:
        if bad_pos == -1 and next_row.POS != row.POS:
            yield row
        else:
            if bad_pos == -1:
                bad_pos = row.POS
            elif bad_pos != next_row.POS:
                bad_pos = -1
        row = next_row
    if row is not None and bad_pos != -1:
        yield row

def vcf_num_rows(vcf_path):
    """A quick way to get the number of sites in a vcf"""
    output = subprocess.check_output(["bcftools", "index", "--nrecords", vcf_path])
    return int(output)


def variants(vcf_path, show_progress=False, ancestral_states=None):
    """
    Yield a tuple of position, alleles, genotypes, metadata
    If ancestral_states is given, it should be a dictionary mapping name to state
    """
    progress = tqdm.tqdm(total=vcf_num_rows(vcf_path), desc="Read sites", disable=not show_progress)

    vcf = cyvcf2.VCF(vcf_path)

    num_diploids = len(vcf.samples)
    num_samples = 2 * num_diploids
    j = 0
    for row in filter_duplicates(vcf):
        progress.update()
        ancestral_state = None
        try:
            if ancestral_states is not None:
                aa = ancestral_states[row.ID]
            else:
                aa = row.INFO["AA"]
            # Format: for indels = AA|REF|ALT|IndelType; for SNPs = AA
            splits = aa.split("|")
            if len(splits[0]) == 1:
                base = splits[0].upper()
                if base in "ACTG":
                    ancestral_state = base
        except KeyError:
            pass
        #only use biallelic sites with data for all samples, where ancestral state is known
        if len(row.ALT)==1 and row.num_called == num_diploids and ancestral_state is not None:
            a = np.zeros(num_samples, dtype=np.uint8)
            all_alleles = set([ancestral_state])
            # Fill in a with genotypes.
            bases = np.array(row.gt_bases)
            for j in range(num_diploids):
                alleles = bases[j].split("|")
                for allele in alleles:
                    all_alleles.add(allele)
                a[2 * j] = alleles[0] != ancestral_state
                a[2 * j + 1] = alleles[1] != ancestral_state
            if len(all_alleles) == 2:
                all_alleles.remove(ancestral_state)
                alleles = [ancestral_state, all_alleles.pop()]
                metadata = {"ID": row.ID, "INFO": dict(row.INFO)}
                yield Site(
                    position=row.POS, alleles=alleles, genotypes=a, metadata=metadata)

    vcf.close()

def ancestors():
    vcf = cyvcf2.VCF(vcf_path)
    

def convert(vcf_file, output_file, max_variants=None, show_progress=False, ancestor_file=None):

    if max_variants is None:
        max_variants = 2**32  # Arbitrary, but > defined max for VCF

    sample_data = tsinfer.SampleData.initialise(path=output_file, num_flush_threads=2)
    vcf = cyvcf2.VCF(vcf_file)
    for sample in tqdm.tqdm(vcf.samples, desc="Read samples", disable=not show_progress):
        metadata = {"name": sample}
        sample_data.add_individual(ploidy=2, metadata)
    vcf.close()

    if ancestor_file is not None:
        ancestors = {}
        vcf = cyvcf2.VCF(ancestor_file)
        for site in tqdm.tqdm(
            vcf, total=vcf_num_rows(ancestor_file), desc="Read ancestors", disable=not show_progress):
            if site.ID and "AA" in site.INFO
                ancestors[site.ID] = site.INFO["AA"]
        vcf.close()
    else:
        ancestors = None

   for index, site in enumerate(variants(vcf_file, show_progress, ancestors)):
        sample_data.add_site(site.position, site.alleles, site.genotypes, site.metadata)
        if index == max_variants:
            break
    sample_data.finalise(command=sys.argv[0], parameters=sys.argv[1:])


def worker(t):
    vcf, output, max_variants = t
    print("Converting", vcf)
    convert(vcf, output, max_variants)


def main():
    parser = argparse.ArgumentParser(
        description="Script to convert VCF files into tsinfer input.")
    parser.add_argument(
        "vcf_pattern",
        help="The input VCF files file pattern, with {} where the autosome number will go. E.g. from http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/")
    parser.add_argument(
        "output_file_pattern",
        help="The tsinfer output file pattern, with {} where the autosome number will go")
    parser.add_argument(
        "-a", "--ancestral_file_pattern", default=None, 
        help="A file pattern for files containing ancestral allele states, with {} where the autosome number will go. This will override ancestral alleles from the main vcf.")
    parser.add_argument(
        "-n", "--max-variants", default=None, type=int,
        help="Keep only the first n variants")
    parser.add_argument(
        "--start", default=1, type=int, help="The first autosome")
    parser.add_argument(
        "--stop", default=22, type=int, help="The last autosome")
    parser.add_argument(
        "-P", "--processes", default=10, type=int, help="The number of worker processes")
    parser.add_argument(
        "-p", "--progress", action="store_true")

    args = parser.parse_args()
    # TODO Fix this up and make it an optional argument.
    chromosomes = list(range(args.start, args.stop + 1))

    # Build the file lists for the 22 autosomes.
    vcf_files = [args.vcf_pattern.format(j) for j in chromosomes]
    output_files = [args.output_file_pattern.format(j) for j in chromosomes]
    if  args.ancestral_file_pattern is None:
        ancestor_files = [None]
    elif "{}" in args.ancestral_file_pattern:
        ancestor_files = [args.ancestral_file_pattern.format(j) for j in chromosomes]
    else:
        ancestor_files = [args.ancestral_file_pattern]

    max_variants = [args.max_variants for _ in chromosomes]

    for vcf_file in vcf_files:
        if not os.path.exists(vcf_file):
            raise ValueError("{} does not exist".format(vcf_file))
    for ancestor_file in ancestor_files:
        if not os.path.exists(ancestor_file):
            raise ValueError("{} does not exist".format(vcf_file))

#     work = reversed(list(zip(
#         vcf_files, genetic_map_files, output_files, max_variants)))
#     with multiprocessing.Pool(args.processes) as pool:
#         pool.map(worker, work)
#     # for t in work:
#     #     worker(t)

    convert(
        vcf_files[0], output_files[0], args.max_variants, show_progress=args.progress, ancestor_files[0])


if __name__ == "__main__":
    main()
