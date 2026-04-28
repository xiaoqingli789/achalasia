#!/usr/bin/env python3
"""
Process genotypes for families to identify recessive and de novo patterns.
Inputs:
  - family_file: tab-separated: child, father, mother (one per line)
  - genotype_file: first row = sample names (including variant ID column), 
                   each row = variant ID + genotypes for all samples
  - output_file: results with counts and patterns
"""

import argparse

def parse_families(family_file):
    """Read family structure: child -> (father, mother)"""
    families = {}
    with open(family_file, 'r') as f:
        for line in f:
            child, father, mother = line.strip().split('\t')
            families[child] = (father, mother)
    return families

def parse_genotypes(genotype_file):
    """Parse genotype file into dict: variant -> {sample: genotype}"""
    genotypes = {}
    with open(genotype_file, 'r') as f:
        header = f.readline().strip('\n').split('\t')
        # First column is variant ID, remaining are sample names
        sample_names = header[1:]
        for line in f:
            fields = line.strip().split('\t')
            variant = fields[0]
            sample_genos = fields[1:]
            # Ensure same number of samples
            if len(sample_genos) != len(sample_names):
                continue  # or raise error, but skip malformed lines
            genotypes[variant] = dict(zip(sample_names, sample_genos))
    return genotypes

def process_variant(variant, sample_genos, families):
    """Analyze one variant across families.
    Returns: tuple (rec_list, child_hom_count, father_hom_count, mother_hom_count,
                   child_list, family_gt_list)
    """
    rec_list = []          # children with recessive pattern (child hom, parents het)
    child_list = []        # all children carrying the variant (hom or het)
    family_gt_list = []    # genotype trios for each family (child|father|mother)
    child_hom = 0
    father_hom = 0
    mother_hom = 0
    AC = 0                 # allele count (not used in output but kept for logic)

    for child, (father, mother) in families.items():
        family_gt = []
        # Check if all three members have genotype data
        if child in sample_genos and father in sample_genos and mother in sample_genos:
            child_gt = sample_genos[child]
            father_gt = sample_genos[father]
            mother_gt = sample_genos[mother]

            # Count child homozygotes
            if child_gt in ['1/1']:
                child_hom += 1
                child_list.append(child)
                family_gt.extend([child_gt, father_gt, mother_gt])
                AC += 2

            # Count father homozygotes
            if father_gt in ['1/1']:
                father_hom += 1
                AC += 2

            # Count mother homozygotes
            if mother_gt in ['1/1']:
                mother_hom += 1
                AC += 2

            # Count child heterozygotes
            if child_gt in ['1/0', '0/1']:
                child_list.append(child)
                family_gt.extend([child_gt, father_gt, mother_gt])
                AC += 1

            # Count parent heterozygotes (for AC)
            if father_gt in ['1/0', '0/1']:
                AC += 1
            if mother_gt in ['1/0', '0/1']:
                AC += 1

            # Recessive pattern: child homozygous, both parents heterozygous
            if (child_gt in ['1/1'] and 
                father_gt in ['1/0', '0/1'] and 
                mother_gt in ['1/0', '0/1']):
                rec_list.append(child)

        # Record this family's trio (may be empty if any member missing)
        family_gt_list.append('|'.join(family_gt))

    return rec_list, child_hom, father_hom, mother_hom, child_list, family_gt_list

def main():
    parser = argparse.ArgumentParser(description='Process genotypes for recessive/de novo patterns')
    parser.add_argument('input_file', help='Family file (child,father,mother per line)')
    parser.add_argument('genotype_file', help='Genotype table (variant ID + sample genotypes)')
    parser.add_argument('output_file', help='Output file')
    args = parser.parse_args()

    families = parse_families(args.input_file)
    genotypes = parse_genotypes(args.genotype_file)

    with open(args.output_file, 'w') as out:
        for variant, sample_genos in genotypes.items():
            rec_list, child_hom, father_hom, mother_hom, child_list, family_gt_list = \
                process_variant(variant, sample_genos, families)

            # Format output fields with leading semicolons for empty cases
            rec_str = ';' + ';'.join(rec_list)
            child_list_str = ';' + ';'.join(child_list)
            family_gt_str = ';' + ';'.join(family_gt_list)

            out.write(f"{variant}\t{rec_str}\t{child_hom}\t{father_hom}\t{mother_hom}\t{child_list_str}\t{family_gt_str}\n")

if __name__ == '__main__':
    main()