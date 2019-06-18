#! /usr/bin/env python

import utils
from Bio import SeqIO
import os
import numpy

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rna_fam_dir_path = asset_dir_path + "/rna_families"
  sampled_rna_fam_dir_path = asset_dir_path + "/sampled_rna_families"
  num_of_samples = 10
  if not os.path.isdir(sampled_rna_fam_dir_path):
    os.mkdir(sampled_rna_fam_dir_path)
  for rna_seq_file in os.listdir(rna_fam_dir_path):
    if not rna_seq_file.endswith(".fa"):
      continue
    rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_seq_file)
    ss_file = "sss_of_" + rna_seq_file[: -3] + ".dat"
    ss_file_path = os.path.join(rna_fam_dir_path, ss_file)
    seq_recs = [rec for rec in SeqIO.parse(rna_seq_file_path, "fasta")]
    ss_recs = [rec for rec in SeqIO.parse(ss_file_path, "fasta")]
    num_of_recs = len(seq_recs)
    rna_seq_file_path = os.path.join(sampled_rna_fam_dir_path, rna_seq_file)
    ss_file_path = os.path.join(sampled_rna_fam_dir_path, ss_file)
    if num_of_recs <= num_of_samples:
      SeqIO.write(seq_recs, rna_seq_file_path, "fasta")
      SeqIO.write(ss_recs, ss_file_path, "fasta")
    else:
      sampled_indexes = numpy.random.choice(numpy.arange(num_of_recs), num_of_samples, replace = False)
      sampled_seq_recs = [rec for (i, rec) in enumerate(seq_recs) if i in sampled_indexes]
      sampled_ss_recs = [rec for (i, rec) in enumerate(ss_recs) if i in sampled_indexes]
      SeqIO.write(sampled_seq_recs, rna_seq_file_path, "fasta")
      SeqIO.write(sampled_ss_recs, ss_file_path, "fasta")

if __name__ == "__main__":
  main()
