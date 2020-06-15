#! /usr/bin/env python

import utils
from Bio import SeqIO
from Bio import AlignIO
import numpy
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time
import datetime
import shutil

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rfam_seed_sta_file_path = asset_dir_path + "/rfam_seed_stas.sth"
  for data_set in ["short", "long"]:
    rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set
    rna_seq_dir_path_4_micro_bench = asset_dir_path + "/compiled_rna_fams_" + data_set + "_4_micro_bench"
    ref_ss_dir_path = asset_dir_path + "/ref_sss_" + data_set
    ref_ss_dir_path_4_micro_bench = asset_dir_path + "/ref_sss_" + data_set + "_4_micro_bench"
    if not os.path.isdir(rna_seq_dir_path):
      os.mkdir(rna_seq_dir_path)
    if not os.path.isdir(rna_seq_dir_path_4_micro_bench):
      os.mkdir(rna_seq_dir_path_4_micro_bench)
    if not os.path.isdir(ref_ss_dir_path):
      os.mkdir(ref_ss_dir_path)
    if not os.path.isdir(ref_ss_dir_path_4_micro_bench):
      os.mkdir(ref_ss_dir_path_4_micro_bench)
    min_sa_len = 0 if data_set == "short" else 201
    max_sa_len = 200 if data_set == "short" else 500
    stas = [sta for sta in AlignIO.parse(rfam_seed_sta_file_path, "stockholm") if min_sa_len <= len(sta[0]) and len(sta[0]) <= max_sa_len and is_valid(sta)]
    num_of_stas = len(stas)
    print("# %s RNA families: %d" % (data_set, num_of_stas))
    sample_rate = 0.01 if data_set == "short" else 0.05
    num_of_samples = int(sample_rate * num_of_stas)
    print("# %s RNA families for micro benchmark: %d" % (data_set, num_of_samples))
    sampled_stas = numpy.random.choice(stas, num_of_samples, replace = False)
    max_seq_num = 10 if data_set == "short" else 5
    for i, sta in enumerate(stas):
      css = convert_css(sta.column_annotations["secondary_structure"])
      rna_seq_file_path = os.path.join(rna_seq_dir_path, "rna_fam_%d.fa" % i)
      ref_ss_file_path = os.path.join(ref_ss_dir_path, "rna_fam_%d.fa" % i)
      rna_seq_file = open(rna_seq_file_path, "w")
      ref_ss_file = open(ref_ss_file_path, "w")
      num_of_seqs = len(sta)
      indexes_are_sampled = True if num_of_seqs > max_seq_num else False
      indexes = numpy.random.choice([k for k in range(num_of_seqs)], max_seq_num, replace = False).tolist() if indexes_are_sampled else range(num_of_seqs)
      recs = [sta[k] for k in indexes]
      for j, rec in enumerate(recs):
        seq_with_gaps = str(rec.seq)
        recovered_ss = recover_ss(css, seq_with_gaps)
        seq = seq_with_gaps.replace("-", "")
        rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, seq))
        ref_ss_file.write(">%d(%s)\n%s\n" % (j, rec.id, recovered_ss))
    for i, sta in enumerate(sampled_stas):
      css = convert_css(sta.column_annotations["secondary_structure"])
      rna_seq_file_path = os.path.join(rna_seq_dir_path_4_micro_bench, "rna_fam_%d.fa" % i)
      ref_ss_file_path = os.path.join(ref_ss_dir_path_4_micro_bench, "rna_fam_%d.fa" % i)
      rna_seq_file = open(rna_seq_file_path, "w")
      ref_ss_file = open(ref_ss_file_path, "w")
      num_of_seqs = len(sta)
      indexes_are_sampled = True if num_of_seqs > max_seq_num else False
      indexes = numpy.random.choice([k for k in range(num_of_seqs)], max_seq_num, replace = False).tolist() if indexes_are_sampled else range(num_of_seqs)
      recs = [sta[k] for k in indexes]
      for j, rec in enumerate(recs):
        seq_with_gaps = str(rec.seq)
        recovered_ss = recover_ss(css, seq_with_gaps)
        seq = seq_with_gaps.replace("-", "")
        rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, seq))
        ref_ss_file.write(">%d(%s)\n%s\n" % (j, rec.id, recovered_ss))

def is_valid(sta):
  for row in sta:
    if any(char in str(row.seq) for char in "RYWSMKHBVDN"):
      return False
  return True

def convert_css(css):
  converted_css = ""
  for char in css:
    if char == "(" or char == "<" or char == "[" or char == "{":
      converted_css += "("
    elif char == ")" or char == ">" or char == "]" or char == "}":
      converted_css += ")"
    else:
      converted_css += "."
  return converted_css

def recover_ss(css, seq_with_gaps):
  pos_map = {}
  pos = 0
  for (i, char) in enumerate(seq_with_gaps):
    if char != "-":
      pos_map[i] = pos
      pos += 1
  recovered_ss = "." * pos
  stack = []
  for (i, char) in enumerate(css):
    if char == "(":
      stack.append(i)
    elif char == ")":
      j = stack.pop()
      if seq_with_gaps[j] == "-" or seq_with_gaps[i] == "-":
        continue
      mapped_j = pos_map[j]
      mapped_i = pos_map[i]
      recovered_ss = recovered_ss[: mapped_j] + "(" + recovered_ss[mapped_j + 1 :]
      recovered_ss = recovered_ss[: mapped_i] + ")" + recovered_ss[mapped_i + 1 :]
  return recovered_ss

if __name__ == "__main__":
  main()
