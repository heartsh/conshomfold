#! /usr/bin/env python

import utils
from Bio import SeqIO
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
  num_of_threads = multiprocessing.cpu_count()
  phylofold_dir_path = asset_dir_path + "/phylofold"
  centroidfold_dir_path = asset_dir_path + "/centroidfold"
  contrafold_dir_path = asset_dir_path + "/contrafold"
  centroidhomfold_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_dir_path = asset_dir_path + "/turbofold"
  rnafold_dir_path = asset_dir_path + "/rnafold"
  temp_dir_path = "/tmp/run_ss_estimation_programs_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(phylofold_dir_path):
    os.mkdir(phylofold_dir_path)
  if not os.path.isdir(centroidfold_dir_path):
    os.mkdir(centroidfold_dir_path)
  if not os.path.isdir(contrafold_dir_path):
    os.mkdir(contrafold_dir_path)
  if not os.path.isdir(centroidhomfold_dir_path):
    os.mkdir(centroidhomfold_dir_path)
  if not os.path.isdir(turbofold_dir_path):
    os.mkdir(turbofold_dir_path)
  if not os.path.isdir(rnafold_dir_path):
    os.mkdir(rnafold_dir_path)
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  rna_dir_path = asset_dir_path + "/sampled_rna_families"
  gammas = [2. ** i for i in range(-7, 11)]
  phylofold_params = []
  turbofold_params = []
  turbofold_params_4_elapsed_time = []
  centroidfold_params = []
  contrafold_params = []
  centroidhomfold_params = []
  centroidfold_params_4_elapsed_time = []
  contrafold_params_4_elapsed_time = []
  centroidhomfold_params_4_elapsed_time = []
  rnafold_params = []
  phylofold_elapsed_time = 0.
  sub_thread_num = 4 if num_of_threads <= 8 else 8
  for rna_file in os.listdir(rna_dir_path):
    if not rna_file.endswith(".fa"):
      continue
    rna_file_path = os.path.join(rna_dir_path, rna_file)
    (rna_family_name, extension) = os.path.splitext(rna_file)
    phylofold_output_dir_path = os.path.join(phylofold_dir_path, "sss_of_" + rna_family_name)
    centroidfold_output_dir_path = os.path.join(centroidfold_dir_path, "sss_of_" + rna_family_name)
    contrafold_output_dir_path = os.path.join(contrafold_dir_path, "sss_of_" + rna_family_name)
    centroidhomfold_output_dir_path = os.path.join(centroidhomfold_dir_path, "sss_of_" + rna_family_name)
    turbofold_output_dir_path = os.path.join(turbofold_dir_path, "sss_of_" + rna_family_name)
    rnafold_output_file_path = os.path.join(rnafold_dir_path, "sss_of_" + rna_family_name + ".dat")
    if not os.path.isdir(phylofold_output_dir_path):
      os.mkdir(phylofold_output_dir_path)
    if not os.path.isdir(centroidfold_output_dir_path):
      os.mkdir(centroidfold_output_dir_path)
    if not os.path.isdir(contrafold_output_dir_path):
      os.mkdir(contrafold_output_dir_path)
    if not os.path.isdir(centroidhomfold_output_dir_path):
      os.mkdir(centroidhomfold_output_dir_path)
    if not os.path.isdir(turbofold_output_dir_path):
      os.mkdir(turbofold_output_dir_path)
    phylofold_command = "phylofold -u -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + phylofold_output_dir_path
    phylofold_params.insert(0, phylofold_command)
    rnafold_params.insert(0, (rna_file_path, rnafold_output_file_path))
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      output_file = "gamma=" + gamma_str + ".dat"
      centroidfold_output_file_path = os.path.join(centroidfold_output_dir_path, output_file)
      centroidfold_params.insert(0, (rna_file_path, centroidfold_output_file_path, gamma_str))
      if gamma == 1:
        centroidfold_params_4_elapsed_time.insert(0, (rna_file_path, centroidfold_output_file_path, gamma_str))
      contrafold_output_file_path = os.path.join(contrafold_output_dir_path, output_file)
      contrafold_params.insert(0, (rna_file_path, contrafold_output_file_path, gamma_str))
      if gamma == 1:
        contrafold_params_4_elapsed_time.insert(0, (rna_file_path, contrafold_output_file_path, gamma_str))
      centroidhomfold_output_file_path = os.path.join(centroidhomfold_output_dir_path, output_file)
      centroidhomfold_params.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
      if gamma == 1:
        centroidhomfold_params_4_elapsed_time.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
      turbofold_output_file_path = os.path.join(turbofold_output_dir_path, output_file)
      turbofold_params.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num))
      if gamma == 1:
        turbofold_params_4_elapsed_time.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num))
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  if False:
    begin = time.time()
    pool.map(utils.run_command, phylofold_params)
    phylofold_elapsed_time = time.time() - begin
    print("The elapsed time of the PhyloFold program for a test set = %f [s]." % phylofold_elapsed_time)
  if False:
    pool.map(run_turbofold, turbofold_params)
    begin = time.time()
    pool.map(run_turbofold, turbofold_params_4_elapsed_time)
    turbofold_elapsed_time = time.time() - begin
    pool = multiprocessing.Pool(num_of_threads)
    begin = time.time()
    pool.map(run_centroidfold, centroidfold_params_4_elapsed_time)
    centroidfold_elapsed_time = time.time() - begin
    begin = time.time()
    pool.map(run_contrafold, contrafold_params_4_elapsed_time)
    contrafold_elapsed_time = time.time() - begin
    begin = time.time()
    pool.map(run_centroidhomfold, centroidhomfold_params_4_elapsed_time)
    centroidhomfold_elapsed_time = time.time() - begin
    pool.map(run_centroidfold, centroidfold_params)
    pool.map(run_contrafold, contrafold_params)
    pool.map(run_centroidhomfold, centroidhomfold_params)
    begin = time.time()
    pool.map(run_rnafold, rnafold_params)
    rnafold_elapsed_time = time.time() - begin
    print("The elapsed time of the CentroidFold program for a test set = %f [s]." % centroidfold_elapsed_time)
    print("The elapsed time of the CONTRAfold program the for a test set = %f [s]." % contrafold_elapsed_time)
    print("The elapsed time of the CentroidHomFold program for a test set = %f [s]." % centroidhomfold_elapsed_time)
    print("The elapsed time of the TurboFold-smp program for a test set = %f [s]." % turbofold_elapsed_time)
    print("The elapsed time of the RNAfold program the for a test set = %f [s]." % rnafold_elapsed_time)
  shutil.rmtree(temp_dir_path)

def run_turbofold(turbofold_params):
  (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num) = turbofold_params
  recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
  rec_seq_len = len(recs)
  turbofold_temp_dir_path = "%s/%s_gamma=%d" % (temp_dir_path, rna_family_name, gamma)
  if not os.path.isdir(turbofold_temp_dir_path):
    os.mkdir(turbofold_temp_dir_path)
  turbofold_config_file_contents = "InSeq = {"
  for i in range(rec_seq_len):
    turbofold_config_file_contents += "%s/%d.fasta;" % (turbofold_temp_dir_path, i)
  turbofold_config_file_contents += "}\nOutCT = {"
  for i in range(rec_seq_len):
    turbofold_config_file_contents += "%s/%d.ct;" % (turbofold_temp_dir_path, i)
  turbofold_config_file_contents += "}\nIterations = 1\nMode = MEA\nMeaGamma = %f\nProcessors = %d" % (gamma, sub_thread_num)
  turbofold_config_file_path = os.path.join(turbofold_temp_dir_path, "turbofold_config.dat")
  turbofold_config_file = open(turbofold_config_file_path, "w")
  turbofold_config_file.write(turbofold_config_file_contents)
  turbofold_config_file.close()
  for (i, rec) in enumerate(recs):
    SeqIO.write([rec], open(os.path.join(turbofold_temp_dir_path, "%d.fasta" % i), "w"), "fasta")
  turbofold_command = "TurboFold-smp " + turbofold_config_file_path
  utils.run_command(turbofold_command)
  turbofold_output_file_contents = ""
  for i in range(rec_seq_len):
    ct_file_path = os.path.join(turbofold_temp_dir_path, "%d.ct" % i)
    ss_string = read_ct_file(ct_file_path)
    turbofold_output_file_contents += ">%d\n%s\n\n" % (i, ss_string)
  turbofold_output_file = open(turbofold_output_file_path, "w")
  turbofold_output_file.write(turbofold_output_file_contents)
  turbofold_output_file.close()

def run_centroidfold(centroidfold_params):
  (rna_file_path, centroidfold_output_file_path, gamma_str) = centroidfold_params
  centroidfold_command = "centroid_fold " + rna_file_path + " -g " + gamma_str
  (output, _, _) = utils.run_command(centroidfold_command)
  lines = [line.split()[0] for (i, line) in enumerate(str(output).split("\\n")) if i % 3 == 2]
  centroidfold_output_file = open(centroidfold_output_file_path, "w+")
  centroidfold_output_buf = ""
  for (i, line) in enumerate(lines):
    centroidfold_output_buf += ">%d\n%s\n\n" % (i, line)
  centroidfold_output_file.write(centroidfold_output_buf)
  centroidfold_output_file.close()

def run_contrafold(contrafold_params):
  (rna_file_path, contrafold_output_file_path, gamma_str) = contrafold_params
  contrafold_command = "centroid_fold -e CONTRAfold --mea " + rna_file_path + " -g " + gamma_str
  (output, _, _) = utils.run_command(contrafold_command)
  lines = [line.split()[0] for (i, line) in enumerate(str(output).split("\\n")) if i % 3 == 2]
  contrafold_output_file = open(contrafold_output_file_path, "w+")
  contrafold_output_buf = ""
  for (i, line) in enumerate(lines):
    contrafold_output_buf += ">%d\n%s\n\n" % (i, line)
  contrafold_output_file.write(contrafold_output_buf)
  contrafold_output_file.close()

def run_rnafold(rnafold_params):
  (rna_file_path, rnafold_output_file_path) = rnafold_params
  rnafold_command = "RNAfold -i " + rna_file_path + " > " + rnafold_output_file_path
  utils.run_command(rnafold_command)

def run_centroidhomfold(centroidhomfold_params):
  (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path) = centroidhomfold_params
  recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
  rec_seq_len = len(recs)
  centroidhomfold_output_file = open(centroidhomfold_output_file_path, "w+")
  centroidhomfold_output_buf = ""
  basename = os.path.basename(rna_file_path)
  seq_file_path = os.path.join(temp_dir_path, "seqs_4_%s_and_gamma=%s.fa" % (basename, gamma_str))
  hom_seq_file_path = os.path.join(temp_dir_path, "hom_seqs_4_%s_and_gamma=%s.fa" % (basename, gamma_str))
  for (i, rec) in enumerate(recs):
    SeqIO.write([rec], open(seq_file_path, "w"), "fasta")
    hom_recs = [rec for (j, rec) in enumerate(recs) if j != i]
    SeqIO.write(recs, open(hom_seq_file_path, "w"), "fasta")
    centroidhomfold_command = "centroid_homfold " + seq_file_path + " -H " + hom_seq_file_path + " -g " + gamma_str
    (output, _, _) = utils.run_command(centroidhomfold_command)
    centroidhomfold_output_buf += ">%d\n%s\n\n" % (i, str(output).split("\\n")[2].split()[0])
  centroidhomfold_output_file.write(centroidhomfold_output_buf)
  centroidhomfold_output_file.close()

def read_ct_file(ct_file_path):
  ct_file = open(ct_file_path, "r")
  lines = ct_file.readlines()
  seq_len = int(lines[0].split()[0])
  ss_string = ["." for i in range(seq_len)]
  num_of_lines = len(lines)
  for line in lines[1 : num_of_lines]:
    if "ENERGY" in line:
      break
    substrings = line.split()
    index_1 = int(substrings[0])
    index_2 = int(substrings[4])
    if index_2 == 0 or index_1 >= index_2:
      continue
    ss_string[index_1 - 1] = "("
    ss_string[index_2 - 1] = ")"
  return "".join(ss_string)

if __name__ == "__main__":
  main()
