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
from os import path

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  temp_dir_path = "/tmp/run_ss_estimation_programs_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  gammas = [2. ** i for i in range(-7, 11)]
  conshomfold_params = []
  conshomfold_params_4_elapsed_time = []
  bpp_conshomfold_params = []
  bpp_conshomfold_params_4_elapsed_time = []
  turbofold_params = []
  turbofold_params_4_elapsed_time = []
  centroidhomfold_params = []
  centroidhomfold_params_4_elapsed_time = []
  contrafold_params = []
  contrafold_params_4_elapsed_time = []
  centroidfold_params = []
  centroidfold_params_4_elapsed_time = []
  rnafold_params = []
  conshomfold_dir_path = asset_dir_path + "/conshomfold"
  centroidhomfold_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_dir_path = asset_dir_path + "/turbofold"
  contrafold_dir_path = asset_dir_path + "/contrafold"
  centroidfold_dir_path = asset_dir_path + "/centroidfold"
  rnafold_dir_path = asset_dir_path + "/rnafold"
  bpp_conshomfold_dir_path = asset_dir_path + "/bpp_conshomfold"
  if not os.path.isdir(conshomfold_dir_path):
    os.mkdir(conshomfold_dir_path)
  if not os.path.isdir(centroidhomfold_dir_path):
    os.mkdir(centroidhomfold_dir_path)
  if not os.path.isdir(turbofold_dir_path):
    os.mkdir(turbofold_dir_path)
  if not os.path.isdir(contrafold_dir_path):
    os.mkdir(contrafold_dir_path)
  if not os.path.isdir(centroidfold_dir_path):
    os.mkdir(centroidfold_dir_path)
  if not os.path.isdir(rnafold_dir_path):
    os.mkdir(rnafold_dir_path)
  if not os.path.isdir(bpp_conshomfold_dir_path):
    os.mkdir(bpp_conshomfold_dir_path)
  # rna_dir_path = asset_dir_path + "/compiled_rna_fams"
  rna_dir_path = asset_dir_path + "/compiled_rna_fams_4_micro_bench"
  sub_thread_num = 4 if num_of_threads <= 8 else 8
  for rna_file in os.listdir(rna_dir_path):
    if not rna_file.endswith(".fa"):
      continue
    rna_file_path = os.path.join(rna_dir_path, rna_file)
    (rna_family_name, extension) = os.path.splitext(rna_file)
    conshomfold_output_dir_path = os.path.join(conshomfold_dir_path, rna_family_name)
    bpp_conshomfold_output_dir_path = os.path.join(bpp_conshomfold_dir_path, rna_family_name)
    centroidhomfold_output_dir_path = os.path.join(centroidhomfold_dir_path, rna_family_name)
    turbofold_output_dir_path = os.path.join(turbofold_dir_path, rna_family_name)
    contrafold_output_dir_path = os.path.join(contrafold_dir_path, rna_family_name)
    centroidfold_output_dir_path = os.path.join(centroidfold_dir_path, rna_family_name)
    rnafold_output_file_path = os.path.join(rnafold_dir_path, rna_family_name + ".fa")
    if not os.path.isdir(conshomfold_output_dir_path):
      os.mkdir(conshomfold_output_dir_path)
    if not os.path.isdir(centroidhomfold_output_dir_path):
      os.mkdir(centroidhomfold_output_dir_path)
    if not os.path.isdir(turbofold_output_dir_path):
      os.mkdir(turbofold_output_dir_path)
    if not os.path.isdir(contrafold_output_dir_path):
      os.mkdir(contrafold_output_dir_path)
    if not os.path.isdir(centroidfold_output_dir_path):
      os.mkdir(centroidfold_output_dir_path)
    conshomfold_command = "conshomfold -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + conshomfold_output_dir_path
    conshomfold_params.insert(0, conshomfold_command)
    conshomfold_command = "conshomfold -b -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + conshomfold_output_dir_path
    conshomfold_params_4_elapsed_time.insert(0, conshomfold_command)
    bpp_conshomfold_command = "conshomfold -u -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + bpp_conshomfold_output_dir_path
    bpp_conshomfold_params.insert(0, bpp_conshomfold_command)
    bpp_conshomfold_command = "conshomfold -b -u -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + bpp_conshomfold_output_dir_path
    bpp_conshomfold_params_4_elapsed_time.insert(0, bpp_conshomfold_command)
    rnafold_params.insert(0, (rna_file_path, rnafold_output_file_path))
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      output_file = "gamma=" + gamma_str + ".fa"
      centroidhomfold_output_file_path = os.path.join(centroidhomfold_output_dir_path, output_file)
      centroidhomfold_params.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
      if gamma == 1:
        centroidhomfold_params_4_elapsed_time.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
      turbofold_output_file_path = os.path.join(turbofold_output_dir_path, output_file)
      turbofold_params.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name))
      if gamma == 1:
        turbofold_params_4_elapsed_time.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name))
      contrafold_output_file_path = os.path.join(contrafold_output_dir_path, output_file)
      contrafold_params.insert(0, (rna_file_path, contrafold_output_file_path, gamma_str, temp_dir_path))
      if gamma == 1:
        contrafold_params_4_elapsed_time.insert(0, (rna_file_path, contrafold_output_file_path, gamma_str, temp_dir_path))
      centroidfold_output_file_path = os.path.join(centroidfold_output_dir_path, output_file)
      centroidfold_params.insert(0, (rna_file_path, centroidfold_output_file_path, gamma_str))
      if gamma == 1:
        centroidfold_params_4_elapsed_time.insert(0, (rna_file_path, centroidfold_output_file_path, gamma_str))
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  pool.map(utils.run_command, conshomfold_params)
  begin = time.time()
  pool.map(utils.run_command, conshomfold_params_4_elapsed_time)
  conshomfold_elapsed_time = time.time() - begin
  pool.map(utils.run_command, bpp_conshomfold_params)
  begin = time.time()
  pool.map(utils.run_command, bpp_conshomfold_params_4_elapsed_time)
  bpp_conshomfold_elapsed_time = time.time() - begin
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_turbofold, turbofold_params)
  begin = time.time()
  pool.map(run_turbofold, turbofold_params_4_elapsed_time)
  turbofold_elapsed_time = time.time() - begin
  pool.map(run_centroidhomfold, centroidhomfold_params)
  begin = time.time()
  pool.map(run_centroidhomfold, centroidhomfold_params_4_elapsed_time)
  centroidhomfold_elapsed_time = time.time() - begin
  pool.map(run_contrafold, contrafold_params)
  begin = time.time()
  pool.map(run_contrafold, contrafold_params_4_elapsed_time)
  contrafold_elapsed_time = time.time() - begin
  pool.map(run_centroidfold, centroidfold_params)
  begin = time.time()
  pool.map(run_centroidfold, centroidfold_params_4_elapsed_time)
  centroidfold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_rnafold, rnafold_params)
  rnafold_elapsed_time = time.time() - begin
  print("The elapsed time of ConsHomfold (Turner) = %f [s]." % conshomfold_elapsed_time)
  print("The elapsed time of ConsHomfold (Posterior) = %f [s]." % bpp_conshomfold_elapsed_time)
  print("The elapsed time of CentroidHomfold = %f [s]." % centroidhomfold_elapsed_time)
  print("The elapsed time of TurboFold = %f [s]." % turbofold_elapsed_time)
  print("The elapsed time of CONTRAfold = %f [s]." % contrafold_elapsed_time)
  print("The elapsed time of CentroidFold = %f [s]." % centroidfold_elapsed_time)
  print("The elapsed time of RNAfold = %f [s]." % rnafold_elapsed_time)
  shutil.rmtree(temp_dir_path)

def run_turbofold(turbofold_params):
  (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name) = turbofold_params
  recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
  rec_lens = [len(rec) for rec in recs]
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
  turbofold_config_file_contents += "}\nIterations = 3\nMode = MEA\nMeaGamma = %f" % gamma
  turbofold_config_file_path = os.path.join(turbofold_temp_dir_path, "turbofold_config.dat")
  turbofold_config_file = open(turbofold_config_file_path, "w")
  turbofold_config_file.write(turbofold_config_file_contents)
  turbofold_config_file.close()
  for (i, rec) in enumerate(recs):
    SeqIO.write([rec], open(os.path.join(turbofold_temp_dir_path, "%d.fasta" % i), "w"), "fasta")
  turbofold_command = "TurboFold " + turbofold_config_file_path
  utils.run_command(turbofold_command)
  turbofold_output_file_contents = ""
  all_files_exist = True
  for i in range(rec_seq_len):
    ct_file_path = os.path.join(turbofold_temp_dir_path, "%d.ct" % i)
    if path.exists(ct_file_path):
      ss_string = read_ct_file(ct_file_path)
      turbofold_output_file_contents += ">%d\n%s\n\n" % (i, ss_string)
    else:
      all_files_exist = False
      turbofold_output_file_contents = ""
      break
  if not all_files_exist:
    print("Some output files are empty. TurboFold is retried with # iterations = 1.")
    turbofold_config_file_contents = "InSeq = {"
    for i in range(rec_seq_len):
      turbofold_config_file_contents += "%s/%d.fasta;" % (turbofold_temp_dir_path, i)
    turbofold_config_file_contents += "}\nOutCT = {"
    for i in range(rec_seq_len):
      turbofold_config_file_contents += "%s/%d.ct;" % (turbofold_temp_dir_path, i)
    turbofold_config_file_contents += "}\nIterations = 1\nMode = MEA\nMeaGamma = %f" % gamma
    turbofold_config_file_path = os.path.join(turbofold_temp_dir_path, "turbofold_config.dat")
    turbofold_config_file = open(turbofold_config_file_path, "w")
    turbofold_config_file.write(turbofold_config_file_contents)
    turbofold_config_file.close()
    utils.run_command(turbofold_command)
    for i in range(rec_seq_len):
      ct_file_path = os.path.join(turbofold_temp_dir_path, "%d.ct" % i)
      ss_string = read_ct_file(ct_file_path)
      turbofold_output_file_contents += ">%d\n%s\n\n" % (i, ss_string)
  turbofold_output_file = open(turbofold_output_file_path, "w")
  turbofold_output_file.write(turbofold_output_file_contents)
  turbofold_output_file.close()

def run_rnafold(rnafold_params):
  (rna_file_path, rnafold_output_file_path) = rnafold_params
  rnafold_command = "RNAfold --noPS -i " + rna_file_path
  (output, _, _) = utils.run_command(rnafold_command)
  lines = [str(line) for line in str(output).strip().split("\\n") if line.startswith(".") or line.startswith("(")]
  buf = ""
  for i, line in enumerate(lines):
    buf += ">%d\n%s\n" % (i, line.split()[0])
  rnafold_output_file = open(rnafold_output_file_path, "w")
  rnafold_output_file.write(buf)
  rnafold_output_file.close()

def run_contrafold(contrafold_params):
  (rna_file_path, contrafold_output_file_path, gamma_str, temp_dir_path) = contrafold_params
  contrafold_output_file = open(contrafold_output_file_path, "w+")
  contrafold_output_buf = ""
  basename = os.path.basename(rna_file_path)
  seq_file_path = os.path.join(temp_dir_path, "seq_4_%s_and_gamma=%s.fa" % (basename, gamma_str))
  recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
  for (i, rec) in enumerate(recs):
    SeqIO.write([rec], open(seq_file_path, "w"), "fasta")
    contrafold_command = "contrafold predict " + seq_file_path + " --gamma " + gamma_str
    (output, _, _) = utils.run_command(contrafold_command)
    contrafold_output_buf += ">%d\n%s\n\n" % (i, str(output).strip().split("\\n")[3])
  contrafold_output_file.write(contrafold_output_buf)
  contrafold_output_file.close()

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
