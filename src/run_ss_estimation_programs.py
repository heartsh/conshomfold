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
  temp_dir_path = "/tmp/run_ss_estimation_programs_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  gammas = [2. ** i for i in range(-7, 11)]
  short_conshomfold_params = []
  short_conshomfold_params_4_elapsed_time = []
  bpp_conshomfold_params = []
  short_turbofold_params = []
  short_turbofold_params_4_elapsed_time = []
  short_centroidhomfold_params = []
  short_centroidhomfold_params_4_elapsed_time = []
  short_rnafold_params = []
  long_conshomfold_params = []
  long_conshomfold_params_4_elapsed_time = []
  long_turbofold_params = []
  long_turbofold_params_4_elapsed_time = []
  long_centroidhomfold_params = []
  long_centroidhomfold_params_4_elapsed_time = []
  long_rnafold_params = []
  bpp_conshomfold_dir_path = asset_dir_path + "/bpp_conshomfold"
  if not os.path.isdir(bpp_conshomfold_dir_path):
    os.mkdir(bpp_conshomfold_dir_path)
  for data_set in ["short", "long"]:
    conshomfold_dir_path = asset_dir_path + "/conshomfold_" + data_set
    centroidhomfold_dir_path = asset_dir_path + "/centroidhomfold_" + data_set
    turbofold_dir_path = asset_dir_path + "/turbofold_" + data_set
    rnafold_dir_path = asset_dir_path + "/rnafold_" + data_set
    if not os.path.isdir(conshomfold_dir_path):
      os.mkdir(conshomfold_dir_path)
    if not os.path.isdir(centroidhomfold_dir_path):
      os.mkdir(centroidhomfold_dir_path)
    if not os.path.isdir(turbofold_dir_path):
      os.mkdir(turbofold_dir_path)
    if not os.path.isdir(rnafold_dir_path):
      os.mkdir(rnafold_dir_path)
    # rna_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set
    rna_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set + "_4_micro_bench"
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
      rnafold_output_file_path = os.path.join(rnafold_dir_path, rna_family_name + ".fa")
      if not os.path.isdir(conshomfold_output_dir_path):
        os.mkdir(conshomfold_output_dir_path)
      if not os.path.isdir(centroidhomfold_output_dir_path):
        os.mkdir(centroidhomfold_output_dir_path)
      if not os.path.isdir(turbofold_output_dir_path):
        os.mkdir(turbofold_output_dir_path)
      conshomfold_command = "conshomfold -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + conshomfold_output_dir_path
      if data_set == "short":
        short_conshomfold_params.insert(0, conshomfold_command)
      else:
        long_conshomfold_params.insert(0, conshomfold_command)
      conshomfold_command = "conshomfold -b -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + conshomfold_output_dir_path
      if data_set == "short":
        short_conshomfold_params_4_elapsed_time.insert(0, conshomfold_command)
      else:
        long_conshomfold_params_4_elapsed_time.insert(0, conshomfold_command)
      if data_set == "short":
        bpp_conshomfold_command = "conshomfold -u -t " + str(sub_thread_num) + " -i " + rna_file_path + " -o " + bpp_conshomfold_output_dir_path
        bpp_conshomfold_params.insert(0, bpp_conshomfold_command)
      if data_set == "short":
        short_rnafold_params.insert(0, (rna_file_path, rnafold_output_file_path))
      else:
        long_rnafold_params.insert(0, (rna_file_path, rnafold_output_file_path))
      for gamma in gammas:
        gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
        output_file = "gamma=" + gamma_str + ".fa"
        centroidhomfold_output_file_path = os.path.join(centroidhomfold_output_dir_path, output_file)
        if data_set == "short":
          short_centroidhomfold_params.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
        else:
          long_centroidhomfold_params.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
        if gamma == 1:
          if data_set == "short":
            short_centroidhomfold_params_4_elapsed_time.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
          else:
            long_centroidhomfold_params_4_elapsed_time.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
        turbofold_output_file_path = os.path.join(turbofold_output_dir_path, output_file)
        if data_set == "short":
          short_turbofold_params.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num))
        else:
          long_turbofold_params.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num))
        if gamma == 1:
          if data_set == "short":
            short_turbofold_params_4_elapsed_time.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num))
          else:
            long_turbofold_params_4_elapsed_time.insert(0, (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name, sub_thread_num))
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  pool.map(utils.run_command, short_conshomfold_params)
  pool.map(utils.run_command, long_conshomfold_params)
  begin = time.time()
  pool.map(utils.run_command, short_conshomfold_params_4_elapsed_time)
  short_conshomfold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(utils.run_command, long_conshomfold_params_4_elapsed_time)
  long_conshomfold_elapsed_time = time.time() - begin
  pool.map(utils.run_command, bpp_conshomfold_params)
  pool.map(run_turbofold, short_turbofold_params)
  pool.map(run_turbofold, long_turbofold_params)
  begin = time.time()
  pool.map(run_turbofold, short_turbofold_params_4_elapsed_time)
  short_turbofold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_turbofold, long_turbofold_params_4_elapsed_time)
  long_turbofold_elapsed_time = time.time() - begin
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_centroidhomfold, short_centroidhomfold_params)
  pool.map(run_centroidhomfold, long_centroidhomfold_params)
  begin = time.time()
  pool.map(run_centroidhomfold, short_centroidhomfold_params_4_elapsed_time)
  short_centroidhomfold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_centroidhomfold, long_centroidhomfold_params_4_elapsed_time)
  long_centroidhomfold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_rnafold, short_rnafold_params)
  short_rnafold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_rnafold, long_rnafold_params)
  long_rnafold_elapsed_time = time.time() - begin
  print("The elapsed time of the ConsHomfold program for test set \"short\" = %f [s]." % short_conshomfold_elapsed_time)
  print("The elapsed time of the CentroidHomfold program for test set \"short\" = %f [s]." % short_centroidhomfold_elapsed_time)
  print("The elapsed time of the TurboFold-smp program for test set \"short\" = %f [s]." % short_turbofold_elapsed_time)
  print("The elapsed time of the RNAfold program for test set \"short\" = %f [s]." % short_rnafold_elapsed_time)
  print("The elapsed time of the ConsHomfold program for test set \"long\" = %f [s]." % long_conshomfold_elapsed_time)
  print("The elapsed time of the CentroidHomfold program for test set \"long\" = %f [s]." % long_centroidhomfold_elapsed_time)
  print("The elapsed time of the TurboFold-smp program for test set \"long\" = %f [s]." % long_turbofold_elapsed_time)
  print("The elapsed time of the RNAfold program for test set \"long\" = %f [s]." % long_rnafold_elapsed_time)
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
  turbofold_config_file_contents += "}\nIterations = 3\nMode = MEA\nMeaGamma = %f\nProcessors = %d" % (gamma, sub_thread_num)
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

def run_rnafold(rnafold_params):
  (rna_file_path, rnafold_output_file_path) = rnafold_params
  rnafold_command = "RNAfold --nops -i " + rna_file_path + " > " + rnafold_output_file_path
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
