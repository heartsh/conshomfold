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
  strap_dir_path = asset_dir_path + "/strap"
  neofold_dir_path = asset_dir_path + "/neofold"
  parasor_dir_path = asset_dir_path + "/parasor"
  centroidhomfold_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_dir_path = asset_dir_path + "/turbofold"
  temp_dir_path = "/tmp/run_ss_estimation_programs_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(strap_dir_path):
    os.mkdir(strap_dir_path)
  if not os.path.isdir(neofold_dir_path):
    os.mkdir(neofold_dir_path)
  if not os.path.isdir(parasor_dir_path):
    os.mkdir(parasor_dir_path)
  if not os.path.isdir(centroidhomfold_dir_path):
    os.mkdir(centroidhomfold_dir_path)
  if not os.path.isdir(turbofold_dir_path):
    os.mkdir(turbofold_dir_path)
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  rna_dir_path = asset_dir_path + "/sampled_rna_families"
  bpp_mat_file = "bpp_mats_on_sta.dat"
  gammas = [2. ** i for i in range(-7, 11)]
  parasor_params = []
  centroidhomfold_params = []
  strap_and_neofold_elapsed_time = 0.
  turbofold_elapsed_time = 0.
  for rna_file in os.listdir(rna_dir_path):
    if not rna_file.endswith(".fa"):
      continue
    rna_file_path = os.path.join(rna_dir_path, rna_file)
    (rna_familiy_name, extension) = os.path.splitext(rna_file)
    strap_output_dir_path = os.path.join(strap_dir_path, rna_familiy_name)
    strap_command = "strap --num_of_times_of_improvements_of_struct_align_prob_mat_quadruples 0 -i " + rna_file_path + " -o " + strap_output_dir_path
    begin = time.time()
    utils.run_command(strap_command)
    elapsed_time = time.time() - begin
    strap_and_neofold_elapsed_time += elapsed_time
    bpp_mat_file_path = os.path.join(strap_output_dir_path, bpp_mat_file)
    neofold_output_dir_path = os.path.join(neofold_dir_path, "sss_of_" + rna_familiy_name)
    parasor_output_dir_path = os.path.join(parasor_dir_path, "sss_of_" + rna_familiy_name)
    centroidhomfold_output_dir_path = os.path.join(centroidhomfold_dir_path, "sss_of_" + rna_familiy_name)
    turbofold_output_dir_path = os.path.join(turbofold_dir_path, "sss_of_" + rna_familiy_name)
    if not os.path.isdir(neofold_output_dir_path):
      os.mkdir(neofold_output_dir_path)
    if not os.path.isdir(parasor_output_dir_path):
      os.mkdir(parasor_output_dir_path)
    if not os.path.isdir(centroidhomfold_output_dir_path):
      os.mkdir(centroidhomfold_output_dir_path)
    if not os.path.isdir(turbofold_output_dir_path):
      os.mkdir(turbofold_output_dir_path)
    for gamma in gammas:
      gamma_str = str(gamma)
      output_file = "gamma=" + gamma_str + ".dat"
      neofold_output_file_path = os.path.join(neofold_output_dir_path, output_file)
      neofold_command = "neofold -f " + rna_file_path + " -p " + bpp_mat_file_path + " -o " + neofold_output_file_path + " --gamma " + gamma_str
      begin = time.time()
      utils.run_command(neofold_command)
      elapsed_time = time.time() - begin
      strap_and_neofold_elapsed_time += elapsed_time
      parasor_output_file_path = os.path.join(parasor_output_dir_path, output_file)
      parasor_params.insert(0, (rna_file_path, gamma, parasor_output_file_path))
      centroidhomfold_output_file_path = os.path.join(centroidhomfold_output_dir_path, output_file)
      centroidhomfold_params.insert(0, (rna_file_path, centroidhomfold_output_file_path, gamma_str))
      recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
      rec_seq_len = len(recs)
      turbofold_config_file_contents = "InSeq = {"
      for i in range(rec_seq_len):
        turbofold_config_file_contents += "%s/%d.fasta;" % (temp_dir_path, i)
      turbofold_config_file_contents += "}\nOutCT = {"
      for i in range(rec_seq_len):
        turbofold_config_file_contents += "%s/%d.ct;" % (temp_dir_path, i)
      turbofold_config_file_contents += "}\nIterations = 1\nMode = MEA\nMeaGamma = %f\nProcessors = %d" % (gamma, num_of_threads)
      turbofold_config_file_path = os.path.join(temp_dir_path, "turbofold_config.dat")
      turbofold_config_file = open(turbofold_config_file_path, "w")
      turbofold_config_file.write(turbofold_config_file_contents)
      turbofold_config_file.close()
      for (i, rec) in enumerate(recs):
        SeqIO.write([rec], open(os.path.join(temp_dir_path, "%d.fasta" % i), "w"), "fasta")
      begin = time.time()
      run_turbofold(turbofold_config_file_path)
      elapsed_time = time.time() - begin
      turbofold_elapsed_time += elapsed_time
      turbofold_output_file_contents = ""
      for i in range(rec_seq_len):
        ct_file_path = os.path.join(temp_dir_path, "%d.ct" % i)
        ss_string = read_ct_file(ct_file_path)
        turbofold_output_file_contents += ">%d\n%s\n\n" % (i, ss_string)
      turbofold_output_file_path = os.path.join(turbofold_output_dir_path, output_file)
      turbofold_output_file = open(turbofold_output_file_path, "w")
      turbofold_output_file.write(turbofold_output_file_contents)
  shutil.rmtree(temp_dir_path)
  pool = multiprocessing.Pool(num_of_threads)
  begin = time.time()
  pool.map(run_parasor, parasor_params)
  parasor_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_centroidhomfold, centroidhomfold_params)
  centroidhomfold_elapsed_time = time.time() - begin
  print("The elapsed time of the STRAP program and NeoFold program for a test set = %f [s]." % strap_and_neofold_elapsed_time)
  print("The elapsed time of the ParasoR program for a test set = %f [s]." % parasor_elapsed_time)
  print("The elapsed time of the CentroidHomFold program for a test set = %f [s]." % centroidhomfold_elapsed_time)
  print("The elapsed time of the TurboFold-smp program for a test set = %f [s]." % turbofold_elapsed_time)

def run_parasor(parasor_params):
  (rna_file_path, gamma, parasor_output_file_path) = parasor_params
  parasor_command = "ParasoR --pre --constraint 0 --struct=%f --input " % gamma + rna_file_path
  (output, _, _) = utils.run_command(parasor_command)
  lines = [line[10 :] for line in str(output).split("\\n") if line.startswith("#structure ")]
  parasor_output_file = open(parasor_output_file_path, "w+")
  parasor_output_buf = ""
  for (i, line) in enumerate(lines):
    parasor_output_buf += ">%d\n%s\n\n" % (i, line)
  parasor_output_file.write(parasor_output_buf)

def run_centroidhomfold(centroidhomfold_params):
  (rna_file_path, centroidhomfold_output_file_path, gamma_str) = centroidhomfold_params
  centroidhomfold_command = "centroid_homfold " + rna_file_path + " -H " + rna_file_path + " -g " + gamma_str
  (output, _, _) = utils.run_command(centroidhomfold_command)
  lines = [line.split()[0] for (i, line) in enumerate(str(output).split("\\n")) if i % 3 == 2]
  centroidhomfold_output_file = open(centroidhomfold_output_file_path, "w+")
  centroidhomfold_output_buf = ""
  for (i, line) in enumerate(lines):
    centroidhomfold_output_buf += ">%d\n%s\n\n" % (i, line)
  centroidhomfold_output_file.write(centroidhomfold_output_buf)

def run_turbofold(turbofold_config_file_path):
  turbofold_command = "TurboFold-smp " + turbofold_config_file_path
  utils.run_command(turbofold_command)

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
