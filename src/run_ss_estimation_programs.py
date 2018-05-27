#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  strap_dir_path = asset_dir_path + "/strap"
  neofold_dir_path = asset_dir_path + "/neofold"
  parasor_dir_path = asset_dir_path + "/parasor"
  centroidhomfold_dir_path = asset_dir_path + "/centroidhomfold"
  if not os.path.isdir(strap_dir_path):
    os.mkdir(strap_dir_path)
  if not os.path.isdir(neofold_dir_path):
    os.mkdir(neofold_dir_path)
  if not os.path.isdir(parasor_dir_path):
    os.mkdir(parasor_dir_path)
  if not os.path.isdir(centroidhomfold_dir_path):
    os.mkdir(centroidhomfold_dir_path)
  rna_dir_path = asset_dir_path + "/sampled_rna_families"
  bpp_mat_file = "bpp_mats_on_sta.dat"
  gammas = [2. ** i for i in range(-7, 11)]
  parasor_params = []
  centroidhomfold_params = []
  strap_and_neofold_elapsed_time = 0.
  for rna_file in os.listdir(rna_dir_path):
    if not rna_file.endswith(".fa"):
      continue
    rna_file_path = os.path.join(rna_dir_path, rna_file)
    (rna_familiy_name, extension) = os.path.splitext(rna_file)
    strap_output_dir_path = os.path.join(strap_dir_path, rna_familiy_name)
    strap_command = "strap --max_base_pairing_span 1000 --num_of_times_of_improvements_of_struct_align_prob_mat_quadruples 0 -i " + rna_file_path + " -o " + strap_output_dir_path
    begin = time.time()
    utils.run_command(strap_command)
    elapsed_time = time.time() - begin
    strap_and_neofold_elapsed_time += elapsed_time
    bpp_mat_file_path = os.path.join(strap_output_dir_path, bpp_mat_file)
    neofold_output_dir_path = os.path.join(neofold_dir_path, "sss_of_" + rna_familiy_name)
    parasor_output_dir_path = os.path.join(parasor_dir_path, "sss_of_" + rna_familiy_name)
    centroidhomfold_output_dir_path = os.path.join(centroidhomfold_dir_path, "sss_of_" + rna_familiy_name)
    if not os.path.isdir(neofold_output_dir_path):
      os.mkdir(neofold_output_dir_path)
    if not os.path.isdir(parasor_output_dir_path):
      os.mkdir(parasor_output_dir_path)
    if not os.path.isdir(centroidhomfold_output_dir_path):
      os.mkdir(centroidhomfold_output_dir_path)
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
  pool = multiprocessing.Pool(multiprocessing.cpu_count())
  parasor_elapsed_time = sum(pool.map(run_parasor, parasor_params))
  centroidhomfold_elapsed_time = sum(pool.map(run_centroidhomfold, centroidhomfold_params))
  print("The elapsed time of the STRAP program and NeoFold program for a test set = %f[s]." % strap_and_neofold_elapsed_time)
  print("The elapsed time of the ParasoR program for a test set = %f[s]." % parasor_elapsed_time)
  print("The elapsed time of the CentroidHomFold program for a test set = %f[s]." % centroidhomfold_elapsed_time)

def run_parasor(parasor_params):
  (rna_file_path, gamma, parasor_output_file_path) = parasor_params
  parasor_command = "ParasoR --pre --constraint 1000 --struct=%f --input " % gamma + rna_file_path
  begin = time.time()
  (output, _, _) = utils.run_command(parasor_command)
  elapsed_time = time.time() - begin
  lines = [line[10 :] for line in str(output).split("\\n") if line.startswith("#structure ")]
  parasor_output_file = open(parasor_output_file_path, "w+")
  parasor_output_buf = ""
  for (i, line) in enumerate(lines):
    parasor_output_buf += ">%d\n%s\n\n" % (i, line)
  parasor_output_file.write(parasor_output_buf)
  return elapsed_time

def run_centroidhomfold(centroidhomfold_params):
  (rna_file_path, centroidhomfold_output_file_path, gamma_str) = centroidhomfold_params
  centroidhomfold_command = "centroid_homfold " + rna_file_path + " -H " + rna_file_path + " -g " + gamma_str
  begin = time.time()
  (output, _, _) = utils.run_command(centroidhomfold_command)
  elapsed_time = time.time() - begin
  lines = [line.split()[0] for (i, line) in enumerate(str(output).split("\\n")) if i % 3 == 2]
  centroidhomfold_output_file = open(centroidhomfold_output_file_path, "w+")
  centroidhomfold_output_buf = ""
  for (i, line) in enumerate(lines):
    centroidhomfold_output_buf += ">%d\n%s\n\n" % (i, line)
  centroidhomfold_output_file.write(centroidhomfold_output_buf)
  return elapsed_time

if __name__ == "__main__":
  main()
