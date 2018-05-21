#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  utils.init_matplotlib()
  strap_dir_path = asset_dir_path + "/strap"
  neofold_dir_path = asset_dir_path + "/neofold"
  if not os.path.isdir(strap_dir_path):
    os.mkdir(strap_dir_path)
  if not os.path.isdir(neofold_dir_path):
    os.mkdir(neofold_dir_path)
  rna_dir_path = asset_dir_path + "/rna_families"
  bpp_mat_file = "bpp_mats_on_sta.dat"
  ss_file = "ss.dat"
  for rna_file in os.listdir(rna_dir_path):
    if not rna_file.endswith(".fa"):
      continue
    rna_file_path = os.path.join(rna_dir_path, rna_file)
    (rna_familiy_name, extension) = os.path.splitext(rna_file)
    strap_output_dir_path = os.path.join(strap_dir_path, rna_familiy_name)
    strap_command = "strap -i " + rna_file_path + " -o " + strap_output_dir_path
    utils.run_command(strap_command)
    bpp_mat_file_path = os.path.join(strap_output_dir_path, bpp_mat_file)
    neofold_output_file_path = os.path.join(neofold_dir_path, "sss_of_" + rna_familiy_name + ".dat")
    neofold_command = "neofold -f " + rna_file_path + " -p " + bpp_mat_file_path + " -o " + neofold_output_file_path
    utils.run_command(neofold_command)


if __name__ == "__main__":
  main()
