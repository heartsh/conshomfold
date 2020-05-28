#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import matplotlib
import os
import multiprocessing
import time
import datetime
import shutil
import community
import networkx
node_size = 150
edge_width = 2
label_font_size = 10 
edge_label_font_size = 4 
(_, axes) = pyplot.subplots(nrows = 2, ncols = 2, figsize = (12, 12))
cmap = pyplot.cm.viridis

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  seq_file_path = asset_dir_path + "/sampled_trnas.fa"
  seqs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  seq = seqs[0]
  seq_lens = [len(seq) for seq in seqs]
  seq_len = len(seq)
  max_upp_mat = utils.get_upp_mats(asset_dir_path + "/sampled_trnas/max_upp_mats.dat", seq_lens)[0]
  access_bpp_mat_4_2l = utils.get_bpp_mats(asset_dir_path + "/sampled_trnas/access_bpp_mats_on_2l.dat", seq_lens)[0]
  access_bpp_mat_4_ml = utils.get_bpp_mats(asset_dir_path + "/sampled_trnas/access_bpp_mats_on_ml.dat", seq_lens)[0]
  bpp_mat_4_el = utils.get_bpp_mats(asset_dir_path + "/sampled_trnas/bpp_mats_on_el.dat", seq_lens)[0]
  xlabels = [str(i + 1) if (i + 1) % 20 == 0 else "" for i in range(seq_len)]
  seaborn.heatmap(access_bpp_mat_4_2l, ax = axes[0][0], xticklabels = xlabels, yticklabels = xlabels, cbar = False, cmap = cmap)
  seaborn.heatmap(access_bpp_mat_4_ml, ax = axes[0][1], xticklabels = False, yticklabels = False, cbar = False, cmap = cmap)
  seaborn.heatmap(bpp_mat_4_el, ax = axes[1][0], xticklabels = False, yticklabels = False, cbar = False, cmap = cmap)
  upp_mat_4_2l = utils.get_upp_mats(asset_dir_path + "/sampled_trnas/upp_mats_on_2l.dat", seq_lens)[0]
  upp_mat_4_ml = utils.get_upp_mats(asset_dir_path + "/sampled_trnas/upp_mats_on_ml.dat", seq_lens)[0]
  upp_mat_4_el = utils.get_upp_mats(asset_dir_path + "/sampled_trnas/upp_mats_on_el.dat", seq_lens)[0]
  upp_mat_4_hl = utils.get_upp_mats(asset_dir_path + "/sampled_trnas/upp_mats_on_hl.dat", seq_lens)[0]
  pyplot.stackplot(range(seq_len), upp_mat_4_2l, upp_mat_4_ml, upp_mat_4_el, upp_mat_4_hl)
  legends = ["Unpairing accessible from 2-loop", "Unpairing accessible from multi-loop", "Unpairing accessible from external loop", "Unpairing accessible from 1-loop"]
  pyplot.legend(legends, loc = "upper right", bbox_to_anchor=(1.05, 1.2))
  pyplot.savefig(image_dir_path + "/trna_prob_dists.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
