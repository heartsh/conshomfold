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
pyplot.figure(figsize = (12, 12))
cmap = pyplot.cm.viridis

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  # Analysis for sampled tRNAs.
  seq_file_path = asset_dir_path + "/sampled_trnas.fa"
  seqs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  seq = seqs[0]
  seq_lens = [len(seq) for seq in seqs]
  seq_len = len(seq)
  stat_dir_path = asset_dir_path + "/sampled_trnas"
  ss_file_path = stat_dir_path + "/gamma=1024.fa"
  sss = [rec for rec in SeqIO.parse(ss_file_path, "fasta")];
  (ss, flat_ss) = utils.get_ss_and_flat_ss(sss[0])
  max_bpp_mat = utils.get_bpp_mats(stat_dir_path + "/max_bpp_mats.dat", seq_lens)[0]
  max_upp_mat = utils.get_upp_mats(stat_dir_path + "/max_upp_mats.dat", seq_lens)[0]
  access_bpp_mat_4_2l = utils.get_bpp_mats(stat_dir_path + "/access_bpp_mats_on_2l.dat", seq_lens)[0]
  access_bpp_mat_4_ml = utils.get_bpp_mats(stat_dir_path + "/access_bpp_mats_on_ml.dat", seq_lens)[0]
  bpp_mat_4_el = utils.get_bpp_mats(stat_dir_path + "/bpp_mats_on_el.dat", seq_lens)[0]
  graph = networkx.Graph()
  graph.add_nodes_from([i for i in range(seq_len + 5)])
  pos = networkx.circular_layout(graph)
  labels = {}
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_2l[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "solid", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_ml[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dashed", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == bpp_mat_4_el[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dotted", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  upp_mat_4_2l = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_2l.dat", seq_lens)[0]
  upp_mat_4_ml = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_ml.dat", seq_lens)[0]
  upp_mat_4_el = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_el.dat", seq_lens)[0]
  upp_mat_4_hl = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_hl.dat", seq_lens)[0]
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    labels[i] = seq[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_2l[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "o", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_ml[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "p", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_el[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "v", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_hl[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "d", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  edges = []
  edge_weights = []
  networkx.draw_networkx_labels(graph, pos, labels = labels, font_color = "r", font_size = label_font_size)
  pyplot.colorbar(pyplot.cm.ScalarMappable(norm = matplotlib.colors.Normalize(), cmap = cmap), fraction = 0.01)
  legends = ["Pairing accessible from 2-loop", "Pairing accessible from multi-loop", "Pairing accessible from external loop", "Unpairing accessible from 2-loop", "Unpairing accessible from multi-loop", "Unpairing accessible from external loop", "Unpairing accessible from 1-loop"]
  pyplot.legend(legends, loc = "center")
  pyplot.axis("off")
  pyplot.savefig(image_dir_path + "/comm_struct_1.eps", bbox_inches = "tight")
  # Analysis for pri miR-16-2.
  pyplot.clf()
  seq_file_path = asset_dir_path + "/homologs_of_pri_miR_16_2.fa"
  seqs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  seq = seqs[0]
  seq_lens = [len(seq) for seq in seqs]
  seq_len = len(seq)
  stat_dir_path = asset_dir_path + "/homologs_of_pri_miR_16_2"
  ss_file_path = stat_dir_path + "/gamma=8.fa"
  sss = [rec for rec in SeqIO.parse(ss_file_path, "fasta")];
  (ss, flat_ss) = utils.get_ss_and_flat_ss(sss[0])
  max_bpp_mat = utils.get_bpp_mats(stat_dir_path + "/max_bpp_mats.dat", seq_lens)[0]
  max_upp_mat = utils.get_upp_mats(stat_dir_path + "/max_upp_mats.dat", seq_lens)[0]
  access_bpp_mat_4_2l = utils.get_bpp_mats(stat_dir_path + "/access_bpp_mats_on_2l.dat", seq_lens)[0]
  access_bpp_mat_4_ml = utils.get_bpp_mats(stat_dir_path + "/access_bpp_mats_on_ml.dat", seq_lens)[0]
  bpp_mat_4_el = utils.get_bpp_mats(stat_dir_path + "/bpp_mats_on_el.dat", seq_lens)[0]
  graph = networkx.Graph()
  graph.add_nodes_from([i for i in range(seq_len + 5)])
  pos = networkx.circular_layout(graph)
  labels = {}
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_2l[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "solid", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_ml[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dashed", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == bpp_mat_4_el[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dotted", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  pyplot.colorbar(pyplot.cm.ScalarMappable(norm = matplotlib.colors.Normalize(), cmap = cmap), fraction = 0.01)
  upp_mat_4_2l = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_2l.dat", seq_lens)[0]
  upp_mat_4_ml = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_ml.dat", seq_lens)[0]
  upp_mat_4_el = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_el.dat", seq_lens)[0]
  upp_mat_4_hl = utils.get_upp_mats(stat_dir_path + "/upp_mats_on_hl.dat", seq_lens)[0]
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    labels[i] = seq[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_2l[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "o", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_ml[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "p", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_el[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "v", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_hl[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "d", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  edges = []
  edge_weights = []
  networkx.draw_networkx_labels(graph, pos, labels = labels, font_color = "r", font_size = label_font_size)
  pyplot.axis("off")
  pyplot.savefig(image_dir_path + "/comm_struct_2.eps", bbox_inches = "tight")
  pyplot.clf()
  stat_dir_path = asset_dir_path + "/homologs_of_pri_miR_16_2"
  ss_file_path = asset_dir_path + "/ref_struct_of_pri_miR_16_2.fa"
  sss = [rec for rec in SeqIO.parse(ss_file_path, "fasta")];
  (ss, flat_ss) = utils.get_ss_and_flat_ss(sss[0])
  labels = {}
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_2l[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "solid", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_ml[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dashed", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in ss:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == bpp_mat_4_el[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dotted", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  pyplot.colorbar(pyplot.cm.ScalarMappable(norm = matplotlib.colors.Normalize(), cmap = cmap), fraction = 0.01)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    labels[i] = seq[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_2l[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "o", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_ml[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "p", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_el[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "v", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(seq_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_ss) and max_upp == upp_mat_4_hl[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "d", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  edges = []
  edge_weights = []
  networkx.draw_networkx_labels(graph, pos, labels = labels, font_color = "r", font_size = label_font_size)
  pyplot.axis("off")
  pyplot.savefig(image_dir_path + "/comm_struct_3.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
