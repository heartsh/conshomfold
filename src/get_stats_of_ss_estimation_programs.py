#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import math

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  utils.init_matplotlib()
  neofold_ss_dir_path = asset_dir_path + "/neofold"
  parasor_ss_dir_path = asset_dir_path + "/parasor"
  rna_fam_dir_path = asset_dir_path + "/sampled_rna_families"
  neofold_ppvs = []
  neofold_senss = []
  parasor_ppvs = []
  parasor_senss = []
  gammas = [2. ** i for i in range(-10, 11)]
  (neofold_total_tp, neofold_total_tn, neofold_total_fp, neofold_total_fn) = (0, 0, 0, 0)
  (parasor_total_tp, parasor_total_tn, parasor_total_fp, parasor_total_fn) = (0, 0, 0, 0)
  for gamma in gammas:
    gamma_str = str(gamma)
    neofold_tp = neofold_tn = neofold_fp = neofold_fn = 0.
    parasor_tp = parasor_tn = parasor_fp = parasor_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_ss_file_path = os.path.join(rna_fam_dir_path, "sss_of_" + rna_fam_name + ".dat")
      ref_sss = utils.get_sss(utils.get_ss_strings(ref_ss_file_path))
      neofold_estimated_ss_dir_path = os.path.join(neofold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(neofold_estimated_ss_dir_path):
        continue
      parasor_estimated_ss_dir_path = os.path.join(parasor_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(parasor_estimated_ss_dir_path):
        continue
      neofold_estimated_ss_file_path = os.path.join(neofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss = utils.get_sss(utils.get_ss_strings(neofold_estimated_ss_file_path))
      for (estimated_ss, ref_ss, rna_seq_len) in zip(estimated_sss, ref_sss, rna_seq_lens):
        for i in range(0, rna_seq_len):
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                neofold_tp += 1
              else:
                neofold_tn += 1
            else:
              if estimated_bin == True:
                neofold_fp += 1
              else:
                neofold_fn += 1
      parasor_estimated_ss_file_path = os.path.join(parasor_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss = utils.get_sss(utils.get_ss_strings(parasor_estimated_ss_file_path))
      for (estimated_ss, ref_ss, rna_seq_len) in zip(estimated_sss, ref_sss, rna_seq_lens):
        for i in range(0, rna_seq_len):
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                parasor_tp += 1
              else:
                parasor_tn += 1
            else:
              if estimated_bin == True:
                parasor_fp += 1
              else:
                parasor_fn += 1
    neofold_total_tp += neofold_tp
    neofold_total_tn += neofold_tn
    neofold_total_fp += neofold_fp
    neofold_total_fn += neofold_fn
    ppv = neofold_tp / (neofold_tp + neofold_fp)
    sens = neofold_tp / (neofold_tp + neofold_fn)
    neofold_ppvs.insert(0, ppv)
    neofold_senss.insert(0, sens)
    parasor_total_tp += parasor_tp
    parasor_total_tn += parasor_tn
    parasor_total_fp += parasor_fp
    parasor_total_fn += parasor_fn
    ppv = parasor_tp / (parasor_tp + parasor_fp)
    sens = parasor_tp / (parasor_tp + parasor_fn)
    parasor_ppvs.insert(0, ppv)
    parasor_senss.insert(0, sens)
  neofold_ppvs = numpy.array(neofold_ppvs) 
  neofold_senss = numpy.array(neofold_senss)
  parasor_ppvs = numpy.array(parasor_ppvs) 
  parasor_senss = numpy.array(parasor_senss)
  lines = pyplot.plot(neofold_ppvs, neofold_senss, "o-", parasor_ppvs, parasor_senss, "v-")
  pyplot.xlabel("Sensitivity")
  pyplot.xlabel("PPV")
  pyplot.legend(handles = lines, loc = 1)
  neofold_mcc = (neofold_total_tp * neofold_total_tn - neofold_total_fp * neofold_total_fn) / math.sqrt((neofold_total_tp + neofold_total_fp) * (neofold_total_tp + neofold_total_fn) * (neofold_total_tn + neofold_total_fp) * (neofold_total_tn + neofold_total_fn))
  parasor_mcc = (parasor_total_tp * parasor_total_tn - parasor_total_fp * parasor_total_fn) / math.sqrt((parasor_total_tp + parasor_total_fp) * (parasor_total_tp + parasor_total_fn) * (parasor_total_tn + parasor_total_fp) * (parasor_total_tn + parasor_total_fn))
  print("The MCC of the NeoFold program for a test dataset = %f." % neofold_mcc)
  print("The MCC of the ParasoR program for a test dataset = %f." % parasor_mcc)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
