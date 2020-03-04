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
  phylofold_ss_dir_path = asset_dir_path + "/phylofold"
  centroidfold_ss_dir_path = asset_dir_path + "/centroidfold"
  contrafold_ss_dir_path = asset_dir_path + "/contrafold"
  centroidhomfold_ss_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_ss_dir_path = asset_dir_path + "/turbofold"
  rna_fam_dir_path = asset_dir_path + "/sampled_rna_families"
  phylofold_ppvs = []
  phylofold_senss = []
  phylofold_fprs = []
  centroidfold_ppvs = []
  centroidfold_senss = []
  centroidfold_fprs = []
  contrafold_ppvs = []
  contrafold_senss = []
  contrafold_fprs = []
  centroidhomfold_ppvs = []
  centroidhomfold_senss = []
  centroidhomfold_fprs = []
  turbofold_ppvs = []
  turbofold_senss = []
  turbofold_fprs = []
  gammas = [2. ** i for i in range(-7, 11)]
  for gamma in gammas:
    gamma_str = str(gamma)
    phylofold_tp = phylofold_tn = phylofold_fp = phylofold_fn = 0.
    centroidfold_tp = centroidfold_tn = centroidfold_fp = centroidfold_fn = 0.
    contrafold_tp = contrafold_tn = contrafold_fp = contrafold_fn = 0.
    centroidhomfold_tp = centroidhomfold_tn = centroidhomfold_fp = centroidhomfold_fn = 0.
    turbofold_tp = turbofold_tn = turbofold_fp = turbofold_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_ss_file_path = os.path.join(rna_fam_dir_path, "sss_of_" + rna_fam_name + ".dat")
      ref_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(ref_ss_file_path))
      phylofold_estimated_ss_dir_path = os.path.join(phylofold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(phylofold_estimated_ss_dir_path):
        continue
      centroidfold_estimated_ss_dir_path = os.path.join(centroidfold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(centroidfold_estimated_ss_dir_path):
        continue
      contrafold_estimated_ss_dir_path = os.path.join(contrafold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(contrafold_estimated_ss_dir_path):
        continue
      centroidhomfold_estimated_ss_dir_path = os.path.join(centroidhomfold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(centroidhomfold_estimated_ss_dir_path):
        continue
      turbofold_estimated_ss_dir_path = os.path.join(turbofold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(turbofold_estimated_ss_dir_path):
        continue
      phylofold_estimated_ss_file_path = os.path.join(phylofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(phylofold_estimated_ss_file_path))
      for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
        estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
        ref_ss, ref_flat_ss = ref_ss_and_flat_ss
        for i in range(0, rna_seq_len):
          estimated_bin = i in estimated_flat_ss
          ref_bin = i in ref_flat_ss
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              phylofold_tn += 1
          else:
            if estimated_bin == True:
              phylofold_fp += 1
            else:
              phylofold_fn += 1
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                phylofold_tp += 1
      centroidfold_estimated_ss_file_path = os.path.join(centroidfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(centroidfold_estimated_ss_file_path))
      for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
        estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
        ref_ss, ref_flat_ss = ref_ss_and_flat_ss
        for i in range(0, rna_seq_len):
          estimated_bin = i in estimated_flat_ss
          ref_bin = i in ref_flat_ss
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              centroidfold_tn += 1
          else:
            if estimated_bin == True:
              centroidfold_fp += 1
            else:
              centroidfold_fn += 1
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                centroidfold_tp += 1
      contrafold_estimated_ss_file_path = os.path.join(contrafold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(contrafold_estimated_ss_file_path))
      for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
        estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
        ref_ss, ref_flat_ss = ref_ss_and_flat_ss
        for i in range(0, rna_seq_len):
          estimated_bin = i in estimated_flat_ss
          ref_bin = i in ref_flat_ss
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              contrafold_tn += 1
          else:
            if estimated_bin == True:
              contrafold_fp += 1
            else:
              contrafold_fn += 1
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                contrafold_tp += 1
      centroidhomfold_estimated_ss_file_path = os.path.join(centroidhomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(centroidhomfold_estimated_ss_file_path))
      for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
        estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
        ref_ss, ref_flat_ss = ref_ss_and_flat_ss
        for i in range(0, rna_seq_len):
          estimated_bin = i in estimated_flat_ss
          ref_bin = i in ref_flat_ss
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              centroidhomfold_tn += 1
          else:
            if estimated_bin == True:
              centroidhomfold_fp += 1
            else:
              centroidhomfold_fn += 1
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                centroidhomfold_tp += 1
      turbofold_estimated_ss_file_path = os.path.join(turbofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(turbofold_estimated_ss_file_path))
      for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
        estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
        ref_ss, ref_flat_ss = ref_ss_and_flat_ss
        for i in range(0, rna_seq_len):
          estimated_bin = i in estimated_flat_ss
          ref_bin = i in ref_flat_ss
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              turbofold_tn += 1
          else:
            if estimated_bin == True:
              turbofold_fp += 1
            else:
              turbofold_fn += 1
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                turbofold_tp += 1
    ppv = phylofold_tp / (phylofold_tp + phylofold_fp)
    sens = phylofold_tp / (phylofold_tp + phylofold_fn)
    fpr = phylofold_fp / (phylofold_tn + phylofold_fp)
    phylofold_ppvs.insert(0, ppv)
    phylofold_senss.insert(0, sens)
    phylofold_fprs.insert(0, fpr)
    ppv = centroidfold_tp / (centroidfold_tp + centroidfold_fp)
    sens = centroidfold_tp / (centroidfold_tp + centroidfold_fn)
    fpr = centroidfold_fp / (centroidfold_tn + centroidfold_fp)
    centroidfold_ppvs.insert(0, ppv)
    centroidfold_senss.insert(0, sens)
    centroidfold_fprs.insert(0, fpr)
    ppv = contrafold_tp / (contrafold_tp + contrafold_fp)
    sens = contrafold_tp / (contrafold_tp + contrafold_fn)
    fpr = contrafold_fp / (contrafold_tn + contrafold_fp)
    contrafold_ppvs.insert(0, ppv)
    contrafold_senss.insert(0, sens)
    contrafold_fprs.insert(0, fpr)
    ppv = centroidhomfold_tp / (centroidhomfold_tp + centroidhomfold_fp)
    sens = centroidhomfold_tp / (centroidhomfold_tp + centroidhomfold_fn)
    fpr = centroidhomfold_fp / (centroidhomfold_tn + centroidhomfold_fp)
    centroidhomfold_ppvs.insert(0, ppv)
    centroidhomfold_senss.insert(0, sens)
    centroidhomfold_fprs.insert(0, fpr)
    ppv = turbofold_tp / (turbofold_tp + turbofold_fp)
    sens = turbofold_tp / (turbofold_tp + turbofold_fn)
    fpr = turbofold_fp / (turbofold_tn + turbofold_fp)
    turbofold_ppvs.insert(0, ppv)
    turbofold_senss.insert(0, sens)
    turbofold_fprs.insert(0, fpr)
  phylofold_ppvs = numpy.array(phylofold_ppvs) 
  phylofold_senss = numpy.array(phylofold_senss)
  phylofold_fprs = numpy.array(phylofold_fprs)
  centroidfold_ppvs = numpy.array(centroidfold_ppvs) 
  centroidfold_senss = numpy.array(centroidfold_senss)
  centroidfold_fprs = numpy.array(centroidfold_fprs)
  contrafold_ppvs = numpy.array(contrafold_ppvs) 
  contrafold_senss = numpy.array(contrafold_senss)
  contrafold_fprs = numpy.array(contrafold_fprs)
  centroidhomfold_ppvs = numpy.array(centroidhomfold_ppvs) 
  centroidhomfold_senss = numpy.array(centroidhomfold_senss)
  centroidhomfold_fprs = numpy.array(centroidhomfold_fprs)
  turbofold_ppvs = numpy.array(turbofold_ppvs) 
  turbofold_senss = numpy.array(turbofold_senss)
  turbofold_fprs = numpy.array(turbofold_fprs)
  line_1, = pyplot.plot(phylofold_senss, phylofold_ppvs, label = "PhyloFold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(turbofold_senss, turbofold_ppvs, label = "TurboFold-smp", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_senss, centroidhomfold_ppvs, label = "CentroidHomFold", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(centroidfold_senss, centroidfold_ppvs, label = "CentroidFold (CentroidFold)", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(contrafold_senss, contrafold_ppvs, label = "CentroidFold (CONTRAfold)", marker = "^", linestyle = "-")
  pyplot.xlabel("Sensitivity")
  pyplot.ylabel("Positive predictive value")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5], loc = 3)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/senss_vs_ppvs_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.figure()
  line_1, = pyplot.plot(phylofold_fprs, phylofold_senss, label = "PhyloFold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(turbofold_fprs, turbofold_senss, label = "TurboFold-smp", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_fprs, centroidhomfold_senss, label = "CentroidHomFold", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(centroidfold_fprs, centroidfold_senss, label = "CentroidFold (CentroidFold)", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(contrafold_fprs, contrafold_senss, label = "CentroidFold (CONTRAfold)", marker = "^", linestyle = "-")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5], loc = 4)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
