#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os
from sklearn.metrics import roc_curve
import math
from math import sqrt

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  consfold_ss_dir_path = asset_dir_path + "/consfold"
  bpp_consfold_ss_dir_path = asset_dir_path + "/bpp_consfold"
  centroidfold_ss_dir_path = asset_dir_path + "/centroidfold"
  centroidhomfold_ss_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_ss_dir_path = asset_dir_path + "/turbofold"
  rnafold_ss_dir_path = asset_dir_path + "/rnafold"
  rna_fam_dir_path = asset_dir_path + "/sampled_rna_families"
  consfold_ppvs = []
  consfold_senss = []
  consfold_fprs = []
  bpp_consfold_ppvs = []
  bpp_consfold_senss = []
  bpp_consfold_fprs = []
  centroidfold_ppvs = []
  centroidfold_senss = []
  centroidfold_fprs = []
  centroidhomfold_ppvs = []
  centroidhomfold_senss = []
  centroidhomfold_fprs = []
  turbofold_ppvs = []
  turbofold_senss = []
  turbofold_fprs = []
  rnafold_ppv = rnafold_sens = 0.
  gammas = [2. ** i for i in range(-7, 11)]
  consfold_f1_score = centroidfold_f1_score = centroidhomfold_f1_score = turbofold_f1_score = rnafold_f1_score = 0.
  consfold_mcc = centroidfold_mcc = centroidhomfold_mcc = turbofold_mcc = rnafold_mcc = 0.
  for gamma in gammas:
    gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
    consfold_tp = consfold_tn = consfold_fp = consfold_fn = 0.
    bpp_consfold_tp = bpp_consfold_tn = bpp_consfold_fp = bpp_consfold_fn = 0.
    centroidfold_tp = centroidfold_tn = centroidfold_fp = centroidfold_fn = 0.
    centroidhomfold_tp = centroidhomfold_tn = centroidhomfold_fp = centroidhomfold_fn = 0.
    turbofold_tp = turbofold_tn = turbofold_fp = turbofold_fn = 0.
    rnafold_tp = rnafold_tn = rnafold_fp = rnafold_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_ss_file_path = os.path.join(rna_fam_dir_path, "sss_of_" + rna_fam_name + ".dat")
      ref_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(ref_ss_file_path))
      consfold_estimated_ss_dir_path = os.path.join(consfold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(consfold_estimated_ss_dir_path):
        continue
      bpp_consfold_estimated_ss_dir_path = os.path.join(bpp_consfold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(bpp_consfold_estimated_ss_dir_path):
        continue
      centroidfold_estimated_ss_dir_path = os.path.join(centroidfold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(centroidfold_estimated_ss_dir_path):
        continue
      centroidhomfold_estimated_ss_dir_path = os.path.join(centroidhomfold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(centroidhomfold_estimated_ss_dir_path):
        continue
      turbofold_estimated_ss_dir_path = os.path.join(turbofold_ss_dir_path, "sss_of_" + rna_fam_name)
      if not os.path.isdir(turbofold_estimated_ss_dir_path):
        continue
      consfold_estimated_ss_file_path = os.path.join(consfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(consfold_estimated_ss_file_path))
      for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
        estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
        ref_ss, ref_flat_ss = ref_ss_and_flat_ss
        for i in range(0, rna_seq_len):
          estimated_bin = i in estimated_flat_ss
          ref_bin = i in ref_flat_ss
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              consfold_tn += 1
          else:
            if estimated_bin == True:
              consfold_fp += 1
            else:
              consfold_fn += 1
          for j in range(i + 1, rna_seq_len):
            estimated_bin = (i, j) in estimated_ss
            ref_bin = (i, j) in ref_ss
            if estimated_bin == ref_bin:
              if estimated_bin == True:
                consfold_tp += 1
      if gamma == 4.:
        ppv = consfold_tp / (consfold_tp + consfold_fp)
        sens = consfold_tp / (consfold_tp + consfold_fn)
        consfold_f1_score = 2 * ppv * sens / (ppv + sens)
        consfold_mcc = (consfold_tp * consfold_tn - consfold_fp * consfold_fn) / sqrt((consfold_tp + consfold_fp) * (consfold_tp + consfold_fn) * (consfold_tn + consfold_fp) * (consfold_tn + consfold_fn))
      if gamma >= 2.:
        bpp_consfold_estimated_ss_file_path = os.path.join(bpp_consfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".dat")
        estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(bpp_consfold_estimated_ss_file_path))
        for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
          estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
          ref_ss, ref_flat_ss = ref_ss_and_flat_ss
          for i in range(0, rna_seq_len):
            estimated_bin = i in estimated_flat_ss
            ref_bin = i in ref_flat_ss
            if estimated_bin == ref_bin:
              if estimated_bin == False:
                bpp_consfold_tn += 1
            else:
              if estimated_bin == True:
                bpp_consfold_fp += 1
              else:
                bpp_consfold_fn += 1
            for j in range(i + 1, rna_seq_len):
              estimated_bin = (i, j) in estimated_ss
              ref_bin = (i, j) in ref_ss
              if estimated_bin == ref_bin:
                if estimated_bin == True:
                  bpp_consfold_tp += 1
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
      if gamma == 4.:
        ppv = centroidfold_tp / (centroidfold_tp + centroidfold_fp)
        sens = centroidfold_tp / (centroidfold_tp + centroidfold_fn)
        centroidfold_f1_score = 2 * ppv * sens / (ppv + sens)
        centroidfold_mcc = (centroidfold_tp * centroidfold_tn - centroidfold_fp * centroidfold_fn) / sqrt((centroidfold_tp + centroidfold_fp) * (centroidfold_tp + centroidfold_fn) * (centroidfold_tn + centroidfold_fp) * (centroidfold_tn + centroidfold_fn))
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
      if gamma == 4.:
        ppv = centroidhomfold_tp / (centroidhomfold_tp + centroidhomfold_fp)
        sens = centroidhomfold_tp / (centroidhomfold_tp + centroidhomfold_fn)
        centroidhomfold_f1_score = 2 * ppv * sens / (ppv + sens)
        centroidhomfold_mcc = (centroidhomfold_tp * centroidhomfold_tn - centroidhomfold_fp * centroidhomfold_fn) / sqrt((centroidhomfold_tp + centroidhomfold_fp) * (centroidhomfold_tp + centroidhomfold_fn) * (centroidhomfold_tn + centroidhomfold_fp) * (centroidhomfold_tn + centroidhomfold_fn))
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
      if gamma == 4.:
        ppv = turbofold_tp / (turbofold_tp + turbofold_fp)
        sens = turbofold_tp / (turbofold_tp + turbofold_fn)
        turbofold_f1_score = 2 * ppv * sens / (ppv + sens)
        turbofold_mcc = (turbofold_tp * turbofold_tn - turbofold_fp * turbofold_fn) / sqrt((turbofold_tp + turbofold_fp) * (turbofold_tp + turbofold_fn) * (turbofold_tn + turbofold_fp) * (turbofold_tn + turbofold_fn))
      if gamma == 4.:
        rnafold_estimated_ss_file_path = os.path.join(rnafold_ss_dir_path, "sss_of_" + rna_fam_name + ".dat")
        estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(rnafold_estimated_ss_file_path))
        for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
          estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
          ref_ss, ref_flat_ss = ref_ss_and_flat_ss
          for i in range(0, rna_seq_len):
            estimated_bin = i in estimated_flat_ss
            ref_bin = i in ref_flat_ss
            if estimated_bin == ref_bin:
              if estimated_bin == False:
                rnafold_tn += 1
            else:
              if estimated_bin == True:
                rnafold_fp += 1
              else:
                rnafold_fn += 1
            for j in range(i + 1, rna_seq_len):
              estimated_bin = (i, j) in estimated_ss
              ref_bin = (i, j) in ref_ss
              if estimated_bin == ref_bin:
                if estimated_bin == True:
                  rnafold_tp += 1
        rnafold_ppv = rnafold_tp / (rnafold_tp + rnafold_fp)
        rnafold_sens = rnafold_tp / (rnafold_tp + rnafold_fn)
        rnafold_fpr = rnafold_fp / (rnafold_tn + rnafold_fp)
        rnafold_f1_score = 2 * rnafold_ppv * rnafold_sens / (rnafold_ppv + rnafold_sens)
        rnafold_mcc = (rnafold_tp * rnafold_tn - rnafold_fp * rnafold_fn) / sqrt((rnafold_tp + rnafold_fp) * (rnafold_tp + rnafold_fn) * (rnafold_tn + rnafold_fp) * (rnafold_tn + rnafold_fn))
    ppv = consfold_tp / (consfold_tp + consfold_fp)
    sens = consfold_tp / (consfold_tp + consfold_fn)
    fpr = consfold_fp / (consfold_tn + consfold_fp)
    consfold_ppvs.insert(0, ppv)
    consfold_senss.insert(0, sens)
    consfold_fprs.insert(0, fpr)
    if gamma >= 2.:
      ppv = bpp_consfold_tp / (bpp_consfold_tp + bpp_consfold_fp)
      sens = bpp_consfold_tp / (bpp_consfold_tp + bpp_consfold_fn)
      fpr = bpp_consfold_fp / (bpp_consfold_tn + bpp_consfold_fp)
      bpp_consfold_ppvs.insert(0, ppv)
      bpp_consfold_senss.insert(0, sens)
      bpp_consfold_fprs.insert(0, fpr)
    ppv = centroidfold_tp / (centroidfold_tp + centroidfold_fp)
    sens = centroidfold_tp / (centroidfold_tp + centroidfold_fn)
    fpr = centroidfold_fp / (centroidfold_tn + centroidfold_fp)
    centroidfold_ppvs.insert(0, ppv)
    centroidfold_senss.insert(0, sens)
    centroidfold_fprs.insert(0, fpr)
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
  consfold_ppvs = numpy.array(consfold_ppvs) 
  consfold_senss = numpy.array(consfold_senss)
  consfold_fprs = numpy.array(consfold_fprs)
  centroidfold_ppvs = numpy.array(centroidfold_ppvs) 
  centroidfold_senss = numpy.array(centroidfold_senss)
  centroidfold_fprs = numpy.array(centroidfold_fprs)
  centroidhomfold_ppvs = numpy.array(centroidhomfold_ppvs) 
  centroidhomfold_senss = numpy.array(centroidhomfold_senss)
  centroidhomfold_fprs = numpy.array(centroidhomfold_fprs)
  turbofold_ppvs = numpy.array(turbofold_ppvs) 
  turbofold_senss = numpy.array(turbofold_senss)
  turbofold_fprs = numpy.array(turbofold_fprs)
  line_1, = pyplot.plot(consfold_ppvs, consfold_senss, label = "ConsFold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(turbofold_ppvs, turbofold_senss, label = "TurboFold-smp", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_ppvs, centroidhomfold_senss, label = "CentroidHomFold", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(centroidfold_ppvs, centroidfold_senss, label = "CentroidFold", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(rnafold_ppv, rnafold_sens, label = "RNAfold", marker = "d", linestyle = "-")
  line_6, = pyplot.plot(bpp_consfold_ppvs, bpp_consfold_senss, label = "ConsFold (Conventional structure scores)", marker = "o", linestyle = "--")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6], loc = 3)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.figure()
  line_1, = pyplot.plot(consfold_fprs, consfold_senss, label = "ConsFold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(turbofold_fprs, turbofold_senss, label = "TurboFold-smp", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_fprs, centroidhomfold_senss, label = "CentroidHomFold", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(centroidfold_fprs, centroidfold_senss, label = "CentroidFold", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(rnafold_fpr, rnafold_sens, label = "RNAfold", marker = "d", linestyle = "-")
  line_6, = pyplot.plot(bpp_consfold_fprs, bpp_consfold_senss, label = "ConsFold (Conventional structure scores)", marker = "o", linestyle = "--")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6], loc = 4)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  print("ConsFold's MCC & F1 Score = %.3f & %.3f" %(consfold_mcc, consfold_f1_score))
  print("TurboFold's MCC & F1 Score = %.3f & %.3f" %(turbofold_mcc, turbofold_f1_score))
  print("CentroidHomFold's MCC & F1 Score = %.3f & %.3f" %(centroidhomfold_mcc, centroidhomfold_f1_score))
  print("CentroidFold's MCC & F1 Score = %.3f & %.3f" %(centroidfold_mcc, centroidfold_f1_score))
  print("RNAfold's MCC & F1 Score = %.3f & %.3f" %(rnafold_mcc, rnafold_f1_score))

if __name__ == "__main__":
  main()
