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
  conshomfold_ppvs = []
  conshomfold_senss = []
  conshomfold_fprs = []
  centroidhomfold_ppvs = []
  centroidhomfold_senss = []
  centroidhomfold_fprs = []
  turbofold_ppvs = []
  turbofold_senss = []
  turbofold_fprs = []
  contrafold_ppvs = []
  contrafold_senss = []
  contrafold_fprs = []
  centroidfold_ppvs = []
  centroidfold_senss = []
  centroidfold_fprs = []
  rnafold_ppv = rnafold_sens = rnafold_fpr = 0.
  bpp_conshomfold_ppvs = []
  bpp_conshomfold_senss = []
  bpp_conshomfold_fprs = []
  gammas = [2. ** i for i in range(-5, 11)]
  # gammas = [2. ** i for i in range(-7, 11)]
  conshomfold_f1_score = centroidhomfold_f1_score = turbofold_f1_score = rnafold_f1_score = 0.
  conshomfold_mcc = centroidhomfold_mcc = turbofold_mcc = rnafold_mcc = 0.
  bpp_conshomfold_ss_dir_path = asset_dir_path + "/bpp_conshomfold"
  conshomfold_ss_dir_path = asset_dir_path + "/conshomfold"
  centroidhomfold_ss_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_ss_dir_path = asset_dir_path + "/turbofold"
  contrafold_ss_dir_path = asset_dir_path + "/contrafold"
  centroidfold_ss_dir_path = asset_dir_path + "/centroidfold"
  rnafold_ss_dir_path = asset_dir_path + "/rnafold"
  # rna_fam_dir_path = asset_dir_path + "/ref_sss"
  rna_fam_dir_path = asset_dir_path + "/ref_sss_4_micro_bench"
  for gamma in gammas:
    gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
    conshomfold_tp = conshomfold_tn = conshomfold_fp = conshomfold_fn = 0.
    bpp_conshomfold_tp = bpp_conshomfold_tn = bpp_conshomfold_fp = bpp_conshomfold_fn = 0.
    centroidhomfold_tp = centroidhomfold_tn = centroidhomfold_fp = centroidhomfold_fn = 0.
    turbofold_tp = turbofold_tn = turbofold_fp = turbofold_fn = 0.
    contrafold_tp = contrafold_tn = contrafold_fp = contrafold_fn = 0.
    centroidfold_tp = centroidfold_tn = centroidfold_fp = centroidfold_fn = 0.
    rnafold_tp = rnafold_tn = rnafold_fp = rnafold_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_ss_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      ref_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(ref_ss_file_path))
      conshomfold_estimated_ss_dir_path = os.path.join(conshomfold_ss_dir_path, rna_fam_name)
      bpp_conshomfold_estimated_ss_dir_path = os.path.join(bpp_conshomfold_ss_dir_path, rna_fam_name)
      centroidhomfold_estimated_ss_dir_path = os.path.join(centroidhomfold_ss_dir_path, rna_fam_name)
      turbofold_estimated_ss_dir_path = os.path.join(turbofold_ss_dir_path, rna_fam_name)
      contrafold_estimated_ss_dir_path = os.path.join(contrafold_ss_dir_path, rna_fam_name)
      centroidfold_estimated_ss_dir_path = os.path.join(centroidfold_ss_dir_path, rna_fam_name)
      conshomfold_estimated_ss_file_path = os.path.join(conshomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      tp, tn, fp, fn = get_pos_neg_counts(conshomfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
      conshomfold_tp += tp
      conshomfold_tn += tn
      conshomfold_fp += fp
      conshomfold_fn += fn
      if gamma >= 2.:
        bpp_conshomfold_estimated_ss_file_path = os.path.join(bpp_conshomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
        tp, tn, fp, fn = get_pos_neg_counts(bpp_conshomfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
        bpp_conshomfold_tp += tp
        bpp_conshomfold_tn += tn
        bpp_conshomfold_fp += fp
        bpp_conshomfold_fn += fn
      centroidhomfold_estimated_ss_file_path = os.path.join(centroidhomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      tp, tn, fp, fn = get_pos_neg_counts(centroidhomfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
      centroidhomfold_tp += tp
      centroidhomfold_tn += tn
      centroidhomfold_fp += fp
      centroidhomfold_fn += fn
      turbofold_estimated_ss_file_path = os.path.join(turbofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      tp, tn, fp, fn = get_pos_neg_counts(turbofold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
      turbofold_tp += tp
      turbofold_tn += tn
      turbofold_fp += fp
      turbofold_fn += fn
      contrafold_estimated_ss_file_path = os.path.join(contrafold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      tp, tn, fp, fn = get_pos_neg_counts(contrafold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
      contrafold_tp += tp
      contrafold_tn += tn
      contrafold_fp += fp
      contrafold_fn += fn
      centroidfold_estimated_ss_file_path = os.path.join(centroidfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      tp, tn, fp, fn = get_pos_neg_counts(centroidfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
      centroidfold_tp += tp
      centroidfold_tn += tn
      centroidfold_fp += fp
      centroidfold_fn += fn
      if gamma == 1.:
        rnafold_estimated_ss_file_path = os.path.join(rnafold_ss_dir_path, rna_fam_name + ".fa")
        tp, tn, fp, fn = get_pos_neg_counts(rnafold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens)
        rnafold_tp += tp
        rnafold_tn += tn
        rnafold_fp += fp
        rnafold_fn += fn
    ppv = get_ppv(conshomfold_tp, conshomfold_fp)
    sens = get_sens(conshomfold_tp, conshomfold_fn)
    fpr = get_fpr(conshomfold_tn, conshomfold_fp)
    conshomfold_ppvs.insert(0, ppv)
    conshomfold_senss.insert(0, sens)
    conshomfold_fprs.insert(0, fpr)
    if gamma == 1.:
      conshomfold_f1_score = get_f1_score(ppv, sens)
      conshomfold_mcc = get_mcc(conshomfold_tp, conshomfold_tn, conshomfold_fp, conshomfold_fn)
    if gamma >= 2.:
      ppv = get_ppv(bpp_conshomfold_tp, bpp_conshomfold_fp)
      sens = get_sens(bpp_conshomfold_tp, bpp_conshomfold_fn)
      fpr = get_fpr(bpp_conshomfold_tn, bpp_conshomfold_fp)
      bpp_conshomfold_ppvs.insert(0, ppv)
      bpp_conshomfold_senss.insert(0, sens)
      bpp_conshomfold_fprs.insert(0, fpr)
    ppv = get_ppv(centroidhomfold_tp, centroidhomfold_fp)
    sens = get_sens(centroidhomfold_tp, centroidhomfold_fn)
    fpr = get_fpr(centroidhomfold_tn, centroidhomfold_fp)
    centroidhomfold_ppvs.insert(0, ppv)
    centroidhomfold_senss.insert(0, sens)
    centroidhomfold_fprs.insert(0, fpr)
    if gamma == 1.:
      centroidhomfold_f1_score = get_f1_score(ppv, sens)
      centroidhomfold_mcc = get_mcc(centroidhomfold_tp, centroidhomfold_tn, centroidhomfold_fp, centroidhomfold_fn)
    ppv = get_ppv(turbofold_tp, turbofold_fp)
    sens = get_sens(turbofold_tp, turbofold_fn)
    fpr = get_fpr(turbofold_tn, turbofold_fp)
    turbofold_ppvs.insert(0, ppv)
    turbofold_senss.insert(0, sens)
    turbofold_fprs.insert(0, fpr)
    if gamma == 1.:
      turbofold_f1_score = get_f1_score(ppv, sens)
      turbofold_mcc = get_mcc(turbofold_tp, turbofold_tn, turbofold_fp, turbofold_fn)
    ppv = get_ppv(contrafold_tp, contrafold_fp)
    sens = get_sens(contrafold_tp, contrafold_fn)
    fpr = get_fpr(contrafold_tn, contrafold_fp)
    contrafold_ppvs.insert(0, ppv)
    contrafold_senss.insert(0, sens)
    contrafold_fprs.insert(0, fpr)
    if gamma == 1.:
      contrafold_f1_score = get_f1_score(ppv, sens)
      contrafold_mcc = get_mcc(contrafold_tp, contrafold_tn, contrafold_fp, contrafold_fn)
    ppv = get_ppv(centroidfold_tp, centroidfold_fp)
    sens = get_sens(centroidfold_tp, centroidfold_fn)
    fpr = get_fpr(centroidfold_tn, centroidfold_fp)
    centroidfold_ppvs.insert(0, ppv)
    centroidfold_senss.insert(0, sens)
    centroidfold_fprs.insert(0, fpr)
    if gamma == 1.:
      centroidfold_f1_score = get_f1_score(ppv, sens)
      centroidfold_mcc = get_mcc(centroidfold_tp, centroidfold_tn, centroidfold_fp, centroidfold_fn)
    if gamma == 1.:
      rnafold_ppv = get_ppv(rnafold_tp, rnafold_fp)
      rnafold_sens = get_sens(rnafold_tp, rnafold_fn)
      rnafold_fpr = get_fpr(rnafold_tn, rnafold_fp)
      rnafold_f1_score = get_f1_score(rnafold_ppv, rnafold_sens)
      rnafold_mcc = get_mcc(rnafold_tp, rnafold_tn, rnafold_fp, rnafold_fn)
  conshomfold_ppvs = numpy.array(conshomfold_ppvs)
  conshomfold_senss = numpy.array(conshomfold_senss)
  conshomfold_fprs = numpy.array(conshomfold_fprs)
  centroidhomfold_ppvs = numpy.array(centroidhomfold_ppvs)
  centroidhomfold_senss = numpy.array(centroidhomfold_senss)
  centroidhomfold_fprs = numpy.array(centroidhomfold_fprs)
  urbofold_ppvs = numpy.array(turbofold_ppvs) 
  turbofold_senss = numpy.array(turbofold_senss)
  turbofold_fprs = numpy.array(turbofold_fprs)
  bpp_conshomfold_ppvs = numpy.array(bpp_conshomfold_ppvs)
  bpp_conshomfold_senss = numpy.array(bpp_conshomfold_senss)
  bpp_conshomfold_fprs = numpy.array(bpp_conshomfold_fprs)
  line_1, = pyplot.plot(conshomfold_ppvs, conshomfold_senss, label = "ConsHomfold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(turbofold_ppvs, turbofold_senss, label = "TurboFold", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_ppvs, centroidhomfold_senss, label = "CentroidHomfold", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(contrafold_ppvs, contrafold_senss, label = "CONTRAfold", marker = "h", linestyle = "-")
  line_5, = pyplot.plot(centroidfold_ppvs, centroidfold_senss, label = "CentroidFold", marker = "v", linestyle = "-")
  line_6, = pyplot.plot(rnafold_ppv, rnafold_sens, label = "RNAfold", marker = "d", linestyle = "-")
  line_7, = pyplot.plot(bpp_conshomfold_ppvs, bpp_conshomfold_senss, label = "ConsHomfold (Conventional)", marker = "o", linestyle = ":")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6, line_7], loc = "lower left")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  line_1, = pyplot.plot(conshomfold_fprs, conshomfold_senss, label = "ConsHomfold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(turbofold_fprs, turbofold_senss, label = "TurboFold", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_fprs, centroidhomfold_senss, label = "CentroidHomfold", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(contrafold_fprs, contrafold_senss, label = "CONTRAfold", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(centroidfold_fprs, centroidfold_senss, label = "Centroidfold", marker = "h", linestyle = "-")
  line_6, = pyplot.plot(rnafold_fpr, rnafold_sens, label = "RNAfold", marker = "d", linestyle = "-")
  line_7, = pyplot.plot(bpp_conshomfold_fprs, bpp_conshomfold_senss, label = "ConsHomfold (Conventional)", marker = "o", linestyle = ":")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  print("ConsHomfold's MCC & F1 score = %.3f & %.3f" %(conshomfold_mcc, conshomfold_f1_score))
  print("TurboFold's MCC & F1 score = %.3f & %.3f" %(turbofold_mcc, turbofold_f1_score))
  print("CentroidHomfold's MCC & F1 score = %.3f & %.3f" %(centroidhomfold_mcc, centroidhomfold_f1_score))
  print("CONTRAfold's MCC & F1 score = %.3f & %.3f" %(contrafold_mcc, contrafold_f1_score))
  print("Centroidfold's MCC & F1 score = %.3f & %.3f" %(centroidfold_mcc, centroidfold_f1_score))
  print("RNAfold's MCC & F1 score = %.3f & %.3f" %(rnafold_mcc, rnafold_f1_score))

def get_pos_neg_counts(estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens):
  tp = tn = fp = fn = 0
  estimated_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(estimated_ss_file_path))
  for (estimated_ss_and_flat_ss, ref_ss_and_flat_ss, rna_seq_len) in zip(estimated_sss_and_flat_sss, ref_sss_and_flat_sss, rna_seq_lens):
    estimated_ss, estimated_flat_ss = estimated_ss_and_flat_ss
    ref_ss, ref_flat_ss = ref_ss_and_flat_ss
    for i in range(0, rna_seq_len):
      estimated_bin = i in estimated_flat_ss
      ref_bin = i in ref_flat_ss
      if estimated_bin == ref_bin:
        if estimated_bin == False:
          tn += 1
      else:
        if estimated_bin == True:
          fp += 1
        else:
          fn += 1
      for j in range(i + 1, rna_seq_len):
        estimated_bin = (i, j) in estimated_ss
        ref_bin = (i, j) in ref_ss
        if estimated_bin == ref_bin:
          if estimated_bin == True:
            tp += 1
  return tp, tn, fp, fn

def get_f1_score(ppv, sens):
  return 2 * ppv * sens / (ppv + sens)

def get_mcc(tp, tn, fp, fn):
  return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

def get_ppv(tp, fp):
  return tp / (tp + fp)

def get_sens(tp, fn):
  return tp / (tp + fn)

def get_fpr(tn, fp):
  return fp / (tn + fp)

if __name__ == "__main__":
  main()
