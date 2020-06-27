#! /usr/bin/env python

import utils
from Bio import SeqIO
import seaborn
from matplotlib import pyplot
import os
from sklearn.metrics import roc_curve
import math
from math import sqrt
import multiprocessing

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  conshomfold_ppvs = []
  conshomfold_senss = []
  conshomfold_fprs = []
  conshomfold_f1_scores = []
  conshomfold_mccs = []
  centroidhomfold_ppvs = []
  centroidhomfold_senss = []
  centroidhomfold_fprs = []
  centroidhomfold_f1_scores = []
  centroidhomfold_mccs = []
  turbofold_ppvs = []
  turbofold_senss = []
  turbofold_fprs = []
  turbofold_f1_scores = []
  turbofold_mccs = []
  contrafold_ppvs = []
  contrafold_senss = []
  contrafold_fprs = []
  contrafold_f1_scores = []
  contrafold_mccs = []
  centroidfold_ppvs = []
  centroidfold_senss = []
  centroidfold_fprs = []
  centroidfold_f1_scores = []
  centroidfold_mccs = []
  bpp_conshomfold_ppvs = []
  bpp_conshomfold_senss = []
  bpp_conshomfold_fprs = []
  bpp_conshomfold_f1_scores = []
  bpp_conshomfold_mccs = []
  rnafold_ppv = rnafold_sens = rnafold_fpr = rnafold_f1_score = rnafold_mcc = 0.
  gammas = [2. ** i for i in range(-7, 11)]
  bpp_conshomfold_ss_dir_path = asset_dir_path + "/bpp_conshomfold"
  conshomfold_ss_dir_path = asset_dir_path + "/conshomfold"
  centroidhomfold_ss_dir_path = asset_dir_path + "/centroidhomfold"
  turbofold_ss_dir_path = asset_dir_path + "/turbofold"
  contrafold_ss_dir_path = asset_dir_path + "/contrafold"
  centroidfold_ss_dir_path = asset_dir_path + "/centroidfold"
  rnafold_ss_dir_path = asset_dir_path + "/rnafold"
  rna_fam_dir_path = asset_dir_path + "/ref_sss"
  # rna_fam_dir_path = asset_dir_path + "/ref_sss_4_micro_bench"
  pool = multiprocessing.Pool(num_of_threads)
  for gamma in gammas:
    conshomfold_count_params = []
    centroidhomfold_count_params = []
    turbofold_count_params = []
    contrafold_count_params = []
    centroidfold_count_params = []
    rnafold_count_params = []
    bpp_conshomfold_count_params = []
    gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
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
      conshomfold_count_params.insert(0, (conshomfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      bpp_conshomfold_estimated_ss_file_path = os.path.join(bpp_conshomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      bpp_conshomfold_count_params.insert(0, (bpp_conshomfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      centroidhomfold_estimated_ss_file_path = os.path.join(centroidhomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      centroidhomfold_count_params.insert(0, (centroidhomfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      turbofold_estimated_ss_file_path = os.path.join(turbofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      turbofold_count_params.insert(0, (turbofold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      contrafold_estimated_ss_file_path = os.path.join(contrafold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      contrafold_count_params.insert(0, (contrafold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      centroidfold_estimated_ss_file_path = os.path.join(centroidfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      centroidfold_count_params.insert(0, (centroidfold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
      if gamma == 1.:
        rnafold_estimated_ss_file_path = os.path.join(rnafold_ss_dir_path, rna_fam_name + ".fa")
        rnafold_count_params.insert(0, (rnafold_estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens))
    results = pool.map(get_pos_neg_counts, conshomfold_count_params)
    conshomfold_tp, conshomfold_tn, conshomfold_fp, conshomfold_fn = final_sum(results)
    ppv = get_ppv(conshomfold_tp, conshomfold_fp)
    sens = get_sens(conshomfold_tp, conshomfold_fn)
    fpr = get_fpr(conshomfold_tn, conshomfold_fp)
    conshomfold_ppvs.insert(0, ppv)
    conshomfold_senss.insert(0, sens)
    conshomfold_fprs.insert(0, fpr)
    conshomfold_f1_scores.append(get_f1_score(ppv, sens))
    conshomfold_mccs.append(get_mcc(conshomfold_tp, conshomfold_tn, conshomfold_fp, conshomfold_fn))
    results = pool.map(get_pos_neg_counts, bpp_conshomfold_count_params)
    bpp_conshomfold_tp, bpp_conshomfold_tn, bpp_conshomfold_fp, bpp_conshomfold_fn = final_sum(results)
    ppv = get_ppv(bpp_conshomfold_tp, bpp_conshomfold_fp)
    sens = get_sens(bpp_conshomfold_tp, bpp_conshomfold_fn)
    fpr = get_fpr(bpp_conshomfold_tn, bpp_conshomfold_fp)
    bpp_conshomfold_ppvs.insert(0, ppv)
    bpp_conshomfold_senss.insert(0, sens)
    bpp_conshomfold_fprs.insert(0, fpr)
    bpp_conshomfold_f1_scores.append(get_f1_score(ppv, sens))
    bpp_conshomfold_mccs.append(get_mcc(bpp_conshomfold_tp, bpp_conshomfold_tn, bpp_conshomfold_fp, bpp_conshomfold_fn))
    results = pool.map(get_pos_neg_counts, centroidhomfold_count_params)
    centroidhomfold_tp, centroidhomfold_tn, centroidhomfold_fp, centroidhomfold_fn = final_sum(results)
    ppv = get_ppv(centroidhomfold_tp, centroidhomfold_fp)
    sens = get_sens(centroidhomfold_tp, centroidhomfold_fn)
    fpr = get_fpr(centroidhomfold_tn, centroidhomfold_fp)
    centroidhomfold_ppvs.insert(0, ppv)
    centroidhomfold_senss.insert(0, sens)
    centroidhomfold_fprs.insert(0, fpr)
    centroidhomfold_f1_scores.append(get_f1_score(ppv, sens))
    centroidhomfold_mccs.append(get_mcc(centroidhomfold_tp, centroidhomfold_tn, centroidhomfold_fp, centroidhomfold_fn))
    results = pool.map(get_pos_neg_counts, turbofold_count_params)
    turbofold_tp, turbofold_tn, turbofold_fp, turbofold_fn = final_sum(results)
    ppv = get_ppv(turbofold_tp, turbofold_fp)
    sens = get_sens(turbofold_tp, turbofold_fn)
    fpr = get_fpr(turbofold_tn, turbofold_fp)
    turbofold_ppvs.insert(0, ppv)
    turbofold_senss.insert(0, sens)
    turbofold_fprs.insert(0, fpr)
    turbofold_f1_scores.append(get_f1_score(ppv, sens))
    turbofold_mccs.append(get_mcc(turbofold_tp, turbofold_tn, turbofold_fp, turbofold_fn))
    results = pool.map(get_pos_neg_counts, contrafold_count_params)
    contrafold_tp, contrafold_tn, contrafold_fp, contrafold_fn = final_sum(results)
    ppv = get_ppv(contrafold_tp, contrafold_fp)
    sens = get_sens(contrafold_tp, contrafold_fn)
    fpr = get_fpr(contrafold_tn, contrafold_fp)
    contrafold_ppvs.insert(0, ppv)
    contrafold_senss.insert(0, sens)
    contrafold_fprs.insert(0, fpr)
    contrafold_f1_scores.append(get_f1_score(ppv, sens))
    contrafold_mccs.append(get_mcc(contrafold_tp, contrafold_tn, contrafold_fp, contrafold_fn))
    results = pool.map(get_pos_neg_counts, centroidfold_count_params)
    centroidfold_tp, centroidfold_tn, centroidfold_fp, centroidfold_fn = final_sum(results)
    ppv = get_ppv(centroidfold_tp, centroidfold_fp)
    sens = get_sens(centroidfold_tp, centroidfold_fn)
    fpr = get_fpr(centroidfold_tn, centroidfold_fp)
    centroidfold_ppvs.insert(0, ppv)
    centroidfold_senss.insert(0, sens)
    centroidfold_fprs.insert(0, fpr)
    centroidfold_f1_scores.append(get_f1_score(ppv, sens))
    centroidfold_mccs.append(get_mcc(centroidfold_tp, centroidfold_tn, centroidfold_fp, centroidfold_fn))
    if gamma == 1.:
      results = pool.map(get_pos_neg_counts, rnafold_count_params)
      rnafold_tp, rnafold_tn, rnafold_fp, rnafold_fn = final_sum(results)
      rnafold_ppv = get_ppv(rnafold_tp, rnafold_fp)
      rnafold_sens = get_sens(rnafold_tp, rnafold_fn)
      rnafold_fpr = get_fpr(rnafold_tn, rnafold_fp)
      rnafold_f1_score = get_f1_score(rnafold_ppv, rnafold_sens)
      rnafold_mcc = get_mcc(rnafold_tp, rnafold_tn, rnafold_fp, rnafold_fn)
  line_1, = pyplot.plot(conshomfold_ppvs, conshomfold_senss, label = "ConsHomfold (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(bpp_conshomfold_ppvs, bpp_conshomfold_senss, label = "ConsHomfold (Posterior)", marker = "o", linestyle = ":")
  line_3, = pyplot.plot(turbofold_ppvs, turbofold_senss, label = "TurboFold", marker = "p", linestyle = "-")
  line_4, = pyplot.plot(centroidhomfold_ppvs, centroidhomfold_senss, label = "CentroidHomfold", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(contrafold_ppvs, contrafold_senss, label = "CONTRAfold", marker = "h", linestyle = "-")
  line_6, = pyplot.plot(centroidfold_ppvs, centroidfold_senss, label = "CentroidFold", marker = "v", linestyle = "-")
  line_7, = pyplot.plot(rnafold_ppv, rnafold_sens, label = "RNAfold", marker = "d", linestyle = "-")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6, line_7], loc = "lower left")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  line_1, = pyplot.plot(conshomfold_fprs, conshomfold_senss, label = "ConsHomfold (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(bpp_conshomfold_fprs, bpp_conshomfold_senss, label = "ConsHomfold (Posterior)", marker = "o", linestyle = ":")
  line_3, = pyplot.plot(turbofold_fprs, turbofold_senss, label = "TurboFold", marker = "p", linestyle = "-")
  line_4, = pyplot.plot(centroidhomfold_fprs, centroidhomfold_senss, label = "CentroidHomfold", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(contrafold_fprs, contrafold_senss, label = "CONTRAfold", marker = "v", linestyle = "-")
  line_6, = pyplot.plot(centroidfold_fprs, centroidfold_senss, label = "Centroidfold", marker = "h", linestyle = "-")
  line_7, = pyplot.plot(rnafold_fpr, rnafold_sens, label = "RNAfold", marker = "d", linestyle = "-")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  gammas = [i for i in range(-7, 11)]
  line_1, = pyplot.plot(gammas, conshomfold_f1_scores, label = "ConsHomfold (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, bpp_conshomfold_f1_scores, label = "ConsHomfold (Posterior)", marker = "o", linestyle = ":")
  line_3, = pyplot.plot(gammas, turbofold_f1_scores, label = "TurboFold", marker = "p", linestyle = "-")
  line_4, = pyplot.plot(gammas, centroidhomfold_f1_scores, label = "CentroidHomfold", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(gammas, contrafold_f1_scores, label = "CONTRAfold", marker = "v", linestyle = "-")
  line_6, = pyplot.plot(gammas, centroidfold_f1_scores, label = "Centroidfold", marker = "h", linestyle = "-")
  line_7, = pyplot.plot(0., rnafold_f1_score, label = "RNAfold", marker = "d", linestyle = "-")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6, line_7], loc = "lower right")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  line_1, = pyplot.plot(gammas, conshomfold_mccs, label = "ConsHomfold (Turner)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, bpp_conshomfold_mccs, label = "ConsHomfold (Posterior)", marker = "o", linestyle = ":")
  line_3, = pyplot.plot(gammas, turbofold_mccs, label = "TurboFold", marker = "p", linestyle = "-")
  line_4, = pyplot.plot(gammas, centroidhomfold_mccs, label = "CentroidHomfold", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(gammas, contrafold_mccs, label = "CONTRAfold", marker = "v", linestyle = "-")
  line_6, = pyplot.plot(gammas, centroidfold_mccs, label = "Centroidfold", marker = "h", linestyle = "-")
  line_7, = pyplot.plot(0., rnafold_mcc, label = "RNAfold", marker = "d", linestyle = "-")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("MCC")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_ss_estimation.eps", bbox_inches = "tight")

def get_pos_neg_counts(params):
  (estimated_ss_file_path, ref_sss_and_flat_sss, rna_seq_lens) = params
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

def final_sum(results):
  final_tp = final_tn = final_fp = final_fn = 0.
  for tp, tn, fp, fn in results:
    final_tp += tp
    final_tn += tn
    final_fp += fp
    final_fn += fn
  return (final_tp, final_tn, final_fp, final_fn)

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
