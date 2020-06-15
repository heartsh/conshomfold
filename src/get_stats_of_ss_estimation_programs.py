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
  short_consfold_ppvs = []
  short_consfold_senss = []
  short_consfold_fprs = []
  short_centroidhomfold_ppvs = []
  short_centroidhomfold_senss = []
  short_centroidhomfold_fprs = []
  short_turbofold_ppvs = []
  short_turbofold_senss = []
  short_turbofold_fprs = []
  short_rnafold_ppv = short_rnafold_sens = 0.
  long_consfold_ppvs = []
  long_consfold_senss = []
  long_consfold_fprs = []
  long_centroidhomfold_ppvs = []
  long_centroidhomfold_senss = []
  long_centroidhomfold_fprs = []
  long_turbofold_ppvs = []
  long_turbofold_senss = []
  long_turbofold_fprs = []
  long_rnafold_ppv = long_rnafold_sens = 0.
  bpp_consfold_ppvs = []
  bpp_consfold_senss = []
  bpp_consfold_fprs = []
  # gammas = [2. ** i for i in range(-7, 11)]
  gammas = [2. ** i for i in range(-5, 11)]
  short_consfold_f1_score = short_centroidhomfold_f1_score = short_turbofold_f1_score = short_rnafold_f1_score = 0.
  short_consfold_mcc = short_centroidhomfold_mcc = short_turbofold_mcc = short_rnafold_mcc = 0.
  long_consfold_f1_score = long_centroidhomfold_f1_score = long_turbofold_f1_score = long_rnafold_f1_score = 0.
  long_consfold_mcc = long_centroidhomfold_mcc = long_turbofold_mcc = long_rnafold_mcc = 0.
  bpp_consfold_ss_dir_path = asset_dir_path + "/bpp_consfold"
  for data_set in ["short", "long"]:
    consfold_ss_dir_path = asset_dir_path + "/consfold_" + data_set
    centroidhomfold_ss_dir_path = asset_dir_path + "/centroidhomfold_" + data_set
    turbofold_ss_dir_path = asset_dir_path + "/turbofold_" + data_set
    rnafold_ss_dir_path = asset_dir_path + "/rnafold_" + data_set
    # rna_fam_dir_path = asset_dir_path + "/ref_sss_" + data_set
    rna_fam_dir_path = asset_dir_path + "/ref_sss_" + data_set + "_4_micro_bench"
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      consfold_tp = consfold_tn = consfold_fp = consfold_fn = 0.
      bpp_consfold_tp = bpp_consfold_tn = bpp_consfold_fp = bpp_consfold_fn = 0.
      centroidhomfold_tp = centroidhomfold_tn = centroidhomfold_fp = centroidhomfold_fn = 0.
      turbofold_tp = turbofold_tn = turbofold_fp = turbofold_fn = 0.
      rnafold_tp = rnafold_tn = rnafold_fp = rnafold_fn = 0.
      for rna_fam_file in os.listdir(rna_fam_dir_path):
        if not rna_fam_file.endswith(".fa"):
          continue
        rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
        rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
        (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
        ref_ss_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
        ref_sss_and_flat_sss = utils.get_sss_and_flat_sss(utils.get_ss_strings(ref_ss_file_path))
        consfold_estimated_ss_dir_path = os.path.join(consfold_ss_dir_path, rna_fam_name)
        if not os.path.isdir(consfold_estimated_ss_dir_path):
          continue
        bpp_consfold_estimated_ss_dir_path = os.path.join(bpp_consfold_ss_dir_path, rna_fam_name)
        if not os.path.isdir(bpp_consfold_estimated_ss_dir_path):
          continue
        centroidhomfold_estimated_ss_dir_path = os.path.join(centroidhomfold_ss_dir_path, rna_fam_name)
        if not os.path.isdir(centroidhomfold_estimated_ss_dir_path):
          continue
        turbofold_estimated_ss_dir_path = os.path.join(turbofold_ss_dir_path, rna_fam_name)
        if not os.path.isdir(turbofold_estimated_ss_dir_path):
          continue
        consfold_estimated_ss_file_path = os.path.join(consfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
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
          f1_score = 2 * ppv * sens / (ppv + sens)
          mcc = (consfold_tp * consfold_tn - consfold_fp * consfold_fn) / sqrt((consfold_tp + consfold_fp) * (consfold_tp + consfold_fn) * (consfold_tn + consfold_fp) * (consfold_tn + consfold_fn))
          if data_set == "short":
            short_consfold_f1_score = f1_score
            short_consfold_mcc = mcc
          else:
            long_consfold_f1_score = f1_score
            long_consfold_mcc = mcc
        if gamma >= 2.:
          bpp_consfold_estimated_ss_file_path = os.path.join(bpp_consfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
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
        centroidhomfold_estimated_ss_file_path = os.path.join(centroidhomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
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
          f1_score = 2 * ppv * sens / (ppv + sens)
          mcc = (centroidhomfold_tp * centroidhomfold_tn - centroidhomfold_fp * centroidhomfold_fn) / sqrt((centroidhomfold_tp + centroidhomfold_fp) * (centroidhomfold_tp + centroidhomfold_fn) * (centroidhomfold_tn + centroidhomfold_fp) * (centroidhomfold_tn + centroidhomfold_fn))
          if data_set == "short":
            short_centroidhomfold_f1_score = f1_score
            short_centroidhomfold_mcc = mcc
          else:
            long_centroidhomfold_f1_score = f1_score
            long_centroidhomfold_mcc = mcc
        turbofold_estimated_ss_file_path = os.path.join(turbofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
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
          f1_score = 2 * ppv * sens / (ppv + sens)
          mcc = (turbofold_tp * turbofold_tn - turbofold_fp * turbofold_fn) / sqrt((turbofold_tp + turbofold_fp) * (turbofold_tp + turbofold_fn) * (turbofold_tn + turbofold_fp) * (turbofold_tn + turbofold_fn))
          if data_set == "short":
            short_turbofold_f1_score = f1_score
            short_turbofold_mcc = mcc
          else:
            long_turbofold_f1_score = f1_score
            long_turbofold_mcc = mcc
        if gamma == 4.:
          rnafold_estimated_ss_file_path = os.path.join(rnafold_ss_dir_path, rna_fam_name + ".fa")
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
          ppv = rnafold_tp / (rnafold_tp + rnafold_fp)
          sens = rnafold_tp / (rnafold_tp + rnafold_fn)
          fpr = rnafold_fp / (rnafold_tn + rnafold_fp)
          f1_score = 2 * ppv * sens / (ppv + sens)
          mcc = (rnafold_tp * rnafold_tn - rnafold_fp * rnafold_fn) / sqrt((rnafold_tp + rnafold_fp) * (rnafold_tp + rnafold_fn) * (rnafold_tn + rnafold_fp) * (rnafold_tn + rnafold_fn))
          if data_set == "short":
            short_rnafold_ppv = ppv
            short_rnafold_sens = sens
            short_rnafold_fpr = fpr 
            short_rnafold_f1_score = f1_score
            short_rnafold_mcc = mcc
          else:
            long_rnafold_ppv = ppv
            long_rnafold_sens = sens
            long_rnafold_fpr = fpr 
            long_rnafold_f1_score = f1_score
            long_rnafold_mcc = mcc
      ppv = consfold_tp / (consfold_tp + consfold_fp)
      sens = consfold_tp / (consfold_tp + consfold_fn)
      fpr = consfold_fp / (consfold_tn + consfold_fp)
      if data_set == "short":
        short_consfold_ppvs.insert(0, ppv)
        short_consfold_senss.insert(0, sens)
        short_consfold_fprs.insert(0, fpr)
      else:
        long_consfold_ppvs.insert(0, ppv)
        long_consfold_senss.insert(0, sens)
        long_consfold_fprs.insert(0, fpr)
      if data_set == "short" and gamma >= 2.:
        ppv = bpp_consfold_tp / (bpp_consfold_tp + bpp_consfold_fp)
        sens = bpp_consfold_tp / (bpp_consfold_tp + bpp_consfold_fn)
        fpr = bpp_consfold_fp / (bpp_consfold_tn + bpp_consfold_fp)
        bpp_consfold_ppvs.insert(0, ppv)
        bpp_consfold_senss.insert(0, sens)
        bpp_consfold_fprs.insert(0, fpr)
      ppv = centroidhomfold_tp / (centroidhomfold_tp + centroidhomfold_fp)
      sens = centroidhomfold_tp / (centroidhomfold_tp + centroidhomfold_fn)
      fpr = centroidhomfold_fp / (centroidhomfold_tn + centroidhomfold_fp)
      if data_set == "short":
        short_centroidhomfold_ppvs.insert(0, ppv)
        short_centroidhomfold_senss.insert(0, sens)
        short_centroidhomfold_fprs.insert(0, fpr)
      else:
        long_centroidhomfold_ppvs.insert(0, ppv)
        long_centroidhomfold_senss.insert(0, sens)
        long_centroidhomfold_fprs.insert(0, fpr)
      ppv = turbofold_tp / (turbofold_tp + turbofold_fp)
      sens = turbofold_tp / (turbofold_tp + turbofold_fn)
      fpr = turbofold_fp / (turbofold_tn + turbofold_fp)
      if data_set == "short":
        short_turbofold_ppvs.insert(0, ppv)
        short_turbofold_senss.insert(0, sens)
        short_turbofold_fprs.insert(0, fpr)
      else:
        long_turbofold_ppvs.insert(0, ppv)
        long_turbofold_senss.insert(0, sens)
        long_turbofold_fprs.insert(0, fpr)
  short_consfold_ppvs = numpy.array(short_consfold_ppvs) 
  short_consfold_senss = numpy.array(short_consfold_senss)
  short_consfold_fprs = numpy.array(short_consfold_fprs)
  short_centroidhomfold_ppvs = numpy.array(short_centroidhomfold_ppvs) 
  short_centroidhomfold_senss = numpy.array(short_centroidhomfold_senss)
  short_centroidhomfold_fprs = numpy.array(short_centroidhomfold_fprs)
  short_turbofold_ppvs = numpy.array(short_turbofold_ppvs) 
  short_turbofold_senss = numpy.array(short_turbofold_senss)
  short_turbofold_fprs = numpy.array(short_turbofold_fprs)
  long_consfold_ppvs = numpy.array(long_consfold_ppvs) 
  long_consfold_senss = numpy.array(long_consfold_senss)
  long_consfold_fprs = numpy.array(long_consfold_fprs)
  long_centroidhomfold_ppvs = numpy.array(long_centroidhomfold_ppvs) 
  long_centroidhomfold_senss = numpy.array(long_centroidhomfold_senss)
  long_centroidhomfold_fprs = numpy.array(long_centroidhomfold_fprs)
  long_turbofold_ppvs = numpy.array(long_turbofold_ppvs) 
  long_turbofold_senss = numpy.array(long_turbofold_senss)
  long_turbofold_fprs = numpy.array(long_turbofold_fprs)
  bpp_consfold_ppvs = numpy.array(bpp_consfold_ppvs) 
  bpp_consfold_senss = numpy.array(bpp_consfold_senss)
  bpp_consfold_fprs = numpy.array(bpp_consfold_fprs)
  line_1, = pyplot.plot(short_consfold_ppvs, short_consfold_senss, label = "ConsFold (Short)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(short_turbofold_ppvs, short_turbofold_senss, label = "TurboFold-smp (Short)", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(short_centroidhomfold_ppvs, short_centroidhomfold_senss, label = "CentroidHomFold (Short)", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(short_rnafold_ppv, short_rnafold_sens, label = "RNAfold (Short)", marker = "d", linestyle = "-")
  line_5, = pyplot.plot(bpp_consfold_ppvs, bpp_consfold_senss, label = "ConsFold (Short & conventional)", marker = "o", linestyle = ":")
  line_6, = pyplot.plot(long_consfold_ppvs, long_consfold_senss, label = "ConsFold (Long)", marker = "o", linestyle = "--")
  line_7, = pyplot.plot(long_turbofold_ppvs, long_turbofold_senss, label = "TurboFold-smp (Long)", marker = "p", linestyle = "--")
  line_8, = pyplot.plot(long_centroidhomfold_ppvs, long_centroidhomfold_senss, label = "CentroidHomFold (Long unaligned)", marker = "s", linestyle = "--")
  line_9, = pyplot.plot(long_rnafold_ppv, long_rnafold_sens, label = "RNAfold (Long)", marker = "d", linestyle = "--")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5], loc = "upper right")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  pyplot.clf()
  line_1, = pyplot.plot(short_consfold_fprs, short_consfold_senss, label = "ConsFold (Short)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(short_turbofold_fprs, short_turbofold_senss, label = "TurboFold-smp (Short)", marker = "p", linestyle = "-")
  line_3, = pyplot.plot(short_centroidhomfold_fprs, short_centroidhomfold_senss, label = "CentroidHomFold (Short)", marker = "s", linestyle = "-")
  line_4, = pyplot.plot(short_rnafold_fpr, short_rnafold_sens, label = "RNAfold (Short)", marker = "d", linestyle = "-")
  line_5, = pyplot.plot(bpp_consfold_fprs, bpp_consfold_senss, label = "ConsFold (Short & conventional)", marker = "o", linestyle = ":")
  line_6, = pyplot.plot(long_consfold_fprs, long_consfold_senss, label = "ConsFold (Long)", marker = "o", linestyle = "--")
  line_7, = pyplot.plot(long_turbofold_fprs, long_turbofold_senss, label = "TurboFold-smp (Long)", marker = "p", linestyle = "--")
  line_8, = pyplot.plot(long_centroidhomfold_fprs, long_centroidhomfold_senss, label = "CentroidHomFold (Long)", marker = "s", linestyle = "--")
  line_9, = pyplot.plot(long_rnafold_fpr, long_rnafold_sens, label = "RNAfold (Long)", marker = "d", linestyle = "--")
  pyplot.legend(handles = [line_6, line_7, line_8, line_9], loc = "lower right")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_ss_estimation.eps", bbox_inches = "tight")
  print("ConsFold's MCC & F1 score for dataset \"short\" = %.3f & %.3f" %(short_consfold_mcc, short_consfold_f1_score))
  print("TurboFold's MCC & F1 score for dataset \"short\" = %.3f & %.3f" %(short_turbofold_mcc, short_turbofold_f1_score))
  print("CentroidHomFold's MCC & F1 score for dataset \"short\" = %.3f & %.3f" %(short_centroidhomfold_mcc, short_centroidhomfold_f1_score))
  print("RNAfold's MCC & F1 score for dataset \"short\" = %.3f & %.3f" %(short_rnafold_mcc, short_rnafold_f1_score))
  print("ConsFold's MCC & F1 score for dataset \"long\" = %.3f & %.3f" %(long_consfold_mcc, long_consfold_f1_score))
  print("TurboFold's MCC & F1 score for dataset \"long\" = %.3f & %.3f" %(long_turbofold_mcc, long_turbofold_f1_score))
  print("CentroidHomFold's MCC & F1 score for dataset \"long\" = %.3f & %.3f" %(long_centroidhomfold_mcc, long_centroidhomfold_f1_score))
  print("RNAfold's MCC & F1 score for dataset \"long\" = %.3f & %.3f" %(long_rnafold_mcc, long_rnafold_f1_score))

if __name__ == "__main__":
  main()
