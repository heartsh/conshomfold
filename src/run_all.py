#! /usr/bin/env python

import sample_rnas
import get_trna_prob_dists
import get_comm_structs
import compile_rna_fams
import run_ss_estimation_programs
import get_stats_of_ss_estimation_programs

def main():
  # sample_rnas.main()
  get_trna_prob_dists.main()
  get_comm_structs.main()
  # compile_rna_fams.main()
  run_ss_estimation_programs.main()
  get_stats_of_ss_estimation_programs.main()

if __name__ == "__main__":
  main()
