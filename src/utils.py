import os
import matplotlib
from matplotlib import pylab
import numpy
import subprocess
from Bio import SeqIO

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

def get_ss_and_flat_ss(ss_string):
  ss = {}
  flat_ss = {}
  stack = []
  for (i, char) in enumerate(ss_string):
    if char == "(":
      stack.append(i)
    elif char == ")":
      pos = stack.pop()
      ss[(pos, i)] = True
      flat_ss[pos] = True
      flat_ss[i] = True
  return ss, flat_ss

def get_ss_strings(ss_file_path):
  ss_strings = [rec.seq for rec in SeqIO.parse(ss_file_path, "fasta")]
  return ss_strings

def get_sss_and_flat_sss(ss_strings):
  return list(map(get_ss_and_flat_ss, ss_strings))

def run_command(command):
  subproc = subprocess.Popen(
    command,
    stdout = subprocess.PIPE,
    shell = True
  )
  (output, error) = subproc.communicate()
  returned_code = subproc.wait()
  return (output, error, returned_code)
