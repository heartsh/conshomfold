extern crate phylofold;
extern crate bio;
extern crate num_cpus;

use phylofold::*;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::BufWriter;
use std::fs::File;

const DEFAULT_GAMMA: Prob = 1.;

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file that contains RNA sequences", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output file", "STR");
  opts.optopt("", "opening_gap_penalty", &format!("An opening-gap penalty (Uses {} by default)", DEFAULT_OPENING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "extending_gap_penalty", &format!("An extending-gap penalty (Uses {} by default)", DEFAULT_EXTENDING_GAP_PENALTY), "FLOAT");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("", "gamma", &format!("An MEA gamma (Uses {} by default)", DEFAULT_GAMMA), "FLOAT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let output_file_path = matches.opt_str("o").unwrap();
  let output_file_path = Path::new(&output_file_path);
  let opening_gap_penalty = if matches.opt_present("opening_gap_penalty") {
    matches.opt_str("opening_gap_penalty").unwrap().parse().unwrap()
  } else {
    DEFAULT_OPENING_GAP_PENALTY
  };
  let extending_gap_penalty = if matches.opt_present("extending_gap_penalty") {
    matches.opt_str("extending_gap_penalty").unwrap().parse().unwrap()
  } else {
    DEFAULT_EXTENDING_GAP_PENALTY
  };
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches.opt_str("min_base_pair_prob").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_BPP
  };
  let offset_4_max_gap_num = if matches.opt_present("offset_4_max_gap_num") {
    matches.opt_str("offset_4_max_gap_num").unwrap().parse().unwrap()
  } else {
    DEFAULT_OFFSET_4_MAX_GAP_NUM
  };
  let gamma = if matches.opt_present("gamma") {
    matches.opt_str("gamma").unwrap().parse().unwrap()
  } else {
    DEFAULT_GAMMA
  };
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut thread_pool = Pool::new(num_of_threads);
  let (bpp_mats, upp_mats) = phyloprob(&mut thread_pool, &fasta_records, opening_gap_penalty, extending_gap_penalty, min_bpp, offset_4_max_gap_num);
  let num_of_fasta_records = fasta_records.len();
  let mut mea_sss = vec![MeaSs::new(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (mea_ss, bpp_mat, upp_mat, fasta_record) in multizip((mea_sss.iter_mut(), bpp_mats.iter(), upp_mats.iter(), fasta_records.iter())) {
      scope.execute(move || {
        *mea_ss = phylofold(bpp_mat, upp_mat, fasta_record.seq.len(), gamma);
      });
    }
  });
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf = String::new();
  for (rna_id, mea_ss) in mea_sss.iter().enumerate() {
    let buf_4_rna_id = format!(">{}\n", rna_id) + &unsafe {String::from_utf8_unchecked(get_mea_ss_str(mea_ss, fasta_records[rna_id].seq.len()))} + if rna_id < num_of_fasta_records - 1 {"\n"} else {""};
    buf.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_output_file.write_all(buf.as_bytes());
}

fn get_mea_ss_str(mea_ss: &MeaSs, seq_len: usize) -> MeaSsStr {
  let mut mea_ss_str = vec![UNPAIRING_BASE; seq_len - 2];
  for &(i, j) in &mea_ss.bp_pos_pairs {
    mea_ss_str[i as usize - 1] = BASE_PAIRING_LEFT_BASE;
    mea_ss_str[j as usize - 1] = BASE_PAIRING_RIGHT_BASE;
  }
  mea_ss_str
}
