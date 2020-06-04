extern crate phylofold;
extern crate bio;
extern crate num_cpus;

use phylofold::*;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::BufWriter;
use std::fs::File;
use std::fs::create_dir;

const DEFAULT_MIN_POW_OF_2: i32 = -7;
const DEFAULT_MAX_POW_OF_2: i32 = 10;
const GAMMA_4_BENCH: Prob = 1.;

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file that contains RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("", "min_pow_of_2", &format!("A minimum power of 2 to calculate a gamma parameter (Uses {} by default)", DEFAULT_MIN_POW_OF_2), "FLOAT");
  opts.optopt("", "max_pow_of_2", &format!("A maximum power of 2 to calculate a gamma parameter (Uses {} by default)", DEFAULT_MAX_POW_OF_2), "FLOAT");
  opts.optflag("u", "uses_bpp_score", "Uses base-pairing probabilities as scores of secondary structures (Not recommended due to poor accuracy)");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("b", "takes_bench", &format!("Compute for only gamma = {} to measure running time", GAMMA_4_BENCH));
  opts.optflag("a", "produces_access_probs", &format!("Also compute accessible probabilities"));
  opts.optflag("p", "outputs_probs", &format!("Output probabilities"));
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
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
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
  let min_pow_of_2 = if matches.opt_present("min_pow_of_2") {
    matches.opt_str("min_pow_of_2").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_POW_OF_2
  };
  let max_pow_of_2 = if matches.opt_present("max_pow_of_2") {
    matches.opt_str("max_pow_of_2").unwrap().parse().unwrap()
  } else {
    DEFAULT_MAX_POW_OF_2
  };
  let uses_bpps = matches.opt_present("u");
  let takes_bench = matches.opt_present("b");
  let produces_access_probs = matches.opt_present("a");
  let outputs_probs = matches.opt_present("p");
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
  let prob_mat_sets = phyloprob(&mut thread_pool, &fasta_records, min_bpp, offset_4_max_gap_num, uses_bpps, produces_access_probs);
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  if takes_bench {
    let num_of_fasta_records = fasta_records.len();
    let mut mea_sss = vec![MeaSs::new(); num_of_fasta_records];
    let output_file_path = output_dir_path.join(&format!("gamma={}.dat", GAMMA_4_BENCH));
    thread_pool.scoped(|scope| {
      for ((rna_id, fasta_record), mea_ss) in fasta_records.iter().enumerate().zip(mea_sss.iter_mut()) {
        let ref prob_mats = prob_mat_sets[rna_id];
        scope.execute(move || {
          *mea_ss = phylofold(&prob_mats.bpp_mat, &prob_mats.upp_mat, fasta_record.seq.len(), GAMMA_4_BENCH);
        });
      }
    });
    let mut buf = String::new();
    let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
    for (rna_id, mea_ss) in mea_sss.iter().enumerate() {
      let buf_4_rna_id = format!(">{}\n", rna_id) + &unsafe {String::from_utf8_unchecked(get_mea_ss_str(&mea_ss, fasta_records[rna_id].seq.len()))} + if rna_id < num_of_fasta_records - 1 {"\n"} else {""};
      buf.push_str(&buf_4_rna_id);
    }
    let _ = writer_2_output_file.write_all(buf.as_bytes());
  } else {
    thread_pool.scoped(|scope| {
      for pow_of_2 in min_pow_of_2 .. max_pow_of_2 + 1 {
        let gamma = (2. as Prob).powi(pow_of_2);
        let ref ref_2_prob_mat_sets = prob_mat_sets;
        let ref ref_2_fasta_records = fasta_records;
        let output_file_path = output_dir_path.join(&format!("gamma={}.dat", gamma));
        scope.execute(move || {
          compute_and_write_mea_sss(ref_2_prob_mat_sets, ref_2_fasta_records, gamma, &output_file_path);
        });
      }
    });
  }
  if outputs_probs {
    write_prob_mat_sets(&output_dir_path, &prob_mat_sets, produces_access_probs);
  }
}

fn compute_and_write_mea_sss(prob_mat_sets: &ProbMatSets, fasta_records: &FastaRecords, gamma: Prob, output_file_path: &Path) {
  let num_of_fasta_records = fasta_records.len();
  let mut buf = String::new();
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  for (rna_id, fasta_record) in fasta_records.iter().enumerate() {
    let ref prob_mats = prob_mat_sets[rna_id];
    let mea_ss = phylofold(&prob_mats.bpp_mat, &prob_mats.upp_mat, fasta_record.seq.len(), gamma);
    let buf_4_rna_id = format!(">{}\n", rna_id) + &unsafe {String::from_utf8_unchecked(get_mea_ss_str(&mea_ss, fasta_records[rna_id].seq.len()))} + if rna_id < num_of_fasta_records - 1 {"\n"} else {""};
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
