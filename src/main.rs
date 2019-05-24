extern crate neofold;
extern crate getopts;
extern crate scoped_threadpool;
extern crate bio;
extern crate itertools;
extern crate num_cpus;

use neofold::*;
use getopts::Options;
use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::File;
use itertools::multizip;

type Arg = String;
type NumOfThreads = u32;
type FastaId = String;
type FastaRecord = (FastaId, Seq, usize);
type FastaRecords = Vec<FastaRecord>;
type Strings = Vec<String>;
type MeaSsChar = u8;
type MeaSsStr = Vec<MeaSsChar>;

const DEFAULT_GAMMA: Prob = 1.;
const NOT_BASE_PAIRING_BASE: MeaSsChar = '.' as MeaSsChar;
const BASE_PAIRING_LEFT_BASE: MeaSsChar = '(' as MeaSsChar;
const BASE_PAIRING_RIGHT_BASE: MeaSsChar = ')' as MeaSsChar;
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("f", "input_fasta_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("p", "input_base_pair_prob_matrix_file_path", "The path to an input file containing base-pairing probability matrices", "STR");
  opts.reqopt("q", "input_unpair_prob_matrix_file_path", "The path to an input file containing unpairing probability matrices", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output file which will contain estimated secondary structures", "STR");
  opts.optopt("", "gamma", &format!("An MEA gamma (Uses {} by default)", DEFAULT_GAMMA), "FLOAT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_fasta_file_path = opts.opt_str("f").expect("Failed to get the path to an input FASTA file containing RNA sequences from command arguments.");
  let input_fasta_file_path = Path::new(&input_fasta_file_path);
  let input_bpp_mat_file_path = opts.opt_str("p").expect("Failed to get the path to an input file containing base-pairing probability matrices from command arguments.");
  let input_bpp_mat_file_path = Path::new(&input_bpp_mat_file_path);
  let input_upp_mat_file_path = opts.opt_str("q").expect("Failed to get the path to an input file containing unpairing probability matrices from command arguments.");
  let input_upp_mat_file_path = Path::new(&input_upp_mat_file_path);
  let output_file_path = opts.opt_str("o").expect("Failed to get the path to an output file which will contain estimated secondary structures from command arguments.");
  let output_file_path = Path::new(&output_file_path);
  let gamma = if opts.opt_present("gamma") {
    opts.opt_str("gamma").expect("Failed to get an MEA gamma from command arguments.").parse().expect("Failed to parse an MEA gamma.")
  } else {
    DEFAULT_GAMMA
  };
  let num_of_threads = if opts.opt_present("t") {
    opts.opt_str("t").expect("Failed to get the number of threads in multithreading from command arguments.").parse().expect("Failed to parse the number of threads in multithreading.")
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_fasta_file_path)).expect("Failed to set a FASTA file reader.");
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.expect("Failed to read a FASTA record.");
    let seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect::<Seq>()};
    let seq_len = seq.len();
    fasta_records.push((String::from(fasta_record.id()), seq, seq_len));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut bpp_mats = vec![SparseProbMat::default(); num_of_fasta_records];
  let mut reader_2_input_bpp_mat_file = BufReader::new(File::open(input_bpp_mat_file_path).expect("Failed to read an input file."));
  let mut buf_4_reader_2_input_bpp_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_bpp_mat_file.read_until(b'>', &mut buf_4_reader_2_input_bpp_mat_file);
  }
  for (i, vec) in reader_2_input_bpp_mat_file.split(b'>').enumerate() {
    if i == num_of_fasta_records {continue;}
    let vec = vec.expect("Failed to read an input file.");
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id = substrings[0].parse::<RnaId>().expect("Failed to parse an RNA ID.");
    let seq_len = fasta_records[rna_id].2;
    bpp_mats[rna_id].insert((0, seq_len + 1), 1.);
    for subsubstring in &substrings[1 ..] {
      let subsubsubstrings = subsubstring.split(",").collect::<Vec<&str>>();
      bpp_mats[rna_id].insert((
        subsubsubstrings[0].parse::<Pos>().expect("Failed to parse an index.") + 1,
        subsubsubstrings[1].parse::<Pos>().expect("Failed to parse an index.") + 1),
        subsubsubstrings[2].parse().expect("Failed to parse a base-pairing probability."),
      );
    }
  }
  let mut upp_mats = vec![Probs::new(); num_of_fasta_records];
  let mut reader_2_input_upp_mat_file = BufReader::new(File::open(input_upp_mat_file_path).expect("Failed to read an input file."));
  let mut buf_4_reader_2_input_upp_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_upp_mat_file.read_until(b'>', &mut buf_4_reader_2_input_upp_mat_file);
  }
  for (i, vec) in reader_2_input_upp_mat_file.split(b'>').enumerate() {
    if i == num_of_fasta_records {continue;}
    let vec = vec.expect("Failed to read an input file.");
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id = substrings[0].parse::<RnaId>().expect("Failed to parse an RNA ID.");
    for subsubstring in &substrings[1 ..] {
      let subsubsubstrings = subsubstring.split(",").collect::<Vec<&str>>();
      upp_mats[rna_id].push(subsubsubstrings[1].parse().expect("Failed to parse a base-pairing probability."));
    }
    upp_mats[rna_id].insert(0, 0.);
    upp_mats[rna_id].push(0.);
  }
  let mut mea_sss = vec![MeaSs::new(); num_of_fasta_records];
  let mut thread_pool = Pool::new(num_of_threads);
  thread_pool.scoped(|scope| {
    for (mea_ss, bpp_mat, upp_mat, fasta_record) in multizip((mea_sss.iter_mut(), bpp_mats.iter(), upp_mats.iter(), fasta_records.iter())) {
      scope.execute(move || {
        *mea_ss = neofold(bpp_mat, upp_mat, fasta_record.2, gamma);
      });
    }
  });
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_output_file = format!("; The version {} of the NeoFold program.\n; The path to the input FASTA file to compute the secondary structures (= SSs) in this file = \"{}\".\n; The path to the input base-pairing matrix file to compute these structures = \"{}\".\n; The values of the parameters used to compute these structures are as follows.\n; \"gamma\" = {}, \"num_of_threads\" = {}.\n; Each row beginning with \">\" is with a pair of the ID of an RNA sequence and expected accuracy of the maximum-expected-accuracy SS computed from this sequence. The row next to this row is with this SS.", VERSION, input_fasta_file_path.display(), input_bpp_mat_file_path.display(), gamma, num_of_threads);
  for (rna_id, mea_ss) in mea_sss.iter().enumerate() {
    let buf_4_rna_id = format!("\n\n>{},{}\n", rna_id, mea_ss.ea) + &unsafe {String::from_utf8_unchecked(get_mea_ss_str(mea_ss, fasta_records[rna_id].2))};
    buf_4_writer_2_output_file.push_str(&buf_4_rna_id);
  }
  let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes());
}

#[inline]
fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

#[inline]
fn get_mea_ss_str(mea_ss: &MeaSs, seq_len: usize) -> MeaSsStr {
  let mut mea_ss_str = vec![NOT_BASE_PAIRING_BASE; seq_len];
  let pseudo_pos_pair = (0, seq_len + 1);
  let mut pos_pair_stack = vec![pseudo_pos_pair];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j) = pos_pair;
    if pos_pair != pseudo_pos_pair {
      mea_ss_str[i - 1] = BASE_PAIRING_LEFT_BASE;
      mea_ss_str[j - 1] = BASE_PAIRING_RIGHT_BASE;
    }
    match mea_ss.bp_pos_pair_seqs_inside_pos_pairs.get(&pos_pair) {
      Some(pos_pairs) => {
        for pos_pair in pos_pairs {
          pos_pair_stack.push(*pos_pair);
        }
      }, None => {},
    }
  }
  mea_ss_str
}
