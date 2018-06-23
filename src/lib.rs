extern crate strap;

pub use strap::*;

pub type Mea = Prob;
#[derive(Clone)]
pub struct MeaSs {
  pub bp_pos_pair_seqs_inside_pos_pairs: PosPairSeqsWithPosPairs,
  pub ea: Mea,
}
type Meas = Vec<Mea>;
type MeaMat = HashMap<PosPair, Mea, Hasher>;

impl MeaSs {
  pub fn new() -> MeaSs {
    MeaSs {
      bp_pos_pair_seqs_inside_pos_pairs: PosPairSeqsWithPosPairs::default(),
      ea: 0.,
    }
  }
}

pub fn neofold(bpp_mat: &SparseProbMat, seq_len: usize, gamma_plus_1: Prob) -> MeaSs {
  let mut mea_mat_4_bp_pos_pairs = MeaMat::default();
  let mut pos_seqs_with_poss_4_forward_bps = PosSeqsWithPoss::default();
  let inversed_gamma_plus_1 = 1. / gamma_plus_1;
  for sub_seq_len in 2 .. seq_len + 3 {
    for i in 0 .. seq_len + 3 - sub_seq_len {
      let j = i + sub_seq_len - 1;
      let pos_pair = (i, j);
      match bpp_mat.get(&pos_pair) {
        Some(&bpp) => {
          if bpp <= inversed_gamma_plus_1 {
            continue;
          }
          let meas_4_bp_pos_pair = get_meas_4_bp_pos_pair(&pos_pair, &mea_mat_4_bp_pos_pairs, &pos_seqs_with_poss_4_forward_bps);
          mea_mat_4_bp_pos_pairs.insert(pos_pair, meas_4_bp_pos_pair[j - i - 1] + gamma_plus_1 * bpp - 1.);
          let poss_exist = match pos_seqs_with_poss_4_forward_bps.get(&j) {
            Some(_) => {true},
            None => {false},
          };
          if poss_exist {
            pos_seqs_with_poss_4_forward_bps.get_mut(&j).expect("Failed to get an element of a hash map.").push(i);
          } else {
            pos_seqs_with_poss_4_forward_bps.insert(j, vec![i]);
          }
        },
        None => {},
      }
    }
  }
  let mut mea_ss = MeaSs::new();
  let pseudo_pos_pair = (0, seq_len + 1);
  let mut pos_pair_stack = vec![pseudo_pos_pair];
  while pos_pair_stack.len() > 0 {
    let pos_pair_1 = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let meas_4_bp_pos_pair = get_meas_4_bp_pos_pair(&pos_pair_1, &mea_mat_4_bp_pos_pairs, &pos_seqs_with_poss_4_forward_bps);
    let (i, j) = pos_pair_1;
    let mea = meas_4_bp_pos_pair[j - i - 1];
    if mea == 0. {continue;}
    let mut n = j - 1;
    while meas_4_bp_pos_pair[n - i] > 0. {
      let mea = meas_4_bp_pos_pair[n - i];
      if mea == meas_4_bp_pos_pair[n - i - 1] {
        n = n - 1;
      } else {
        match pos_seqs_with_poss_4_forward_bps.get(&n) {
          Some(poss) => {
            for &m in poss {
              if m <= i {continue;}
              let pos_pair_2 = (m, n);
              if mea == meas_4_bp_pos_pair[m - i - 1] + mea_mat_4_bp_pos_pairs[&pos_pair_2] {
                let bp_pos_pairs_exist = match mea_ss.bp_pos_pair_seqs_inside_pos_pairs.get(&pos_pair_1) {
                  Some(_) => {true},
                  None => {false},
                };
                if bp_pos_pairs_exist {
                  mea_ss.bp_pos_pair_seqs_inside_pos_pairs.get_mut(&pos_pair_1).expect("Failed to get an element of a hash map.").push(pos_pair_2);
                } else {
                  mea_ss.bp_pos_pair_seqs_inside_pos_pairs.insert(pos_pair_1, vec![pos_pair_2]);
                }
                pos_pair_stack.push(pos_pair_2);
                n = m - 1;
                break;
              }
            }
          },
          None => {},
        }
      }
    }
  }
  mea_ss.ea = mea_mat_4_bp_pos_pairs[&pseudo_pos_pair];
  mea_ss
}

fn get_meas_4_bp_pos_pair(pos_pair: &PosPair, mea_mat_4_bp_pos_pairs: &MeaMat, pos_seqs_with_poss_4_forward_bps: &PosSeqsWithPoss) -> Meas {
  let (i, j) = *pos_pair;
  let sub_seq_len = j - i + 1;
  let mut meas_4_bp_pos_pair = vec![0.; sub_seq_len - 1];
  for n in i + 2 .. j {
    let mut mea = 0.;
    match pos_seqs_with_poss_4_forward_bps.get(&n) {
      Some(poss) => {
        for &m in poss {
          if m <= i {continue;}
          let ea = meas_4_bp_pos_pair[m - i - 1] + mea_mat_4_bp_pos_pairs[&(m, n)];
          if ea > mea {
            mea = ea;
          }
        }
      },
      None => {},
    };
    let ea = meas_4_bp_pos_pair[n - i - 1];
    if ea > mea {
      mea = ea;
    }
    meas_4_bp_pos_pair[n - i] = mea;
  }
  meas_4_bp_pos_pair
}
