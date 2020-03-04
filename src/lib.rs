extern crate phyloprob;

pub use phyloprob::*;

pub type Mea = Prob;
#[derive(Clone)]
pub struct MeaSs {
  pub bp_pos_pairs: PosPairs,
  pub ea: Mea,
}
pub type Meas = Vec<Mea>;
pub type MeaMat = FxHashMap<PosPair, Mea>;
pub type Poss = Vec<Pos>;
pub type PosSeqsWithPoss = FxHashMap<Pos, Poss>;
pub type PosPairs = Vec<PosPair>;
pub type PosPairSeqsWithPosPairs = FxHashMap<PosPair, PosPairs>;
pub type MeaSsChar = u8;
pub type MeaSsStr = Vec<MeaSsChar>;

impl MeaSs {
  pub fn new() -> MeaSs {
    MeaSs {
      bp_pos_pairs: PosPairs::new(),
      ea: 0.,
    }
  }
}

pub const UNPAIRING_BASE: MeaSsChar = '.' as MeaSsChar;
pub const BASE_PAIRING_LEFT_BASE: MeaSsChar = '(' as MeaSsChar;
pub const BASE_PAIRING_RIGHT_BASE: MeaSsChar = ')' as MeaSsChar;

pub fn phylofold(bpp_mat: &SparseProbMat, upp_mat: &Probs, seq_len: usize, gamma: Prob) -> MeaSs {
  let mut mea_mat = vec![vec![0.; seq_len]; seq_len];
  let seq_len = seq_len as Pos;
  for sub_seq_len in 1 .. seq_len + 1 {
    for i in 0 .. seq_len + 1 - sub_seq_len {
      let j = i + sub_seq_len - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      if i == j {
        mea_mat[long_i][long_j] = upp_mat[long_i];
        continue;
      }
      let mut mea = mea_mat[long_i + 1][long_j] + upp_mat[long_i];
      let ea = mea_mat[long_i][long_j - 1] + upp_mat[long_j];
      if ea > mea {
        mea = ea;
      }
      let pos_pair = (i, j);
      if bpp_mat.contains_key(&pos_pair) {
        let ea = mea_mat[long_i + 1][long_j - 1] + gamma * bpp_mat[&pos_pair];
        if ea > mea {
          mea = ea;
        }
      }
      for k in long_i .. long_j {
        let ea = mea_mat[long_i][k] + mea_mat[k + 1][long_j];
        if ea > mea {
          mea = ea;
        }
      }
      mea_mat[long_i][long_j] = mea;
    }
  }
  let mut mea_ss = MeaSs::new();
  let mut pos_pair_stack = vec![(0, seq_len - 1)];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j) = pos_pair;
    if j <= i {continue;}
    let (long_i, long_j) = (i as usize, j as usize);
    let mea = mea_mat[long_i][long_j];
    if mea == mea_mat[long_i + 1][long_j] + upp_mat[long_i] {
      pos_pair_stack.push((i + 1, j));
    } else if mea == mea_mat[long_i][long_j - 1] + upp_mat[long_j] {
      pos_pair_stack.push((i, j - 1));
    } else if bpp_mat.contains_key(&pos_pair) && mea == mea_mat[long_i + 1][long_j - 1] + gamma * bpp_mat[&pos_pair] {
      pos_pair_stack.push((i + 1, j - 1));
      mea_ss.bp_pos_pairs.push(pos_pair);
    } else {
      for k in i .. j {
        let long_k = k as usize;
        if mea == mea_mat[long_i][long_k] + mea_mat[long_k + 1][long_j] {
          pos_pair_stack.push((i, k));
          pos_pair_stack.push((k + 1, j));
          break;
        }
      }
    }
  }
  mea_ss.ea = mea_mat[0][seq_len as usize - 1];
  mea_ss
}
