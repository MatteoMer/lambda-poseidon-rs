mod utils;

use crate::utils::vector_matrix_multiply;
use lambdaworks_math::elliptic_curve::short_weierstrass::curves::bls12_381::default_types::FrElement as FE;

pub type FrElement = FE;

pub struct Constants {
    t: usize,
    alpha: u32,
    m: Vec<Vec<FrElement>>,
    c: Vec<FrElement>,
    r_f: u32,
    r_p: u32,
}

pub struct Poseidon {
    constants: Constants,
    round_ct: usize,
}

impl Poseidon {
    pub fn new(constants: Constants) -> Self {
        Poseidon {
            constants,
            round_ct: 0,
        }
    }

    fn full_round(&mut self, input_vec: &Vec<FrElement>, nb_round: u32) -> Vec<FrElement> {
        let mut round_vec = input_vec.clone();

        for _i in 0..nb_round {
            //s-box
            round_vec = self.s_box(&round_vec);

            //mds matrix
            round_vec = vector_matrix_multiply(&round_vec, &self.constants.m);

            //round constant
            round_vec = self.round_constant(&round_vec);
        }
        round_vec
    }

    fn partial_round(&mut self, input_vec: &Vec<FrElement>, nb_round: u32) -> Vec<FrElement> {
        let mut round_vec = input_vec.clone();

        for _i in 0..nb_round {
            //s-box
            let first_elem = self.s_box(&vec![round_vec[0].clone()]);
            round_vec[0] = first_elem[0].clone();

            //mds matrix
            round_vec = vector_matrix_multiply(&round_vec, &self.constants.m);

            //round constant
            round_vec = self.round_constant(&round_vec);
        }
        round_vec
    }

    fn s_box(&self, input: &Vec<FrElement>) -> Vec<FrElement> {
        input.iter().map(|x| x.pow(self.constants.alpha)).collect()
    }

    fn round_constant(&mut self, input: &Vec<FrElement>) -> Vec<FrElement> {
        let mut result = Vec::with_capacity(input.len());
        for elem in input {
            let sum = elem + &self.constants.c[self.round_ct];
            self.round_ct += 1;
            result.push(sum);
        }
        result
    }

    pub fn hash(&mut self, input_vec: Vec<FrElement>) -> Vec<FrElement> {
        if input_vec.len() != self.constants.t {
            panic!("The input vector doesn't match the parameters.");
        }
        // First full rounds
        let mut result_vec = self.full_round(&input_vec, self.constants.r_f / 2);

        // Partial rounds
        result_vec = self.partial_round(&result_vec, self.constants.r_p);

        // Final full rounds
        result_vec = self.full_round(&result_vec, self.constants.r_f / 2);

        result_vec
    }
}

#[cfg(test)]
mod tests;
