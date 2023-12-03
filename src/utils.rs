use lambdaworks_math::elliptic_curve::short_weierstrass::curves::bls12_381::default_types::FrElement;

pub fn vector_matrix_multiply(
    vector: &Vec<FrElement>,
    matrix: &Vec<Vec<FrElement>>,
) -> Vec<FrElement> {
    if vector.len() != matrix.len() {
        panic!("Dimensions do not match for vector-matrix multiplication");
    }

    let mut result = Vec::with_capacity(matrix[0].len());
    for col in 0..matrix[0].len() {
        let mut sum = FrElement::from_hex_unchecked("0");
        for (row, val) in vector.iter().enumerate() {
            sum = sum + (val.clone() * matrix[row][col].clone());
        }
        result.push(sum);
    }

    result
}
