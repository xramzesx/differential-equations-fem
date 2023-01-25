// Gaussian Elimination of Quadratic Matrices
// Takes an augmented matrix as input, returns vector of results
// Wikipedia reference: augmented matrix: https://en.wikipedia.org/wiki/Augmented_matrix
// Wikipedia reference: algorithm: https://en.wikipedia.org/wiki/Gaussian_elimination
// Source : https://github.com/TheAlgorithms/Rust/blob/master/src/math/gaussian_elimination.rs

pub fn gaussian_elimination(matrix: &mut [Vec<f64>]) -> Vec<f64> {
    let size = matrix.len();
    assert_eq!(size, matrix[0].len() - 1);

    for i in 0..size - 1 {
        for j in i..size - 1 {
            echelon(matrix, i, j);
        }
    }

    for i in (1..size).rev() {
        eliminate(matrix, i);
    }

    // Disable cargo clippy warnings about needless range loops.
    // Checking the diagonal like this is simpler than any alternative.
    #[allow(clippy::needless_range_loop)]
    for i in 0..size {
        if matrix[i][i] == 0f64 {
            println!("Infinitely many solutions");
        }
    }

    let mut result: Vec<f64> = vec![0f64; size];
    for i in 0..size {
        result[i] = matrix[i][size] / matrix[i][i];
    }
    result
}

fn echelon(matrix: &mut [Vec<f64>], i: usize, j: usize) {
    let size = matrix.len();
    if matrix[i][i] == 0f64 {
    } else {
        let factor = matrix[j + 1][i] / matrix[i][i];
        (i..size + 1).for_each(|k| {
            matrix[j + 1][k] -= factor * matrix[i][k];
        });
    }
}

fn eliminate(matrix: &mut [Vec<f64>], i: usize) {
    let size = matrix.len();
    if matrix[i][i] == 0f64 {
    } else {
        for j in (1..i + 1).rev() {
            let factor = matrix[j - 1][i] / matrix[i][i];
            for k in (0..size + 1).rev() {
                matrix[j - 1][k] -= factor * matrix[i][k];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::gaussian_elimination;

    #[test]
    fn test_gauss() {
        let mut matrix: Vec<Vec<f64>> = vec![
            vec![1.5, 2.0, 1.0, -1.0, -2.0, 1.0, 1.0],
            vec![3.0, 3.0, -1.0, 16.0, 18.0, 1.0, 1.0],
            vec![1.0, 1.0, 3.0, -2.0, -6.0, 1.0, 1.0],
            vec![1.0, 1.0, 99.0, 19.0, 2.0, 1.0, 1.0],
            vec![1.0, -2.0, 16.0, 1.0, 9.0, 10.0, 1.0],
            vec![1.0, 3.0, 1.0, -5.0, 1.0, 1.0, 95.0],
        ];
        let result: Vec<f64> = vec![
            -264.0590678357175, 159.63207198892474, 
            -6.156922011998154, 35.310383017997204, 
            -18.806691278264864, 81.67838024919239
        ];
        assert_eq!(gaussian_elimination(&mut matrix), result);
    }
}
