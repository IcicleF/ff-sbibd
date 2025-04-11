//! Symmetric Balanced Incomplete Block Design (SBIBD) finder.

use std::collections::HashSet;

use crate::extfield::ExtField;
use crate::field::Field;
use crate::irr::find_irr;
use crate::poly::Poly;

/// Find a perfect difference set of order q.
pub fn find_diffset(q: usize) -> Vec<usize> {
    // Construct (F_q)^3.
    let f = ExtField::over_prime(q);
    let irr = find_irr(f.clone(), 3);
    let f3 = unsafe { ExtField::with_irr(f.clone(), q * q * q, irr) };
    let alpha = f3.primitive();

    // Find a line in the finite projective plane of order q.
    // We use z = 0 as this line.
    let line = {
        let mut line = HashSet::new();
        for j in 1..q {
            line.insert(f3.poly2elem(Poly::with_coeffs(f.clone(), [0, j, 0])));
            line.insert(f3.poly2elem(Poly::with_coeffs(f.clone(), [j, 0, 0])));
            for i in 1..q {
                line.insert(f3.poly2elem(Poly::with_coeffs(f.clone(), [j, f.mul(i, j), 0])));
            }
        }
        line
    };
    assert!(line.len() == (q + 1) * (q - 1));

    let mut x = 1;
    let mut ret = Vec::new();
    for i in 0..(q * q + q + 1) {
        // Determine whether the line contains a point equivalent to x.
        if line.contains(&x) {
            ret.push(i);
        }
        x = f3.mul(x, alpha);
    }

    assert!(
        ret.len() == q + 1,
        "ret.len() = {} for q = {}",
        ret.len(),
        q
    );
    ret
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_2() {
        let ds = find_diffset(2);
        for i in 0..7 {
            let mut shift = ds.iter().map(|x| (x + i) % 7).collect::<Vec<_>>();
            shift.sort();
            if shift == vec![0, 1, 3] || shift == vec![0, 1, 5] {
                return;
            }
        }
        panic!("failed to find a perfect difference set");
    }

    #[test]
    fn orders() {
        for i in [3, 4, 5, 6] {
            find_diffset(i - 1);
        }
    }
}
