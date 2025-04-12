//! Brute-force irreducible polynomial searcher.

use crate::field::Field;
use crate::poly::Poly;

/// Enumerate all monic polynomials of the given order over the given field,
/// and use them to try to divide the given polynomial.
fn try_divide<F: Field>(order: usize, poly: &mut Poly<F>, record: &mut Option<Poly<F>>) {
    fn dfs<F: Field>(
        mut divisor: Poly<F>,
        i: usize,
        order: usize,
        poly: &mut Poly<F>,
        record: &mut Option<Poly<F>>,
    ) {
        if i == order {
            if (poly.clone() % divisor.clone()).degree().is_none() {
                // Found a divisor, divide the polynomial by it.
                if record.is_none() {
                    *record = Some(divisor.clone());
                }
                *poly = poly.clone() / divisor;
            }
            return;
        }
        for j in 0..poly.field().order() {
            divisor.set(i, j);
            dfs(divisor.clone(), i + 1, order, poly, record);
        }
    }
    let divisor = {
        let mut divisor = Poly::new(poly.field().clone());
        divisor.set(order, 1);
        divisor
    };
    dfs(divisor, 0, order, poly, record);
}

/// Find an irreducible polynomial of the given order over the given field.
pub(crate) fn find_irr<F: Field>(f: F, order: usize) -> Poly<F> {
    if order == 0 {
        panic!("order must be greater than 0");
    }

    let q = f.order();
    let t = q.pow(order as u32);

    // Irreducible polynomials of order `n` are factors of `x^q^n - x`.
    // But we should subtract the factor `x` out first, since it's not we want.
    let mut poly = Poly::new(f.clone());
    poly.set(t - 1, 1);
    poly.set(0, f.neg(1));

    // Remove all divisors of orders lower than `n`.
    for i in 1..order {
        if order % i != 0 {
            continue;
        }
        let mut record = None;
        try_divide(i, &mut poly, &mut record);
    }

    // Find a divisor of order `n`.
    let mut record = None;
    try_divide(order, &mut poly, &mut record);
    debug_assert!(
        record.is_some(),
        "failed to find an irreducible polynomial of order {}",
        order,
    );
    record.unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prime::PrimeField;

    #[test]
    fn find_irr_2() {
        let f = PrimeField::new(2);
        let irr = find_irr(f.clone(), 2);

        debug_assert!(irr.degree().is_some());
        assert_eq!(irr.degree().unwrap(), 2);
        assert_eq!(irr[0], 1);
        assert_eq!(irr[2], 1);
        assert_eq!(irr[1], 1);
    }

    #[test]
    fn find_irr_3() {
        let f = PrimeField::new(2);
        let irr = find_irr(f.clone(), 3);

        debug_assert!(irr.degree().is_some());
        assert_eq!(irr.degree().unwrap(), 3);
        assert_eq!(irr[0], 1);
        assert_eq!(irr[3], 1);
        debug_assert!(
            (irr[1] == 0 && irr[2] == 1) || (irr[1] == 1 && irr[2] == 0),
            "irreducible polynomial should be of the form x^3 + x^2 + 1 or x^3 + x + 1",
        );
    }
}
