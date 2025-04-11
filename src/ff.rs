//! General finite field arithmetic.

use crate::field::Field;
use crate::irr::find_irr;
use crate::poly::Poly;
use crate::prime::PrimeField;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Order {
    q: usize,
    p: usize,
    k: usize,
}

impl Order {
    /// Take a prime power q = p^k and make it a finite field order.
    fn new(q: usize) -> Self {
        let mut p = 2;
        let mut q2 = q;

        while q2 > 1 {
            if q2 % p == 0 {
                break;
            }
            p += 1;
        }

        let mut k = 0;
        while q2 > 1 {
            if q2 % p == 0 {
                k += 1;
                q2 /= p;
            }
        }
        if q2 != 1 {
            panic!("{} is not a prime power", q);
        }
        Order { q, p, k }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteField {
    /// The order of the field, q = p^k.
    order: Order,

    /// The prime field of order p.
    pf: PrimeField,

    /// The irreducible polynomial of degree k that generates the field.
    irr: Poly<PrimeField>,
}

impl FiniteField {
    /// Construct a finite field of order q = p^k with the given irreducible polynomial.
    ///
    /// # Safety
    ///
    /// `irr` must be an irreducible polynomial of degree k over the prime field of order p.
    pub unsafe fn with_irr(q: usize, irr: Poly<PrimeField>) -> Self {
        let order = Order::new(q);
        let pf = PrimeField::new(order.p);
        assert_eq!(
            pf,
            *irr.field(),
            "irreducible polynomial must be over the prime field of order p = {} (found {})",
            order.p,
            irr.field().order()
        );
        assert!(
            irr.degree().unwrap_or(0) == order.k,
            "irreducible polynomial must be of degree k = {} (found {})",
            order.k,
            match irr.degree() {
                Some(d) => format!("{}", d),
                None => "-inf".to_owned(),
            }
        );
        assert!(irr[order.k] == 1, "irreducible polynomial must be monic");
        FiniteField { order, pf, irr }
    }

    /// Construct a finite field of order q = p^k.
    pub fn new(q: usize) -> Self {
        let order = Order::new(q);
        let pf = PrimeField::new(order.p);
        let irr = find_irr(pf, order.k);
        FiniteField { order, pf, irr }
    }

    /// Get the irreducible polynomial of the field.
    pub fn irr(&self) -> &Poly<PrimeField> {
        &self.irr
    }

    /// Regard an element as a vector over the prime field.
    fn int2poly(&self, a: usize) -> Poly<PrimeField> {
        assert!(a < self.order.q, "element out of range");

        let mut v = vec![0; self.order.k];
        let mut a = a;
        for i in 0..self.order.k {
            v[i] = a % self.order.p;
            a /= self.order.p;
        }
        Poly::with_coeffs(self.pf, v)
    }

    /// Regard a vector over the prime field as an element.
    fn poly2int(&self, a: Poly<PrimeField>) -> usize {
        assert!(
            a.degree().unwrap_or(0) < self.order.k,
            "element out of range"
        );

        let mut n = 0;
        for i in (0..self.order.k).rev() {
            n = n * self.order.p + a[i];
        }
        n
    }
}

impl Field for FiniteField {
    fn order(&self) -> usize {
        self.order.q
    }

    fn add(&self, a: usize, b: usize) -> usize {
        assert!(a < self.order.q && b < self.order.q, "element out of range");
        let va = self.int2poly(a);
        let vb = self.int2poly(b);
        self.poly2int(va + vb)
    }

    fn sub(&self, a: usize, b: usize) -> usize {
        assert!(a < self.order.q && b < self.order.q, "element out of range");
        let va = self.int2poly(a);
        let vb = self.int2poly(b);
        self.poly2int(va - vb)
    }

    fn mul(&self, a: usize, b: usize) -> usize {
        assert!(a < self.order.q && b < self.order.q, "element out of range");
        let va = self.int2poly(a);
        let vb = self.int2poly(b);
        self.poly2int((va * vb) % self.irr.clone())
    }

    fn inv(&self, a: usize) -> usize {
        assert!(a < self.order.q, "element out of range");
        assert!(a != 0, "division by zero");
        for i in 1..self.order.q {
            if self.mul(a, i) == 1 {
                return i;
            }
        }
        panic!("no inverse found");
    }

    fn primitive(&self) -> usize {
        let mut irr = self.irr.clone();
        irr.set(self.order.k, 0);
        assert!(irr.degree().unwrap() == self.order.k - 1);
        self.poly2int(irr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn custom_irr_4() {
        let pf = PrimeField::new(2);

        // Irreducible polynomial: x^2 + x + 1
        let irr = Poly::with_coeffs(pf, [1, 1, 1]);
        let ff = unsafe { FiniteField::with_irr(4, irr) };

        let add_table = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]];
        let mul_table = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2]];

        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(
                    ff.add(i, j),
                    add_table[i][j],
                    "expect {} + {} = {}, got {}",
                    i,
                    j,
                    add_table[i][j],
                    ff.add(i, j)
                );
                assert_eq!(
                    ff.mul(i, j),
                    mul_table[i][j],
                    "expect {} * {} = {}, got {}",
                    i,
                    j,
                    mul_table[i][j],
                    ff.mul(i, j)
                );
            }
        }

        for i in 0..4 {
            for j in 1..4 {
                assert_eq!(i, ff.mul(ff.div(i, j), j));
            }
        }

        assert_eq!(ff.inv(1), 1);
        assert_eq!(ff.inv(2), 3);
        assert_eq!(ff.inv(3), 2);
    }

    #[test]
    fn nine() {
        let ff = FiniteField::new(9);
        assert_eq!(ff.order(), 9);

        // (F, +) is Abelian group
        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_eq!(
                        ff.add(i, ff.add(j, k)),
                        ff.add(ff.add(i, j), k),
                        "associativity failed for {} + {} + {}",
                        i,
                        j,
                        k
                    );
                }
            }
        }
        for i in 0..9 {
            assert_eq!(ff.add(i, 0), i, "identity failed for {} + 0", i);
            assert_eq!(
                ff.add(i, ff.neg(i)),
                0,
                "inverse failed for {} + {}",
                i,
                ff.neg(i)
            );
        }
        for i in 0..9 {
            for j in 0..9 {
                assert_eq!(
                    ff.add(i, j),
                    ff.add(j, i),
                    "commutativity failed for {} + {}",
                    i,
                    j
                );
            }
        }

        // (F, *) is Abelian monoid
        for i in 1..9 {
            for j in 1..9 {
                for k in 1..9 {
                    assert_eq!(
                        ff.mul(i, ff.mul(j, k)),
                        ff.mul(ff.mul(i, j), k),
                        "associativity failed for {} * {} * {}",
                        i,
                        j,
                        k
                    );
                }
            }
        }
        for i in 0..9 {
            assert_eq!(ff.mul(i, 1), i, "identity failed for {} * 1", i);
            if i != 0 {
                assert_eq!(
                    ff.mul(i, ff.inv(i)),
                    1,
                    "inverse failed for {} * {}",
                    i,
                    ff.inv(i)
                );
            }
        }
        for i in 1..9 {
            for j in 1..9 {
                assert_eq!(
                    ff.mul(i, j),
                    ff.mul(j, i),
                    "commutativity failed for {} * {}",
                    i,
                    j
                );
            }
        }

        // Distributivity for `*` over `+`
        for i in 0..9 {
            for j in 0..9 {
                for k in 0..9 {
                    assert_eq!(
                        ff.mul(i, ff.add(j, k)),
                        ff.add(ff.mul(i, j), ff.mul(i, k)),
                        "distributivity failed for {} * ({} + {})",
                        i,
                        j,
                        k
                    );
                }
            }
        }
    }
}
