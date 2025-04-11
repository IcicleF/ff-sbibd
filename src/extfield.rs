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
    fn new_prime(q: usize) -> Self {
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

    /// Take a power q = p^k and make it a finite field order.
    fn new(q: usize, p: usize) -> Self {
        let mut k = 0;
        let mut q2 = q;
        while q2 > 1 {
            if q2 % p == 0 {
                k += 1;
                q2 /= p;
            }
        }
        if q2 != 1 {
            panic!("{} is not a power of {}", q, p);
        }
        Order { q, p, k }
    }
}

/// A finite field that is an extension of another field.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExtField<F: Field> {
    /// The order of the field, q = p^k.
    order: Order,

    /// The prime field of order p.
    sub: F,

    /// The irreducible polynomial of degree k that generates the field.
    irr: Poly<F>,
}

impl ExtField<PrimeField> {
    /// Construct a extended field of order q = p^k, where p is a prime number.
    pub fn over_prime(q: usize) -> ExtField<PrimeField> {
        let order = Order::new_prime(q);
        ExtField::new(PrimeField::new(order.p), q)
    }
}

impl<F: Field> ExtField<F> {
    /// Construct a finite field of order q = p^k with the given irreducible polynomial,
    /// where p is the order of the subfield.
    ///
    /// # Safety
    ///
    /// `irr` must be an irreducible polynomial of degree k over the prime field of order p.
    pub unsafe fn with_irr(sub: F, q: usize, irr: Poly<F>) -> Self {
        let order = Order::new(q, sub.order());
        assert!(
            sub == *irr.field(),
            "irreducible polynomial must be over the same subfield"
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
        ExtField { order, sub, irr }
    }

    /// Construct a extended field of order q = p^k, where p is the order of the subfield.
    pub fn new(sub: F, q: usize) -> Self {
        let order = Order::new(q, sub.order());
        let irr = find_irr(sub.clone(), order.k);
        ExtField { order, sub, irr }
    }

    /// Get the irreducible polynomial of the field.
    pub fn irr(&self) -> &Poly<F> {
        &self.irr
    }

    /// Determine if the field is prime.
    pub fn is_prime(&self) -> bool {
        self.order.k == 1
    }

    /// Regard an element as a vector over the subfield.
    pub fn elem2poly(&self, a: usize) -> Poly<F> {
        assert!(a < self.order.q, "element out of range");

        let mut v = vec![0; self.order.k];
        let mut a = a;
        for i in 0..self.order.k {
            v[i] = a % self.order.p;
            a /= self.order.p;
        }
        Poly::with_coeffs(self.sub.clone(), v)
    }

    /// Regard a vector over the subfield as an element.
    pub fn poly2elem(&self, a: Poly<F>) -> usize {
        assert!(a.field() == &self.sub, "field mismatch");

        let mut n = 0;
        for i in (0..self.order.k).rev() {
            n = n * self.sub.order() + a[i];
        }
        n
    }
}

impl<F: Field> Field for ExtField<F> {
    fn order(&self) -> usize {
        self.order.q
    }

    fn add(&self, a: usize, b: usize) -> usize {
        assert!(a < self.order.q && b < self.order.q, "element out of range");
        let va = self.elem2poly(a);
        let vb = self.elem2poly(b);
        self.poly2elem(va + vb)
    }

    fn sub(&self, a: usize, b: usize) -> usize {
        assert!(a < self.order.q && b < self.order.q, "element out of range");
        let va = self.elem2poly(a);
        let vb = self.elem2poly(b);
        self.poly2elem(va - vb)
    }

    fn mul(&self, a: usize, b: usize) -> usize {
        assert!(a < self.order.q && b < self.order.q, "element out of range");
        let va = self.elem2poly(a);
        let vb = self.elem2poly(b);
        self.poly2elem((va * vb) % self.irr.clone())
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transparent() {
        let ff = ExtField::over_prime(5);
        let pf = PrimeField::new(5);

        for i in 0..5 {
            if i != 0 {
                assert_eq!(ff.inv(i), pf.inv(i));
            }

            for j in 0..5 {
                assert_eq!(ff.add(i, j), pf.add(i, j));
                assert_eq!(ff.sub(i, j), pf.sub(i, j));
                assert_eq!(ff.mul(i, j), pf.mul(i, j));
                if j != 0 {
                    assert_eq!(ff.div(i, j), pf.div(i, j));
                }
            }
        }
    }

    #[test]
    fn custom_irr_4() {
        let pf = PrimeField::new(2);

        // Irreducible polynomial: x^2 + x + 1
        let irr = Poly::with_coeffs(pf, [1, 1, 1]);
        let ff = unsafe { ExtField::with_irr(pf, 4, irr) };

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
        let ff = ExtField::over_prime(9);
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
