//! Prime finite fields.

use crate::field::Field;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PrimeField {
    p: usize,
}

impl PrimeField {
    /// Construct a prime field of order p.
    pub fn new(p: usize) -> Self {
        debug_assert!(p > 1, "invalid order {}", p);

        // Ensure p is a prime.
        for i in 2..(p as f64).sqrt() as usize + 1 {
            if p % i == 0 {
                panic!("invalid order {}", p);
            }
        }
        Self { p }
    }
}

impl Field for PrimeField {
    fn order(&self) -> usize {
        self.p
    }

    fn add(&self, a: usize, b: usize) -> usize {
        (a + b) % self.p
    }

    fn sub(&self, a: usize, b: usize) -> usize {
        (a + self.p - b) % self.p
    }

    fn mul(&self, a: usize, b: usize) -> usize {
        (a * b) % self.p
    }

    fn inv(&self, a: usize) -> usize {
        debug_assert!(a != 0, "division by zero");

        fn mod_pow(mut base: usize, mut exp: usize, modulus: usize) -> usize {
            let mut result = 1;
            base = base % modulus;
            while exp > 0 {
                if exp % 2 == 1 {
                    result = result * base % modulus;
                }
                base = base * base % modulus;
                exp /= 2;
            }
            result
        }

        mod_pow(a, self.p - 2, self.p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        let field = PrimeField { p: 7 };

        assert_eq!(field.order(), 7);
        assert_eq!(field.add(3, 4), 0);
        assert_eq!(field.sub(3, 4), 6);
        assert_eq!(field.mul(3, 4), 5);
        assert_eq!(field.inv(3), 5);
        assert_eq!(field.div(3, 4), 6);
    }

    #[test]
    fn primitive() {
        assert_eq!(PrimeField::new(2).primitive(), 1);
        assert_eq!(PrimeField::new(3).primitive(), 2);
        assert_eq!(PrimeField::new(5).primitive(), 2);
        assert_eq!(PrimeField::new(7).primitive(), 3);
    }

    #[test]
    #[should_panic(expected = "division by zero")]
    fn div_by_zero() {
        let field = PrimeField { p: 7 };
        field.div(3, 0);
    }

    #[test]
    #[should_panic(expected = "invalid order 6")]
    fn invalid_order() {
        PrimeField::new(6);
    }
}
