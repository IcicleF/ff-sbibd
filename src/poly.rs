use crate::field::Field;
use std::fmt::Debug;
use std::ops::{Add, Div, Index, Mul, Neg, Rem, Sub};
use std::vec;

/// A unary polynomial over a field.
#[derive(Clone, PartialEq, Eq)]
pub struct Poly<F: Field> {
    f: F,
    coeffs: Vec<usize>,
}

impl<F: Field> Poly<F> {
    fn normalize(&mut self) {
        // Remove trailing zero coefficients
        while self.coeffs.len() > 0 && self.coeffs[self.coeffs.len() - 1] == 0 {
            self.coeffs.pop();
        }
    }

    /// Create a new polynomial with zero coefficients.
    pub fn new(f: F) -> Self {
        Self { f, coeffs: vec![] }
    }

    /// Create a new polynomial with the given coefficients.
    /// The coefficients are in increasing order of degree.
    pub fn with_coeffs(f: F, coeffs: impl IntoIterator<Item = usize>) -> Self {
        let coeffs = coeffs.into_iter().collect::<Vec<_>>();
        debug_assert!(
            coeffs.iter().all(|&c| c < f.order()),
            "coefficients must be in the field",
        );

        // Remove trailing zero coefficients
        let mut p = Self { f, coeffs };
        p.normalize();
        p
    }

    /// Return the degree of the polynomial, or `None` if the polynomial is zero.
    pub fn degree(&self) -> Option<usize> {
        return if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        };
    }

    /// Return the coefficient field.
    pub fn field(&self) -> &F {
        &self.f
    }

    /// Set the coefficient at the given index.
    pub fn set(&mut self, index: usize, value: usize) {
        if index >= self.coeffs.len() {
            self.coeffs.resize(index + 1, 0);
        }
        self.coeffs[index] = value;
        self.normalize();
    }
}

impl<F: Field> Debug for Poly<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.is_empty() {
            return write!(f, "0");
        }
        let mut s = String::new();
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if *coeff != 0 {
                if !s.is_empty() {
                    s.push_str(" + ");
                }
                match (i, *coeff) {
                    (0, _) => s.push_str(&format!("{}", coeff)),
                    (1, 1) => s.push_str("x"),
                    (1, _) => s.push_str(&format!("{}x", coeff)),
                    (_, 1) => s.push_str(&format!("x^{}", i)),
                    _ => s.push_str(&format!("{}x^{}", coeff, i)),
                }
            }
        }
        write!(f, "{}", s)
    }
}

impl<F: Field> Add for Poly<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        debug_assert!(self.f == rhs.f, "polynomials must be over the same field");

        let mut coeffs = vec![0; self.coeffs.len().max(rhs.coeffs.len())];
        for i in 0..coeffs.len() {
            coeffs[i] = self.f.add(
                if i < self.coeffs.len() {
                    self.coeffs[i]
                } else {
                    0
                },
                if i < rhs.coeffs.len() {
                    rhs.coeffs[i]
                } else {
                    0
                },
            );
        }
        Poly { f: self.f, coeffs }
    }
}

impl<F: Field> Sub for Poly<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        debug_assert!(self.f == rhs.f, "polynomials must be over the same field");

        let mut coeffs = vec![0; self.coeffs.len().max(rhs.coeffs.len())];
        for i in 0..coeffs.len() {
            coeffs[i] = self.f.sub(
                if i < self.coeffs.len() {
                    self.coeffs[i]
                } else {
                    0
                },
                if i < rhs.coeffs.len() {
                    rhs.coeffs[i]
                } else {
                    0
                },
            );
        }
        Poly::with_coeffs(self.f, coeffs)
    }
}

impl<F: Field> Neg for Poly<F> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let coeffs = self
            .coeffs
            .iter()
            .map(|&c| self.f.neg(c))
            .collect::<Vec<_>>();
        Poly::with_coeffs(self.f, coeffs)
    }
}

impl<F: Field> Mul for Poly<F> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        debug_assert!(self.f == rhs.f, "polynomials must be over the same field");

        if self.degree().is_none() || rhs.degree().is_none() {
            return Poly::new(self.f);
        }

        let mut coeffs = vec![0; self.coeffs.len() + rhs.coeffs.len() - 1];
        for i in 0..self.coeffs.len() {
            for j in 0..rhs.coeffs.len() {
                coeffs[i + j] = self
                    .f
                    .add(coeffs[i + j], self.f.mul(self.coeffs[i], rhs.coeffs[j]));
            }
        }
        Poly::with_coeffs(self.f, coeffs)
    }
}

impl<F: Field> Div for Poly<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        debug_assert!(self.f == rhs.f, "polynomials must be over the same field");
        debug_assert!(rhs.degree().is_some(), "division by zero");

        if self.degree().is_none() {
            return self;
        }

        let d1 = self.degree().unwrap();
        let d2 = rhs.degree().unwrap();
        if d1 < d2 {
            return Poly::new(self.f);
        }

        let mut coeffs = self.coeffs.clone();
        let mut ret = vec![0; d1 - d2 + 1];

        for d in (d2..=d1).rev() {
            let c = self.f.div(coeffs[d], rhs.coeffs[d2]);
            ret[d - d2] = c;
            for i in 0..=d2 {
                coeffs[d - i] = self.f.sub(coeffs[d - i], self.f.mul(c, rhs.coeffs[d2 - i]));
            }
        }
        Poly::with_coeffs(self.f, ret)
    }
}

impl<F: Field> Rem for Poly<F> {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        debug_assert!(self.f == rhs.f, "polynomials must be over the same field");
        debug_assert!(rhs.degree().is_some(), "division by zero");

        if self.degree().is_none() {
            return self;
        }

        let d1 = self.degree().unwrap();
        let d2 = rhs.degree().unwrap();
        if d1 < d2 {
            return self;
        }

        let mut coeffs = self.coeffs.clone();
        for d in (d2..=d1).rev() {
            let c = self.f.div(coeffs[d], rhs.coeffs[d2]);
            for i in 0..=d2 {
                coeffs[d - i] = self.f.sub(coeffs[d - i], self.f.mul(c, rhs.coeffs[d2 - i]));
            }
        }
        Poly::with_coeffs(self.f, coeffs)
    }
}

impl<F: Field> Index<usize> for Poly<F> {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        if index < self.coeffs.len() {
            &self.coeffs[index]
        } else {
            &0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn debug() {
        use crate::prime::PrimeField;
        let pf = PrimeField::new(2);

        // x^2 + 1
        let p1 = Poly::with_coeffs(pf, [1, 0, 1]);

        // x^3 + x^2 + x + 1
        let p2 = Poly::with_coeffs(pf, [1, 1, 1, 1]);

        assert_eq!(format!("{:?}", p1), "1 + x^2");
        assert_eq!(format!("{:?}", p2), "1 + x + x^2 + x^3");
        assert_eq!(format!("{:?}", Poly::new(pf)), "0");
        assert_eq!(format!("{:?}", Poly::with_coeffs(pf, [0])), "0");
        assert_eq!(format!("{:?}", Poly::with_coeffs(pf, [1])), "1");
    }

    #[test]
    fn pf() {
        use crate::prime::PrimeField;
        let pf = PrimeField::new(2);

        // x^2 + 1
        let p1 = Poly::with_coeffs(pf, [1, 0, 1]);

        // x^3 + x^2 + x + 1
        let p2 = Poly::with_coeffs(pf, [1, 1, 1, 1]);

        // Degree
        assert_eq!(Poly::new(pf).degree(), None);
        assert_eq!(p1.degree(), Some(2));
        assert_eq!(p2.degree(), Some(3));

        // Indexing
        assert_eq!(p1[0], 1);
        assert_eq!(p1[1], 0);
        assert_eq!(p1[2], 1);
        assert_eq!(Poly::new(pf)[3], 0);

        // Arithmetic
        assert_eq!(p1.clone() + p2.clone(), Poly::with_coeffs(pf, [0, 1, 0, 1]));
        assert_eq!(p1.clone() - p2.clone(), Poly::with_coeffs(pf, [0, 1, 0, 1]));
        assert_eq!(
            p1.clone() * p2.clone(),
            Poly::with_coeffs(pf, [1, 1, 0, 0, 1, 1])
        );
        assert_eq!(p2.clone() / p1.clone(), Poly::with_coeffs(pf, [1, 1]));
        assert_eq!(p2.clone() % p1.clone(), Poly::new(pf));
    }
}
