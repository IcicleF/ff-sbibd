/// Trait for *finite* fields.
pub trait Field: Send + Sync + Eq + Clone {
    /// Return the order of the field.
    fn order(&self) -> usize;

    /// Add two elements in the field.
    fn add(&self, a: usize, b: usize) -> usize;

    /// Subtract two elements in the field.
    fn sub(&self, a: usize, b: usize) -> usize;

    /// Negate an element in the field.
    fn neg(&self, a: usize) -> usize {
        self.sub(0, a)
    }

    /// Multiply two elements in the field.
    fn mul(&self, a: usize, b: usize) -> usize;

    /// Get the multiplicative inverse of an element in the field.
    fn inv(&self, a: usize) -> usize;

    /// Divide two elements in the field.
    /// This is equivalent to multiplying the first element by the inverse of the second.
    fn div(&self, a: usize, b: usize) -> usize {
        self.mul(a, self.inv(b))
    }

    /// Get a primitive element of the field.
    fn primitive(&self) -> usize {
        if self.order() == 2 {
            return 1;
        }
        'outer: for i in 2..self.order() {
            let mut x = 1;
            for _ in 1..(self.order() - 1) {
                x = self.mul(x, i);
                if x == 1 {
                    continue 'outer;
                }
            }
            return i;
        }
        unreachable!("impossible to find no primitive elements in a prime field");
    }
}
