//! Arithmetic table.

use crate::field::Field;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct ArithmeticTable {
    /// The field order.
    pub order: usize,

    /// Addition table.
    pub add: Vec<Vec<usize>>,

    /// Negation table.
    pub neg: Vec<usize>,

    /// Multiplication table.
    pub mul: Vec<Vec<usize>>,

    /// Inversion table.
    pub inv: Vec<usize>,
}

impl ArithmeticTable {
    pub fn new<F: Field>(f: &F) -> Arc<Self> {
        let order = f.order();

        let mut add = vec![vec![0; order]; order];
        let mut neg = vec![0; order];
        let mut mul = vec![vec![0; order]; order];
        let mut inv = vec![0; order];

        for i in 0..order {
            for j in 0..order {
                add[i][j] = f.add(i, j);
                mul[i][j] = f.mul(i, j);
            }
            neg[i] = f.neg(i);
            if i != 0 {
                inv[i] = f.inv(i);
            }
        }

        Arc::new(Self {
            order,
            add,
            neg,
            mul,
            inv,
        })
    }
}
