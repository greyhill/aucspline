extern crate num;

use ::num::*;
use std::ops::*;

#[derive(Clone, Debug)]
pub struct Polynomial<F>
where F: Float {
    coeffs: Vec<F>,
    left: F,
    right: F,
    shift: F,
    after: F,
}

#[derive(Clone, Debug)]
pub struct PolynomialSum<F>
where F: Float {
    poly: Vec<Polynomial<F>>,
}

impl<F> Polynomial<F>
where F: Float {
    pub fn new(left: F, right: F, coeffs: Vec<F>) -> Polynomial<F> {
        Polynomial{
            left: left,
            right: right,
            coeffs: coeffs,
            shift: F::zero(),
            after: F::zero(),
        }
    }

    fn eval_inner(self: &Self, v: F) -> F {
        let mut tr = F::zero();
        for (n, a) in self.coeffs.iter().enumerate() {
            tr = tr + *a * (v + self.shift).powi(n as i32);
        }
        tr
    }

    pub fn support(self: &Self, v: F) -> bool {
        if v + self.shift < self.left {
            false
        } else if v + self.shift >= self.right {
            false
        } else {
            true
        }
    }

    pub fn eval(self: &Self, v: F) -> F {
        if v + self.shift < self.left {
            F::zero()
        } else if v + self.shift >= self.right {
            self.after
        } else {
            self.eval_inner(v)
        }
    }

    pub fn integrate_indefinite(self: &Self) -> Polynomial<F> {
        let mut new_coeffs = Vec::<F>::with_capacity(self.coeffs.len() + 1);
        new_coeffs.push(F::zero());
        for (n, &c) in self.coeffs.iter().enumerate() {
            new_coeffs.push(c / F::from(n + 1).unwrap());
        }
        let mut tr = Polynomial{
            coeffs: new_coeffs,
            left: self.left,
            right: self.right,
            shift: self.shift,
            after: F::zero()
        };
        let f0 = tr.eval_inner(tr.left - tr.shift);
        let f1 = tr.eval_inner(tr.right - tr.shift);
        tr.coeffs[0] = tr.coeffs[0] - f0;
        tr.after = f1 - f0;

        tr
    }

    pub fn shift(self: &Self, shift: F) -> Polynomial<F> {
        Polynomial{
            coeffs: self.coeffs.clone(),
            left: self.left,
            right: self.right,
            shift: self.shift + shift,
            after: self.after,
        }
    }

    pub fn convolve_box(self: &Self, width: F, height: F) -> PolynomialSum<F> {
        let int = self.integrate_indefinite();
        let r = int.shift(width/F::from(2).unwrap());
        let l = int.shift(-width/F::from(2).unwrap());
        (r - l)*height
    }
}

impl<F> PolynomialSum<F>
where F: Float {
    pub fn eval(self: &Self, v: F) -> F {
        self.poly.iter().fold(F::zero(), |l, r| l + r.eval(v))
    }

    pub fn integrate_indefinite(self: &Self) -> PolynomialSum<F> {
        PolynomialSum{
            poly: self.poly.iter().map(|p| p.integrate_indefinite()).collect(),
        }
    }

    pub fn convolve_box(self: &Self, width: F, height: F) -> PolynomialSum<F> {
        let mut poly = Vec::new();
        for p in self.poly.iter() {
            poly.extend(p.convolve_box(width, height).poly);
        }
        PolynomialSum{
            poly: poly,
        }
    }
}

impl<'a, 'b, F> Add<&'b Polynomial<F>> for &'a Polynomial<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn add(self, rhs: &'b Polynomial<F>) -> PolynomialSum<F> {
        PolynomialSum{
            poly: vec![self.clone(), rhs.clone()],
        }
    }
}

impl<F> Add<Polynomial<F>> for Polynomial<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn add(self, rhs: Polynomial<F>) -> PolynomialSum<F> {
        PolynomialSum{
            poly: vec![self.clone(), rhs.clone()],
        }
    }
}

impl<'a, F> Neg for &'a Polynomial<F>
where F: Float {
    type Output = Polynomial<F>;

    fn neg(self) -> Polynomial<F> {
        Polynomial{
            left: self.left,
            right: self.right,
            coeffs: self.coeffs.iter().map(|&c| -c).collect(),
            shift: self.shift,
            after: -self.after,
        }
    }
}

impl<F> Neg for Polynomial<F>
where F: Float {
    type Output = Polynomial<F>;

    fn neg(self) -> Polynomial<F> {
        Polynomial{
            left: self.left,
            right: self.right,
            coeffs: self.coeffs.iter().map(|&c| -c).collect(),
            shift: self.shift,
            after: -self.after,
        }
    }
}

impl<'a, 'b, F> Sub<&'b Polynomial<F>> for &'a Polynomial<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn sub(self, rhs: &'b Polynomial<F>) -> PolynomialSum<F> {
        self + &-rhs
    }
}

impl<F> Sub<Polynomial<F>> for Polynomial<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn sub(self, rhs: Polynomial<F>) -> PolynomialSum<F> {
        self + -rhs
    }
}

impl<'a, 'b, F> Add<&'b PolynomialSum<F>> for &'a PolynomialSum<F>
where F: Float {
    type Output = PolynomialSum<F>;

    fn add(self, rhs: &'b PolynomialSum<F>) -> PolynomialSum<F> {
        let mut poly = self.poly.clone();
        poly.extend(rhs.poly.clone());
        PolynomialSum{
            poly: poly,
        }
    }
}

impl<F> Add<PolynomialSum<F>> for PolynomialSum<F>
where F: Float {
    type Output = PolynomialSum<F>;

    fn add(self, rhs: PolynomialSum<F>) -> PolynomialSum<F> {
        let mut poly = self.poly.clone();
        poly.extend(rhs.poly.clone());
        PolynomialSum{
            poly: poly,
        }
    }
}

impl<'a, F> Mul<F> for &'a Polynomial<F>
where F: Float {
    type Output = Polynomial<F>;

    fn mul(self, rhs: F) -> Polynomial<F> {
        Polynomial{
            coeffs: self.coeffs.iter().map(|&c| c*rhs).collect(),
            left: self.left,
            right: self.right,
            shift: self.shift,
            after: rhs*self.after,
        }
    }
}

impl<F> Mul<F> for Polynomial<F>
where F: Float {
    type Output = Polynomial<F>;

    fn mul(self, rhs: F) -> Polynomial<F> {
        Polynomial{
            coeffs: self.coeffs.iter().map(|&c| c*rhs).collect(),
            left: self.left,
            right: self.right,
            shift: self.shift,
            after: rhs*self.after,
        }
    }
}

impl<'a, F> Mul<F> for &'a PolynomialSum<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn mul(self, rhs: F) -> PolynomialSum<F> {
        PolynomialSum{
            poly: self.poly.iter().map(|p| p*rhs).collect()
        }
    }
}

impl<F> Mul<F> for PolynomialSum<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn mul(self, rhs: F) -> PolynomialSum<F> {
        PolynomialSum{
            poly: self.poly.iter().map(|p| p*rhs).collect()
        }
    }
}

impl<'a, F> Neg for &'a PolynomialSum<F> 
where F: Float {
    type Output = PolynomialSum<F>;

    fn neg(self) -> PolynomialSum<F> {
        PolynomialSum{
            poly: self.poly.iter().map(|p| -p).collect()
        }
    }
}

#[test]
fn test_integrate() {
    let f_1 = Polynomial::new(0f32, 1f32, vec![1f32]);
    let f_x = f_1.integrate_indefinite();
    let f_x2 = f_x.integrate_indefinite();
    let f_x3 = f_x2.integrate_indefinite();
    assert_eq!(f_x.coeffs, vec![0f32, 1f32]);
    assert_eq!(f_x2.coeffs, vec![0f32, 0f32, 0.5f32]);
    assert_eq!(f_x3.coeffs, vec![0f32, 0f32, 0f32, 1f32 / 6f32]);
}

#[test]
fn test_trap() {
    let offset = 2f32;
    let f_box = Polynomial::new(offset - 1f32, offset + 1f32, vec![1f32]);
    let f_trap = f_box.convolve_box(1f32, 1f32);
    assert_eq!(f_trap.eval(offset + -2f32), 0f32);
    assert_eq!(f_trap.eval(offset + -1.5f32), 0f32);
    assert_eq!(f_trap.eval(offset + -1f32), 0.5f32);
    assert_eq!(f_trap.eval(offset + -0.5f32), 1f32);
    assert_eq!(f_trap.eval(offset + 0f32), 1f32);
    assert_eq!(f_trap.eval(offset + 0.5f32), 1f32);
    assert_eq!(f_trap.eval(offset + 1f32), 0.5f32);
    assert_eq!(f_trap.eval(offset + 1.5f32), 0f32);
    assert_eq!(f_trap.eval(offset + 2f32), 0f32);
}

