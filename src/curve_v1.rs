extern crate num;
use num::{Integer, Signed, abs};
use std::str::FromStr;


use cosmwasm_std::{Decimal, StdResult, Uint128, Decimal256, Uint256};


// const A: i128 = 85;
// const RIGHT_LIMIT: str = "100000000000";
const LEFT_LIMIT: Decimal256 = Decimal256::zero();
const A: &str = "85";
const prec1: &str = "1000";

// Aruthmetics 
/////////////////////////////////////////////////////////////
// const DECIMAL_FRACTIONAL: Uint128 = Uint128::from(1000000_000u128);
// reverse decimal
#[warn(non_snake_case)]
pub fn r_d(decimal: Decimal256) -> Decimal256 {
    let decimal_fractional: Uint256 = Uint256::from(1000000_000u128);
    if decimal.is_zero() {
        return Decimal256::zero();
    }

    Decimal256::from_ratio(decimal_fractional, decimal * decimal_fractional)
}

// decimal subtraction
#[warn(non_snake_case)]
pub fn d_s(a: Decimal256, b: Decimal256) -> StdResult<Decimal256> {
    let decimal_fractional: Uint256 = Uint256::from(1000000_000u128);
    Ok(Decimal256::from_ratio(
        a * decimal_fractional - b * decimal_fractional,
        decimal_fractional,
    ))
}


// decimal multiplication
#[warn(non_snake_case)]
pub fn d_m(a: Decimal256, b: Decimal256) -> Decimal256 {
    let decimal_fractional: Uint256 = Uint256::from(1000000_000u128);
    Decimal256::from_ratio(a * decimal_fractional * b, decimal_fractional)
}
/////////////////////////////////////////////////////////////

pub   fn abs_value(n1: & Decimal256, n2: & Decimal256) -> Decimal256 {
    //let val: Decimal = d_s(*n1, *n2).unwrap();
    
    if PartialOrd::le(n1, n2) {
        return d_s(*n2, *n1).unwrap();
    }
    else {
        return d_s(*n1, *n2).unwrap();
    }
}

/* ----------------------------------------------------------- */
/// A trait for things that can be approximately equal.

pub trait Epsilon {
    type RHS;
    type Precision;

    /// Return true if self and `other` differ no more than by a given amount.
    fn close(&self, other: Self::RHS, precision: Self::Precision) -> bool;

    /// Return true if self is close to zero.
    fn near_zero(&self, precision: Self::Precision) -> bool;
}
// impl<T> std::marker::Copy for &T {}

/*
impl<T> Clone for T {
    fn clone(&self) -> Self {
        *self
    }
}
*/

impl Epsilon for Decimal256 {
    type RHS = Decimal256;
    type Precision = Decimal256;

    fn close(&self, other: Decimal256, precision: Decimal256) -> bool {
       // let val: Decimal = abs_value(self, & mut precision);
        Decimal256::le(&abs_value( &  other, self), &precision)
    }

    fn near_zero(&self, precision: Decimal256) -> bool {
        let t: Decimal256 = abs_value(& self, & Decimal256::zero());
        return Decimal256::lt(&t, &precision);
        // abs(*self) < abs(precision)
    }
}


/// Configuration structure for the Newton's method (one root version).
#[derive(Debug, Clone, Copy)]
pub struct OneRootNewtonCfg {
    /// The real root, if any, is most likely to be within this distance from
    /// the reported root, but this is not guaranteed.
    pub precision: Decimal256,
    /// A limit on the number of iterations to perform. Pass `None` if you
    /// don't want a limit.
    pub max_iters: Option<u32>
}

pub fn get_function_value(arr:[Decimal256;3], val:Decimal256, indx: u32) -> Decimal256 {
    if indx == 0 {
        
        // return d_m(val, arr[0]) + d_m(d_m(d_m(val, val), val) , r_d(arr[1])) - arr[2];
        let v1: Decimal256 = d_m(val, r_d(arr[1]));
        println!("v1=  {}", v1);
        let v2: Decimal256 = d_m(v1, val);
        println!("v2=  {}", v2);
        let v3: Decimal256 = d_m(v2, val);
        println!("v3=  {}", v3);
        return d_s(v3, arr[2]).unwrap();
    }
    else  {
        return arr[0] + d_m(arr[1], r_d( val)) - d_m(arr[2], val);
    }
}


pub fn get_deriv_value(arr:[Decimal256;2], val:Decimal256, indx: u32) -> Decimal256 {
    if indx == 0 {
        return arr[0] +  d_m(val, val);
    }
    else {
        return d_s(Decimal256::one(), d_m(arr[0], r_d(d_m(val, val)))).unwrap();
        //return -arr[0]/ (val * val) - arr[1];
    }
}




pub fn newton_one(config: OneRootNewtonCfg,
    first_approx: Decimal256,
    a_func: [Decimal256;3],
    a_der: [Decimal256; 2],
    indx: u32
) -> Option<Decimal256>
{
    let mut left: Decimal256 = LEFT_LIMIT;
    let mut right: Decimal256 =  Decimal256::from_str("100000000000").unwrap();
    let mut left_val: Decimal256 = get_function_value(a_func, left, indx);
    let mut right_val: Decimal256 = get_function_value(a_func, right, indx);
    let mut root = first_approx;
    let mut prev_root = None;
    let mut iter = 0;

    while prev_root.map_or(true, |old| !root.close(old, config.precision))
            && config.max_iters.map_or(true, |max| iter < max) {
        iter += 1;
        if let Some(next) = next_newton_iter(config.precision,
                          left, 
                          right, 
                          root, 
                          a_func, 
                          a_der,
                        indx) {
            prev_root = Some(root);
            root = next;
        } else if let Some(fallback_root) = linear_fallback(left, right, left_val, right_val) {
            prev_root = Some(root);
            root = fallback_root;
        } else {
            return None
        }
        let val_at_root = get_deriv_value(a_der, root, indx);
       // if left_val * val_at_root <= 0 {
        if Decimal256::le(&d_m( left_val, val_at_root), &Decimal256::zero()) {
            right = root;
            right_val = val_at_root;
        } else {
            left = root;
            left_val = val_at_root;
        }
    }
    return Some(root);
}


fn next_newton_iter(prec: Decimal256,
        left: Decimal256,
        right: Decimal256,
        old: Decimal256,
        a_func: [Decimal256;3],
        a_dir: [Decimal256;2],
        indx: u32) -> Option<Decimal256>
{
    let d = get_deriv_value(a_dir, old, indx);
    if d.near_zero(prec) {
        return None
    }
    //let res = old - get_function_value(a_func, old, indx) / d;
    let res: Decimal256 = d_s(old, d_m(get_function_value(a_func, old, indx), r_d(d))).unwrap();
    if res < left {
        None
    } else if res > right {
        None
    } else {
    Some(res)
    }
}

fn linear_fallback(x1: Decimal256 , x2: Decimal256, y1: Decimal256, y2: Decimal256) -> Option<Decimal256>
{
    // let res = ((y2 - y1) * x1 - (x2 - x1) * y1) / (y2 - y1);
    let t1: Decimal256 = d_s(y2, y1).unwrap();
    let t2: Decimal256 = d_s(x2, x1).unwrap();
    let t3: Decimal256 = d_s(d_m(t1, x1), d_m(t2, y1)).unwrap();
    let res: Decimal256 = d_m(t3, r_d(t1));
    if res < x1 {
        None
    } else if res > x2 {
        None
    } else {
        Some(res)
    }
}

pub fn compute_d(op: Decimal256, ap: Decimal256) -> Option<Decimal256>
{
    let n4: Decimal256 = Decimal256::from_str("4").unwrap();
    let n3: Decimal256 = Decimal256::from_str("3").unwrap();

    println!("Start newton b 1 n3, n4 = {}, {}", n3, n4);
    let sum: Decimal256 = op + ap;
    let prod: Decimal256 =  d_m(op, ap);

    println!("Start newton b 1 sum, prod = {}, {}", sum, prod);

    let a4: Decimal256 = d_m(Decimal256::from_str(A).unwrap(), n4);
    let prod4: Decimal256 = d_m(n4, prod);
    let a4_1: Decimal256 = d_s(a4, Decimal256::one()).unwrap();
    let a4_sum: Decimal256 = d_m(a4 , sum);
   // let prod4_3: i128 = 3.0 / prod4;
    let prod4_3: Decimal256 = d_m(n3, r_d(prod4));
   let d0_0: Decimal256 = op + ap;
    //let prec = 10;
    let _cfg = OneRootNewtonCfg {
        precision: Decimal256::from_str(prec1).unwrap(),
        max_iters: None
    };


    let arr_func:[Decimal256; 3] = [a4_1, prod4, a4_sum];
    let arr_deriv:[Decimal256; 2] = [a4_1, prod4_3];

    println!(" arr func = {}, {}, {}", arr_func[0], arr_func[1], arr_func[2]);

    println!(" arr deriv = {}, {}", arr_deriv[0], arr_deriv[1]);

    println!("Start newton");

    let sol = newton_one(_cfg, d0_0, arr_func, arr_deriv, 0);

    return sol;
}

/*
pub fn get_ask_amount(op: i128, of: i128, d: i128) -> i128
{
    let a4: i128 = 4 * A;
    let x1: i128 = op + of;
    let _target_y = |x: i128| a4 * d + d * d * d / (4 * x1 * x) -  a4 * ( x1 + x) - d;
    let _der_y = |x: i128| (- d) * d * d / (4 * x1 * x *x) - a4;

    let prec = 1000;
    let _cfg = OneRootNewtonCfg {
        precision: prec,
        max_iters: None
    };

   // let sol = newton_one(_cfg, 0.0, 10e25, of, &_target_y, &_der_y);
   let sol = newton_one(_cfg, 0, i128::pow(10, 25), d, , &_der_y);

    let ask_pool: i128;

    return sol;
}
*/

/*
pub fn curve_v1(_offer_pool: u128, _ask_pool: u128, _offer: u128)  -> u128
{
    let op: i128 = _offer_pool as i128;
    let ap: i128 = _ask_pool as i128;
    let of: i128 = _offer as i128;
    //let d0 = 1000.0;
    let d0: i128 = op + ap as i128;
    let d = compute_d(op, ap, d0);

    //println!("d version 2 = {0}", d);

    let ask_f = get_ask_amount(op, of, d);
    let ask_amnt: u128 = ask_f as u128;
    return _ask_pool - ask_amnt;
}
*/