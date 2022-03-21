use cosmwasm_std::Uint256;
use cosmwasm_std::Uint128;
 use cosmwasm_std::Decimal;
 use cosmwasm_std::Decimal256;

pub mod curve_v1;

fn main() {
    let op: Uint256 = Uint256::from(40000000000u128);
    let ap: Uint256 = Uint256::from(25000000000u128);
    let of: Uint256 = Uint256::from(2500000000u128);

    let op1: Decimal256 = Decimal256::from_ratio(op, 1u128);
    let ap1: Decimal256 = Decimal256::from_ratio(ap, 1u128);
    let of1: Decimal256 = Decimal256::from_ratio(of, 1u128);

    let sol: Decimal256 = curve_v1::compute_d(op1, ap1).unwrap();

   //let sol1 = curve_v1::curve_v1(op, ap, of);

   //let sol = curve_v1::compute_d(op, ap);
 
    println!("Solution d = {:?}", sol);

    //let res = of / ap;

    //println!("res = {}", res);

}
