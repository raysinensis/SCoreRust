extern crate nalgebra as na;
use rand_distr::StandardNormal;
use rand::prelude::*;
use std::collections::HashMap;

type MatrixBig = na::SMatrix<f64, 2, 7>;
type VectorBig = na::SVector<f64, 7>;
fn main() {
    let m = MatrixBig::from_vec(vec![0.0, 20.0, 3.0, 4.0, 5.0, 6.0, 7.0,
        1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0]
    );
    let m2 = m.row_mean_tr();
    let ran: Vec<f64> = StandardNormal.sample_iter(&mut StdRng::seed_from_u64(42)).take(m2.len()).collect();
    let mran = na::SVector::from_vec(ran.clone());
    let m3 = m2 + mran;
    let mut i: [f64; 7] = m3.into();
    let mut genes = vec!["a", "b", "c", "d", "e", "f", "g"];
    let keys: HashMap<_, _> = genes.iter().cloned().zip(i.iter()).collect();
    genes.sort_by(|a, b| keys[a].partial_cmp(keys[b]).unwrap());
    
    let n = 4;
    let groups = genes.len() / n;
    let remains = genes.len() % n;
    let split = (groups + 1)*remains;
    let mut iter = genes[..split].chunks(groups + 1).chain(genes[split..].chunks(groups));
    let res: Vec<&[&str]> = iter.collect();
    println!("{:?}", res);
    println!("{:?}", m3)
}

//fn cutn(x: Vec<String>, n: usize) -> Vec<&[&str]> {
//    let groups = x.len() / n;
//    let remains = x.len() % n;//

//    let split = (groups + 1)*remains;
//    let mut iter = x[..split].chunks(groups + 1).chain(x[split..].chunks(groups));
//    let res: Vec<&[&str]> = iter.collect();
//    println!("{:?}", res);
//    return res;
//}

    //let m2: i32 = 5;
    //let m2: Vec<f64> = vec![1.0, 20.0, 3.0];
    //let dst: Vec<&[&str]> = genes.chunks(2).collect();
    //m4.sort_by(|a, b| a.partial_cmp(b).unwrap());
    //let cuts = cutn(genes.clone(), 4);