extern crate nalgebra as na;
use rand_distr::StandardNormal;
use rand::prelude::*;
use std::collections::{HashMap, HashSet};

//setup
type MatrixBig = na::SMatrix<f64, 2, 7>;
type VectorBig = na::SVector<f64, 7>;
type DMatrixf64 = na::OMatrix<f64, na::Dynamic, na::Dynamic>;
fn main() {
    // settings
    let n = 2;
    let nsamp = 2;
    let seed: u64 = 34;

    //mock data
    let m = MatrixBig::from_vec(vec![-100.0, 20.0, 3.0, 4.0, 5.0, 6.0, 7.0,
        1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0]
    );
    //get mean, add random noise, order genes
    let m2 = m.row_mean_tr();
    let ran: Vec<f64> = StandardNormal.sample_iter(&mut StdRng::seed_from_u64(seed)).take(m2.len()).collect();
    let mran = na::SVector::from_vec(ran.clone()) / 1.0e+30;
    let m3 = m2 + mran;
    let mut i: [f64; 7] = m3.into();
    let origgenes = vec!["a", "b", "c", "d", "e", "f", "g"];
    let mut genes = origgenes.clone();
    let keys: HashMap<_, _> = genes.iter().cloned().zip(i.iter()).collect();
    genes.sort_by(|a, b| keys[a].partial_cmp(keys[b]).unwrap());
    //put into bins
    
    let groups = genes.len() / n;
    let remains = genes.len() % n;
    let split = (groups + 1)*remains;
    let mut iter = genes[..split].chunks(groups + 1).chain(genes[split..].chunks(groups));
    let res: Vec<&[&str]> = iter.collect();
    println!("{:?}", res);
    println!("{:?}", m3);
    //sample matching bins as controls
    let target = vec!["a", "b", "f"];
    let mut controls: Vec<&str> = Vec::new();
    for i in 0..res.len() {
        let mut a: HashSet<_> = res[i].iter().cloned().collect();
        let mut b: HashSet<_> = target.iter().cloned().collect();
        let mut intersection = a.intersection(&b);
        let mut intv: Vec<&str> = intersection.cloned().collect();
        for j in 1..=intv.len() {
            let sample = res[i].iter().choose_multiple(&mut StdRng::seed_from_u64(seed + j as u64), nsamp);
            println!("{:?}", sample);
            controls.extend(sample);
        }
    }
    let mut uniqcontrols = HashSet::new();
    controls.retain(|e| uniqcontrols.insert(*e));
    println!("{:?}", uniqcontrols);
    let origind : Vec<usize> = (0..origgenes.len()).collect();
    let origkeys: HashMap<_, _> = origgenes.iter().cloned().zip(origind.iter()).collect();
    let mut indcontrols: Vec<&usize> = Vec::new();
    for x in uniqcontrols {
        indcontrols.push(origkeys[x])
    }
    let mut scorecontrols = m.select_columns(indcontrols.iter().cloned()).column_mean();
    println!("{:?}", scorecontrols);
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
    //let mut scorecontrols = DMatrixf64::zeros(1, uniqcontrols.len());