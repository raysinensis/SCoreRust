use ndarray::prelude::*;
use rand_distr::StandardNormal;
use rand::prelude::*;
use std::collections::{HashMap, HashSet};
use extendr_api::prelude::*;
/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn pass_mat(mat: ArrayView<f64, Ix2>) {
    println!("{:?}", mat);
}
/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn pass_features(features: Vec<String>) {
    println!("{:?}", features);
}
/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn calc_modulescore() -> Vec<f64> {
    // settings
    let n = 3;
    let nsamp = 2;
    let seed: u64 = 34;
    let origgenes = vec!["a", "b", "c", "d", "e", "f", "g"];
    let target = vec!["a", "b", "c"];
    //mock data
    let m: Array<f64, _> = array![
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 200.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 400.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, -100.0, 1.0]
    ];
    //get mean, add random noise, order genes
    let m2 = m.mean_axis(Axis(1)).unwrap();
    let ran: Vec<f64> = StandardNormal.sample_iter(&mut StdRng::seed_from_u64(seed)).take(m2.len()).collect();
    let mran = Array::from_vec(ran) / 1.0e+30;
    let m3 = m2 + mran;
    let mut genes = origgenes.clone();
    let keys: HashMap<_, _> = genes.iter().cloned().zip(m3.iter()).collect();
    genes.sort_by(|a, b| keys[a].partial_cmp(keys[b]).unwrap());
    //put into bins
    let groups = genes.len() / n;
    let remains = genes.len() % n;
    let split = (groups + 1)*remains;
    let mut iter = genes[..split].chunks(groups + 1).chain(genes[split..].chunks(groups));
    let res: Vec<&[&str]> = iter.collect();
    //sample matching bins as controls
    let mut controls: Vec<&str> = Vec::new();
    for i in 0..res.len() {
        let mut a: HashSet<_> = res[i].iter().cloned().collect();
        let mut b: HashSet<_> = target.iter().cloned().collect();
        let mut intersection = a.intersection(&b);
        let mut intv: Vec<&str> = intersection.cloned().collect();
        for j in 1..=intv.len() {
            let sample = res[i].iter().choose_multiple(&mut StdRng::seed_from_u64(seed + j as u64), nsamp);
            controls.extend(sample);
        }
    }
    let mut uniqcontrols = HashSet::new();
    controls.retain(|e| uniqcontrols.insert(*e));
    //calculate controls
    let origind : Vec<usize> = (0..origgenes.len()).collect();
    let origkeys: HashMap<_, _> = origgenes.iter().cloned().zip(origind.iter()).collect();
    let mut indcontrols: Vec<&usize> = Vec::new();
    for x in uniqcontrols {
        indcontrols.push(origkeys[x])
    }
    let mut indcontrols2: Vec<usize> = Vec::new();
    for &x in indcontrols {
        indcontrols2.push(x)
    }
    let mut scorecontrols = m.select(Axis(0), &indcontrols2).mean_axis(Axis(0)).unwrap();
    //calculate tagets
    let mut indtarget: Vec<&usize> = Vec::new();
    for x in target {
        indtarget.push(origkeys[x])
    }
    let mut indtarget2: Vec<usize> = Vec::new();
    for &x in indtarget {
        indtarget2.push(x)
    }
    let mut scoretarget = m.select(Axis(0), &indtarget2).mean_axis(Axis(0)).unwrap();
    let mut scores = scoretarget - scorecontrols;
    return scores.to_vec();
}

fn main() {
    let res = calc_modulescore();
    println!("{:?}", res);
    let res2 = pass_features(vec!(String::from("ZFP36")));
    let res3 = pass_mat(array![
        [1.0, 2.0, 4.0, 1.0],
        [-1.0, 0.0, 4.0, 1.0],
        [-1.0, 0.0, 4.0, 1.0]
    ].view());
}

extendr_module! {
    mod SCorerustR;
    fn calc_modulescore;
    fn pass_features;
    fn pass_mat;
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
//fn fred(a: Vec<f64>) -> Vec<f64> {
//println!("{:?}", a);
//a
//fred(vec![4.0, 6.0])
//}

//fn pass_mat(mat: ArrayView2<f64>) {
//    println!("{:?}", mat);
//}
//let res3 = pass_mat(array![
//    [1.0, 2.0],
//    [-1.0, 0.0]
//].view());