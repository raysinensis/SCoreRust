use ndarray::prelude::*;
use rand_distr::StandardNormal;
use rand::prelude::*;
use std::collections::{HashMap, HashSet};
use extendr_api::prelude::*;
//use std::time::Instant;

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
    let v8: Vec<&str> = features.iter().map(AsRef::as_ref).collect();
    println!("{:?}", v8);
}
/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn calc_modulescore(mat: ArrayView<f64, Ix2>, features: Vec<String>, allfeatures: Vec<String>, nbin: i32, nsample:i32) -> Vec<f64> {
    // settings
    //let start = Instant::now();
    let n = nbin as usize;
    let nsamp = nsample as usize;
    let seed: u64 = 34;
    //let origgenes = vec!["a", "b", "c", "d", "e", "f", "g"];
    let origgenes: Vec<&str> = allfeatures.iter().map(AsRef::as_ref).collect();
    //let target = vec!["a", "b", "c"];
    let target: Vec<&str> = features.iter().map(AsRef::as_ref).collect();
    //mock data
    let m = mat;
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
    //println!("{:?}",start.elapsed());
    //sample matching bins as controls
    //initiate backgrounds
    let origind : Vec<usize> = (0..origgenes.len()).collect();
    let origkeys: HashMap<_, _> = origgenes.iter().cloned().zip(origind.iter()).collect();
    let mut bcv = Vec::new();
    for i in 0..res.len() {
        let temp = res[i].iter().choose_multiple(&mut StdRng::seed_from_u64(seed as u64), nsamp);
        let mut indtemp: Vec<&usize> = Vec::new();
        for x in temp {
            indtemp.push(origkeys[x])
        };
        let mut indtemp2: Vec<usize> = Vec::new();
        for &x in indtemp {
        indtemp2.push(x)
    }
        let mut scoretemp = m.select(Axis(0), &indtemp2).sum_axis(Axis(0));
        bcv.push(vec![scoretemp]);
    };
    //retrieve, and if more background is needed, add more
    let (bcv,bcs) = sort_backgrounds(bcv);
    println!("{:?}", bcv);
    println!("{:?}", bcs);
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
    let mut scores = scoretarget;
    //println!("{:?}", start.elapsed());
    return scores.to_vec();
}

fn sort_backgrounds(bc: Vec<Vec<Array1<f64>>>) -> (Vec<Vec<Array1<f64>>>, Array1<f64>) {
    let mut bcv = bc;
    bcv[1].push(array![1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 200.0, 1.0]);
    bcv[1].push(array![1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 200.0, 3.0]);

    let temp = (&bcv[1][0..=2]);
    let t2 = &temp[0] + &temp[2];
    println!("{:?}", t2);
    let t3: Array1<f64> = temp.iter().cloned().reduce(|a, b| (a + b)).unwrap();
    println!("{:?}", t3);

    return((bcv,t3));
}

fn main() {
    let m: Array<f64, _> = array![
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 200.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 400.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, -100.0, 1.0]
    ];
    let origgenes = vec!["a", "b", "c", "d", "e", "f", "g"];
    let target = vec!["a", "b", "c"];
    //let res = calc_modulescore(m.view(), 3, 2);
    //println!("{:?}", res);
    //let res2 = pass_features(vec!(String::from("ZFP36")));
    //let res3 = pass_mat(array![
    //    [1.0, 2.0, 4.0, 1.0],
    //    [-1.0, 0.0, 4.0, 1.0],
    //    [-1.0, 0.0, 4.0, 1.0]
    //].view());
    //let res4: () = sort_backgrounds();
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