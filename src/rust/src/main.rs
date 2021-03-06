use ndarray::prelude::*;
use rand_distr::StandardNormal;
use rand::prelude::*;
use std::collections::{HashMap, HashSet};
use extendr_api::prelude::*;
use std::time::Instant;
use ndarray::parallel::prelude::*;
use hdf5::types::*;
use hdf5::{File, Reader};
use sprs::CsMatBase;
extern crate rayon;

fn append_h5ad(fo: &str, f: &str, m: ArrayView<f64, Ix2>, newvars: Vec<String>) {
    //copy over original file
    //except for var

    //to sparse
    let sp_mat: CsMatBase<_, usize, Vec<usize>, Vec<usize>, Vec<_>, usize> = CsMatBase::csr_from_dense(m.view(), -999999.0);
    let (indptr, indices, data) = sp_mat.clone().into_raw_storage();
    //get loc
    let file = hdf5::File::open_rw(f).expect("Unable to read file");
    let indicesloc = file.dataset("X/indices").expect("Unable to read file");
    let indptrloc = file.dataset("X/indptr").expect("Unable to read file");
    let dataloc = file.dataset("X/data").expect("Unable to read file");
    let shapeloc = file.group("X").expect("Unable to read file").attr("shape").expect("Unable to read file");
    //append vars
    let var_index_name = file
        .group("var").expect("Unable to read file")
        .attr("_index").expect("Unable to read file")
        .read_scalar::<VarLenUnicode>().expect("Unable to read file");
    let varloc = file.dataset(&format!("var/{}", var_index_name.as_str())).expect("Unable to read file");
    let var_vec = file
        .dataset(&format!("var/{}", var_index_name.as_str())).expect("Unable to read file")
        .read_1d::<VarLenUnicode>().expect("Unable to read file")
        .to_vec();
    let mut vars: Vec<String> = var_vec.iter().map(ToString::to_string).collect();
    vars.extend(newvars);
    let vars_str: Vec<&str> = vars.iter().map(AsRef::as_ref).collect();
    let mut vars2_: Vec<VarLenUnicode> = Vec::new();
    for x in vars_str {
        let mut temp : VarLenUnicode = x.parse().unwrap();
        vars2_.push(temp)
    }
    let vars_ = arr1(&vars2_);
    println!("{:?}",vars_);
    //re-shape
    let ipd = indptrloc.resize(indptr.len()).expect("Unable to write file");
    let ind = indicesloc.resize(indices.len()).expect("Unable to write file");
    let red = dataloc.resize(data.len()).expect("Unable to write file");
    //write
    let ipr = indptrloc.write(&indptr.to_owned()).expect("Unable to write file");
    let inr = indicesloc.write(&indices.to_owned()).expect("Unable to write file");
    let rer = dataloc.write(&data.to_owned()).expect("Unable to write file");
    let var = varloc.write(&vars_.slice(s![2..232]).to_owned()).expect("Unable to write file");
    let sha = arr1(&[m.ncols(), m.nrows()]);
    //let shr = shapeloc.write(&sha).expect("Unable to write file");
}

fn read_h5ad(f: &str) -> (ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>, Vec<String>) {
    println!("{:?}", "reading...");
    let file = hdf5::File::open(f).expect("Unable to read file");
    let indptr = file.dataset("X/indptr").expect("Unable to read file").read_1d::<usize>().expect("Unable to read file").to_vec();
    let indices = file.dataset("X/indices").expect("Unable to read file").read_1d::<usize>().expect("Unable to read file").to_vec();
    let data = file.dataset("X/data").expect("Unable to read file").read_1d::<f64>().expect("Unable to read file").to_vec();
    let var_index_name = file
        .group("var").expect("Unable to read file")
        .attr("_index").expect("Unable to read file")
        .read_scalar::<VarLenUnicode>().expect("Unable to read file");
    let obs_index_name = file
        .group("obs").expect("Unable to read file")
        .attr("_index").expect("Unable to read file")
        .read_scalar::<VarLenUnicode>().expect("Unable to read file");
    let var_vec = file
        .dataset(&format!("var/{}", var_index_name.as_str())).expect("Unable to read file")
        .read_1d::<VarLenUnicode>().expect("Unable to read file")
        .to_vec();
    let obs_vec = file
        .dataset(&format!("obs/{}", obs_index_name.as_str())).expect("Unable to read file")
        .read_1d::<VarLenUnicode>().expect("Unable to read file")
        .to_vec();
    let counts_mtx = CsMatBase::new((var_vec.len(), obs_vec.len()), indptr, indices, data);
    //println!("{:?}", counts_mtx.shape());
    //println!("{:?}", counts_mtx.indptr());
    let (indptr, indices, data) = counts_mtx.clone().into_raw_storage();
    //println!("{:?}", indptr);
    //println!("{:?}", indices);
    //println!("{:?}", data);
    let m2 = counts_mtx.to_dense();
    return (m2, var_vec.iter().map(ToString::to_string).collect());
}

/// Test passing expr matrix from R.
/// @export
#[extendr]
fn pass_mat(mat: ArrayView<f64, Ix2>) {
    println!("{:?}", mat);
}
/// Test passing gene list from R.
/// @export
#[extendr]
fn pass_features(features: Vec<String>) {
    let v8: Vec<&str> = features.iter().map(AsRef::as_ref).collect();
    println!("{:?}", v8);
}

/// Precalculate background nrep number of times.
fn calc_background(mat: ArrayView<f64, Ix2>, allfeatures: Vec<String>, nbin: i32, nsample:i32, nrep: i32, nthread: i32, order: Vec<String>) -> Array<f64, Ix2> {
    std::env::set_var("RAYON_NUM_THREADS", nthread.to_string());
    let n = nbin as usize;
    let nsamp = nsample as usize;
    let nreps: usize = nrep as usize;
    let seed: u64 = 34;
    //get mean, add random noise, order genes
    let g: Vec<String> = order;
    let genes: Vec<&str> = g.iter().map(AsRef::as_ref).collect();
    let origgenes: Vec<&str> = allfeatures.iter().map(AsRef::as_ref).collect();
    //put into bins
    let groups = genes.len() / n;
    let remains = genes.len() % n;
    let split = (groups + 1)*remains;
    let iter = genes[..split].chunks(groups + 1).chain(genes[split..].chunks(groups));
    let res: Vec<&[&str]> = iter.collect();
    //println!("{:?}",start.elapsed());
    //sample matching bins as controls
    let mut controls = Array2::zeros((nreps*res.len(), mat.shape()[1]));
    println!("{:?}", controls.shape());
    for i in 0..res.len() {
        for j in 1..=nreps {
            let sample = res[i].iter().choose_multiple(&mut StdRng::seed_from_u64(seed + j as u64), nsamp);
            let origind : Vec<usize> = (0..origgenes.len()).collect();
            let origkeys: HashMap<_, _> = origgenes.iter().cloned().zip(origind.iter()).collect();
            let mut indcontrols: Vec<&usize> = Vec::new();
            for x in sample {
                indcontrols.push(origkeys[x])
            }
            let mut indcontrols2: Vec<usize> = Vec::new();
            for &x in indcontrols {
                indcontrols2.push(x)
            }
            let mut scorecontrols = Vec::new();
            mat.select(Axis(0), &indcontrols2).axis_iter(Axis(1)).into_par_iter().map(|row| row.mean().unwrap()).collect_into_vec(&mut scorecontrols);   
            controls.slice_mut(s![i * nreps + j - 1,..]).assign(&Array::from_vec(scorecontrols));
        }
    }
    return controls;
}

/// order genes by expression, random breaking ties
/// @export
#[extendr]
fn order_expr(mat: ArrayView<f64, Ix2>, allfeatures: Vec<String>, nthread: i32) -> Vec<String> {
    std::env::set_var("RAYON_NUM_THREADS", nthread.to_string());
    let seed: u64 = 34;
    let origgenes = allfeatures;
    //get mean, add random noise, order genes
    let m2 = mat.mean_axis(Axis(1)).unwrap();
    let ran: Vec<f64> = StandardNormal.sample_iter(&mut StdRng::seed_from_u64(seed)).take(m2.len()).collect();
    let mran = Array::from_vec(ran) / 1.0e+30;
    let m3 = m2 + mran;
    let mut genes = origgenes.clone();
    let keys: HashMap<_, _> = genes.iter().cloned().zip(m3.iter()).collect();
    genes.sort_by(|a, b| keys[a].partial_cmp(keys[b]).unwrap());
    return genes;
}
/// Calculate pathway scoring, similar to Seurat::AddModuleScore
/// @export
#[extendr]
fn calc_modulescore(mat: ArrayView<f64, Ix2>, features: Vec<String>, allfeatures: Vec<String>, nbin: i32, nsample:i32, nthread: i32) -> Vec<f64> {
    // settings
    //let start = Instant::now();
    std::env::set_var("RAYON_NUM_THREADS", nthread.to_string());
    let n = nbin as usize;
    let nsamp = nsample as usize;
    let seed: u64 = 34;
    //get mean, add random noise, order genes
    let g: Vec<String> = order_expr(mat, allfeatures.clone(), nthread);
    let genes: Vec<&str> = g.iter().map(AsRef::as_ref).collect();
    let origgenes: Vec<&str> = allfeatures.iter().map(AsRef::as_ref).collect();
    let target: Vec<&str> = features.iter().map(AsRef::as_ref).collect();
    //put into bins
    let groups = genes.len() / n;
    let remains = genes.len() % n;
    let split = (groups + 1)*remains;
    let iter = genes[..split].chunks(groups + 1).chain(genes[split..].chunks(groups));
    let res: Vec<&[&str]> = iter.collect();
    //println!("{:?}",start.elapsed());
    //sample matching bins as controls
    let mut controls: Vec<&str> = Vec::new();
    for i in 0..res.len() {
        let a: HashSet<_> = res[i].iter().cloned().collect();
        let b: HashSet<_> = target.iter().cloned().collect();
        let intersection = a.intersection(&b);
        let intv: Vec<&str> = intersection.cloned().collect();
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
    //println!("{:?}", start.elapsed());
    let mut scorecontrols = Vec::new();
    mat.select(Axis(0), &indcontrols2).axis_iter(Axis(1)).into_par_iter().map(|row| row.mean().unwrap()).collect_into_vec(&mut scorecontrols);
    //calculate tagets
    let mut indtarget: Vec<&usize> = Vec::new();
    for x in target {
        indtarget.push(origkeys[x])
    }
    let mut indtarget2: Vec<usize> = Vec::new();
    for &x in indtarget {
        indtarget2.push(x)
    }
    let mut scoretarget = Vec::new();
    mat.select(Axis(0), &indtarget2).axis_iter(Axis(1)).into_par_iter().map(|row| row.mean().unwrap()).collect_into_vec(&mut scoretarget);
    let scores = Array::from_vec(scoretarget) - Array::from_vec(scorecontrols);
    //println!("{:?}", start.elapsed());
    //println!("{:?}", scorecontrols);
    return scores.to_vec();
}

/// Calculate pathway scoring, similar to Seurat::AddModuleScore, with precalculated order
/// @export
#[extendr]
fn calc_modulescore_orderin(mat: ArrayView<f64, Ix2>, features: Vec<String>, allfeatures: Vec<String>, nbin: i32, nsample:i32, nthread: i32, order: Vec<String>) -> Vec<f64> {
    // settings
    //let start = Instant::now();
    std::env::set_var("RAYON_NUM_THREADS", nthread.to_string());
    let n = nbin as usize;
    let nsamp = nsample as usize;
    let seed: u64 = 34;
    //get mean, add random noise, order genes
    let g: Vec<String> = order;
    let genes: Vec<&str> = g.iter().map(AsRef::as_ref).collect();
    let origgenes: Vec<&str> = allfeatures.iter().map(AsRef::as_ref).collect();
    let target: Vec<&str> = features.iter().map(AsRef::as_ref).collect();
    //put into bins
    let groups = genes.len() / n;
    let remains = genes.len() % n;
    let split = (groups + 1)*remains;
    let iter = genes[..split].chunks(groups + 1).chain(genes[split..].chunks(groups));
    let res: Vec<&[&str]> = iter.collect();
    //println!("{:?}",start.elapsed());
    //sample matching bins as controls
    let mut controls: Vec<&str> = Vec::new();
    for i in 0..res.len() {
        let a: HashSet<_> = res[i].iter().cloned().collect();
        let b: HashSet<_> = target.iter().cloned().collect();
        let intersection = a.intersection(&b);
        let intv: Vec<&str> = intersection.cloned().collect();
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
    //println!("{:?}", start.elapsed());
    let mut scorecontrols = Vec::new();
    mat.select(Axis(0), &indcontrols2).axis_iter(Axis(1)).into_par_iter().map(|row| row.mean().unwrap()).collect_into_vec(&mut scorecontrols);
    //calculate tagets
    let mut indtarget: Vec<&usize> = Vec::new();
    for x in target {
        indtarget.push(origkeys[x])
    }
    let mut indtarget2: Vec<usize> = Vec::new();
    for &x in indtarget {
        indtarget2.push(x)
    }
    let mut scoretarget = Vec::new();
    mat.select(Axis(0), &indtarget2).axis_iter(Axis(0)).into_par_iter().map(|row| row.mean().unwrap()).collect_into_vec(&mut scoretarget);
    let scores = Array::from_vec(scoretarget) - Array::from_vec(scorecontrols);
    //println!("{:?}", start.elapsed());
    //println!("{:?}", scorecontrols);
    return scores.to_vec();
}

fn main() {
    //let nt = rayon::current_num_threads();
    std::env::set_var("RUST_BACKTRACE", "1");
    //rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    let a = Array::linspace(0., 6400., 6400).into_shape((40, 160)).unwrap();
    let start = Instant::now();
    let mut sums: Vec<f64> = Vec::new();
    a.axis_iter(Axis(0)).into_par_iter().map(|row| row.sum()).collect_into_vec(&mut sums);
    println!("{:?}",start.elapsed());
    let start2 = Instant::now();
    let scorecontrols = a.mean_axis(Axis(0)).unwrap();
    println!("{:?}",start2.elapsed());
    //println!("{:?}", nt);
    let m: Array<f64, _> = array![
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 200.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 400.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, -100.0, 1.0]
    ];
    let origgenes: Vec<String> = vec!["a", "b", "c", "d", "e", "f", "g"].iter().map(|&s|s.into()).collect();
    let target: Vec<String> = vec!["a", "b", "c"].iter().map(|&s|s.into()).collect();
    let res = calc_modulescore(m.clone().view(), target.clone(), origgenes.clone(), 2, 2, 2);
    println!("{:?}", res);
    let m: Array<f64, _> = array![
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 200.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 400.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, 100.0, 1.0],
        [1.0, 2.0, -10.0, 1.0, 1.0, 2.0, -100.0, 1.0]
    ];
    let res2 = order_expr(m.view(), origgenes.clone(), 2);
    //let res3 = calc_modulescore_orderin(m.clone().view(), target.clone(), origgenes.clone(), 2, 2, 2, res2);
    let res4 = calc_background(m.clone().view(), origgenes.clone(), 2, 2, 2, 2, res2);
    println!("{:?}", res4);
    let (res4, g4) = read_h5ad("/Users/rfu/SCore-rust-dev/inst/test.h5ad");
    let res5 = res4.view();
    println!("{:?}", res5.shape());
    let res6 = ndarray::concatenate(Axis(0),&[res5,res5]).unwrap();
    println!("{:?}", res6.shape());
    let res7 = append_h5ad("/Users/rfu/SCore-rust-dev/inst/test.h5ad", "/Users/rfu/SCore-rust-dev/inst/test2.h5ad", res6.view(), g4);
    
   

    //println!("{:?}", data);
    
    //let res5 = order_expr(res4.view(), g4.clone(), 2); 
    //println!("{:?}", res5[1]);
    //let res2 = pass_features(vec!(String::from("ZFP36")));
    //let res3 = pass_mat(array![
    //    [1.0, 2.0, 4.0, 1.0],
    //    [-1.0, 0.0, 4.0, 1.0],
    //    [-1.0, 0.0, 4.0, 1.0]
    //].view());
}

extendr_module! {
    mod SCoreRust;
    fn calc_modulescore;
    fn pass_features;
    fn pass_mat;
    fn order_expr;
    fn calc_modulescore_orderin;
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