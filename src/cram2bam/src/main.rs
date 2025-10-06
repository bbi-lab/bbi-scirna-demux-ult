#![allow(unused_parens)]

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;


/*
** Maximum number of microtiter plates of any type, e.g. RT, Ligation, PCR.
*/
const MAX_NUM_PLATES: usize = 16;
const MAX_READ_LENGTH: usize = 2048;


use std::collections::{HashMap, HashSet};

extern crate clap;
use clap::{Arg, Command};
use rust_htslib::errors::Error;
use rust_htslib::bam::{Read, Reader};
use rust_htslib::bam::record::Aux;
use rust_htslib::tpool::ThreadPool;
use csv::{ReaderBuilder, Trim};
use serde::Deserialize;
use base_custom::BaseCustom;
use regex::Regex;
use chrono;


/*
** Notes on processing:
**   o  convert barcode indices to base 4 sequences
**   o  convert PCR wells to base 4 sequences
**   o  make filenames using encoded indices
**   o  write BAM file (no need to convert seq/qual/umi)
**   o  construct read name from bits of information...
**   o  review bbi-rnaseq-demux rna_rt_lig demux program
*/

/*
** Define command line arguments.
*/
fn set_cl_options() -> Result<clap::Command, Box<dyn std::error::Error>> {
  let cl_options = Command::new("cram2bam")
        .version(env!("CARGO_PKG_VERSION"))
        .about("Convert Ultima ucram files to BBI ubam files.")
        .arg(Arg::new("ucram_file")  // required=true, no default
                  .required(true)
                  .short('i')
                  .long("ucram_file")
                  .help("Input Ultima ucram file path."))
        .arg(Arg::new("ubam_file")  // required=true, no default
                  .required(true)
                  .short('o')
                  .long("ubam_file")
                  .help("Output ubam file path."))
        .arg(Arg::new("sample_sheet")   // required=true, no default
                  .required(true)
                  .short('s')
                  .long("sample_sheet")
                  .help("Sample sheet file path."))
        .arg(Arg::new("rt_barcode_file")   // required=true, no default
                  .required(true)
                  .short('r')
                  .long("rt_barcodes")
                  .help("RT barcode file path."))
        .arg(Arg::new("ligation_barcode_file")   // required=true, no default
                  .required(true)
                  .short('l')
                  .long("ligation_barcodes")
                  .help("Ligation barcode file path."))
        .arg(Arg::new("lane_index")   // required=true, no default
                  .required(true)
                  .long("lane_index")
                  .help("Lane index."))
        .arg(Arg::new("pcr7_index")   // required=true, no default
                  .required(true)
                  .short('7')
                  .long("pcr7_index")
                  .help("PCR 7 index."))
        .arg(Arg::new("pcr5_index")   // required=true, no default
                  .required(true)
                  .short('5')
                  .long("pcr5_index")
                  .help("PCR 5 index."))
        .arg(Arg::new("number_threads")   //required=false
                  .required(false)
                  .short('t')
                  .long("ncpu")
                  .number_of_values(1)
                  .default_value("2")
                  .value_parser(clap::value_parser!(u32))
                  .help("Number of threads to use for reading and writing CRAM/BAM file."));

    Ok(cl_options)
}


pub fn index_to_well(mut well_index: usize, across_row_first: bool) -> Result<(usize, String), Error> {

  /*
  ** well_index == 0 is none.
  */
  if(well_index < 1) {
    return( Ok((0_usize, String::from("none"))) );
  } 

  well_index = well_index - 1;
  let nrow: usize = 8;
  let ncol: usize = 12;
  let ipl: usize = well_index / 96;
  let i96: usize = well_index - ipl * 96;
  let well_row: char;
  let well_col: usize;
  if(across_row_first) {
    well_row = char::from_u32((65 + i96 / ncol).try_into().unwrap()).expect("unable to convert usize to u8");
    well_col = (i96 % ncol) + 1;
  }
  else {
    well_row = char::from_u32((65 + (i96 % nrow)).try_into().unwrap()).expect("unable to convert usize to u8");
    well_col = i96 / nrow + 1;
  }

  let well_id: String = format!("{}{:02}", well_row, well_col);

  /*
  ** Return tuple of ipl and well_id
  */
  Ok((ipl + 1, well_id))
}


/*
** Make a HashMap that maps utiter plate wells to
** well indices. Well A01 has index 1.
*/
fn make_well_index_map(max_index: usize, across_row_first: bool, with_plate: bool) -> Result<HashMap<String, usize>, Error> {
  let mut well_index_map: HashMap<String, usize> = HashMap::with_capacity(max_index);
  let mut ipl: usize;
  let mut well: String;
  for i in 1..(max_index+1) {
    (ipl, well) = index_to_well(i, across_row_first).unwrap();
    if(with_plate) {
      let plate_well: String = format!("P{:02}-{}", ipl, well);
      let _ = well_index_map.insert(plate_well, i);
    }
    else {
      let _ = well_index_map.insert(well, i);
    }
  }
  Ok(well_index_map)
}


fn make_rt_well_to_index_map(max_num_plates: usize) -> Result<HashMap<String, usize>, Error> {
  let plate_well_index_map_by_row = make_well_index_map(max_num_plates * 96 + 1, true, true).unwrap();
  Ok(plate_well_index_map_by_row)
}


fn make_rt_index_to_well_map(max_num_plates: usize) -> Result<Vec<String>, Error> {
  let mut rt_index_to_well_map: Vec<String> = Vec::with_capacity(max_num_plates * 96 + 1);
  rt_index_to_well_map.push("Undetermined".to_string());
  for well_index in (1..(max_num_plates * 96 + 1)) {
    let (ipl, well_string) = index_to_well(well_index, true).unwrap();
    let well_name = format!("P{:02}_{}", ipl, well_string);
    rt_index_to_well_map.push(well_name);
  }

  Ok(rt_index_to_well_map)
}


fn make_lig_well_to_index_map(max_num_plates: usize) -> Result<HashMap<String, usize>, Error> {
  let mut lig_well_to_index_map: HashMap<String, usize> = HashMap::with_capacity(max_num_plates * 96 + 1);

  for i in 0..(max_num_plates * 96 + 1) {
    let well_name = format!("LIG{}", i);
    lig_well_to_index_map.insert(well_name, i);
  }

  Ok(lig_well_to_index_map)
}


/*
** Make a vector of 'sample range' structs for the samples
** in the specified lane. Separate the barcode ranges and
** expand them into sorted lists of individual indices, not
** ranges.
**
** sample ranges hashmap
**   {
**       "sample_id": "sample.2",
**       "ranges": "32-35,100-104:1-12:17-24",
**       "lanes": "5-8",
**       "tissue": "lung",
**       "genome": "mouse",
**       "hash_file": "",
**       "sample_flags": "S",
**       "external_sample_name": "extsample.2",
**       "wrap_group": "Smyth",
**       "rt_file": "",
**       "ligation_file": "",
**       "p7_file": "/net/bbi/barcodes2/p7_file.txt",
**       "p5_file": "/net/bbi/barcodes2/p5_file.txt",
**       "library": "RNA3-088"
**      }
*/

#[allow(dead_code)]
#[derive(Deserialize, Debug, Clone)]
struct SampleMap {
  sample_id: String,
  ranges: String,
  lanes: String,
  tissue: String,
  genome: String,
  hash_file: String,
  sample_flags: String,
  external_sample_name: String,
  wrap_group: String,
  rt_file: String,
  ligation_file: String,
  p7_file: String,
  p5_file: String,
  library: String,
  process_group: String
}


/*
** Make a sorted vector of distinct indices as usize values from
** a String slice of index 'ranges' from the samplesheet.json file.
** For example, &str = "1-4,5,7".
*/
fn make_index_vec(index_str: &str) -> Result<Vec<usize>, Error> {
  let re = Regex::new(r"([0-9]+)([-]([0-9]+))?").unwrap();
  let mut index_vec: Vec<usize> = Vec::new();
  for index_range in index_str.split(",") { 
    let re_mats = re.captures(index_range).unwrap();
    let index1: usize = re_mats.get(1).unwrap().as_str().parse().unwrap();
    let mut index2: usize = index1;
    if(re_mats.get(3) != None) {
      index2 = re_mats.get(3).unwrap().as_str().parse().unwrap();
    }
    for idx in index1..(index2+1) {
      index_vec.push(idx);
    }
  }

  /*
  ** Find distinct indices.
  */
  let mut index_set: HashSet<usize> = HashSet::new();
  for index in index_vec.iter() {
    index_set.insert(*index);
  } 

  let mut index_vec_distinct: Vec<usize> = Vec::new();
  for index in index_set {
    index_vec_distinct.push(index);
  }

  /*
  ** Sort indices.
  */
  index_vec_distinct.sort();

  Ok(index_vec_distinct)
}


fn get_sample_index_vecs(sample_map: &SampleMap) -> Result<(Vec<usize>, Vec<usize>, Vec<usize>, Vec<usize>), Error> {
  let lane_index_vec: Vec<usize> = make_index_vec(&sample_map.lanes).unwrap();

  let ranges_parts: Vec<&str> = sample_map.ranges.split(":").collect();
  let rt_index_vec: Vec<usize> = make_index_vec(ranges_parts[0]).unwrap();
  let p7_index_vec: Vec<usize> = make_index_vec(ranges_parts[1]).unwrap();
  let p5_index_vec: Vec<usize> = make_index_vec(ranges_parts[2]).unwrap();

  Ok((lane_index_vec, rt_index_vec, p7_index_vec, p5_index_vec))
}


/*   
** Deserialize a SampleMap.
**
** See
**   https://docs.rs/serde_json/latest/serde_json/fn.from_value.html
** and
**   https://docs.rs/serde_json/latest/serde_json/value/enum.Value.html
**
*/
fn deserialize_sample_map(serialized_sample_map: serde_json::Value) -> Result<SampleMap, Error> {

  let sample_map: SampleMap = serde_json::from_value(serialized_sample_map).expect("Error: deserialize_sample_map: unable to deserialize the sample map.\n  Perhaps the samplesheet sample_index_list map changed.\n  If so, update the SampleMapStrings struct in this program.");

  Ok(sample_map)
}


fn get_input_file_samples(samplesheet_filename: &str, lane_index: usize, p7_index: usize, p5_index: usize) -> Result<Vec<SampleMap>, Error> {
  /*
  ** Read samplesheet JSON file.
  */
  let reader = std::fs::File::open(samplesheet_filename).expect(&format!("Error: unable to open file {}", samplesheet_filename));
  let samplesheet_json: serde_json::Value = serde_json::from_reader(reader).unwrap();

  /*
  ** Get sample maps from sample_index_list in JSON file.
  */
  let mut sample_map_vec_all: Vec<SampleMap> = Vec::new();
  for sample_index_map in samplesheet_json["sample_index_list"].as_array().unwrap() {
    sample_map_vec_all.push(deserialize_sample_map(sample_index_map.clone()).unwrap());
  } 

  let mut sample_map_vec_file: Vec<SampleMap> = Vec::new();

  /*
  ** Loop over each sample dictionary in samplesheet_json["sample_index_list"].
  */
  for sample_map in sample_map_vec_all.iter() {
    let (lane_index_vec, rt_index_vec, p7_index_vec, p5_index_vec) = get_sample_index_vecs(sample_map).unwrap();

    /*
    ** Skip if this sample entry does not include 'lane_index'.
    */
    if(!lane_index_vec.iter().any(|i| *i == lane_index)) {
      println!("skip sample/not in lanes:\n{:#?}", sample_map);
      continue;
    }

    /*
    ** Skip if this sample entry does not include 'p7_index'.
    */
    if(!p7_index_vec.clone().iter().any(|i| *i == p7_index)) {
      println!("skip sample/not in p7:\n{:#?}", sample_map);
      continue;
    }

    /*
    ** Skip if this sample entry does not include 'p5_index'.
    */
    if(!p5_index_vec.clone().iter().any(|i| *i == p5_index)) {
      println!("skip sample/not in p5:\n{:#?}", sample_map);
      continue;
    }

    sample_map_vec_file.push((*sample_map).clone());
  }

  Ok(sample_map_vec_file)
}


/// Read a barcode TSV file that has the format
///
/// barcode name\tbarcode sequence
///
/// Arguments:
///- file_path: a reference to the input file path
///
/// Return:
///
/// A hash map: the key is the barcode sequence and the value is
/// the barcode name. The barcode sequence may be String or Vec<u8>,
/// choose below.
///
/// Based on URL: https://stackoverflow.com/questions/78639668/fast-reading-from-a-tsv-file-in-rust
#[derive(Deserialize)]
struct RecordCsv {
  hash_name: String,
  hash_barcode: String,
}

pub fn read_barcode_file(file_path: &str) -> Result<HashMap<String, String>, std::io::Error> {
  let reader = std::fs::File::open(file_path).expect(&format!("Error: unable to open file {}", file_path));
  let mut tsv_reader = ReaderBuilder::new()
                         .has_headers(false)
                         .trim(Trim::Fields)
                         .delimiter(b'\t')
                         .comment(Some(b'#'))
                         .from_reader(reader);

  let mut hash_map: HashMap<String, String> = HashMap::new();

  /*
  ** Read the file and store the names and sequences.
  */
  for result in tsv_reader.deserialize() {
     let record: RecordCsv = result?;
     if(hash_map.contains_key(&record.hash_barcode)) {
       panic!("Error: read_barcode_file: duplicate barcode string in file {:#?}", file_path);
     }
     hash_map.entry(record.hash_barcode.to_owned()).or_insert(record.hash_name.to_string());
  }

  Ok(hash_map)
}


fn make_barcode_map(sample_map_vec: &Vec<SampleMap>, barcode_type: &str, default_filename: &str) -> Result<HashMap<String, usize>, Error> {
  /* 
  ** Find barcode file path.
  ** Check the samplesheet json file for a file path. If zero length, use default path.
  */
  let mut file_name = "".to_string();

  /*
  **  Set barcode file path using either the value in
  **  the samplesheet or the default.
  */
  if(barcode_type == "rt_file") {
    for sample_map in sample_map_vec.iter() {
      let tstr = sample_map.rt_file.trim().to_string();
      if(tstr.len() > 0) {
        file_name = tstr;
        break
      }
    }
  } 
  else if(barcode_type == "ligation_file") {
    for sample_map in sample_map_vec.iter() {
      let tstr = sample_map.ligation_file.trim().to_string();
      if(tstr.len() > 0) {
        file_name = tstr;
        break
      }
    }
  }
  else {
    panic!("Error: make_barcode_id_map: unrecognized barcode_type value.");
  }

  if(file_name.len() == 0) {
    file_name = default_filename.to_string();
  }               

  let barcode_well_map = read_barcode_file(&file_name).unwrap();

  /*
  ** Convert wells to indices.
  */
  let mut barcode_index_map: HashMap<String, usize> = HashMap::with_capacity(MAX_NUM_PLATES * 96 + 1);

  if(barcode_type == "rt_file") {
    let well_index_map = make_well_index_map(MAX_NUM_PLATES * 96 + 1, true, true).unwrap();
    for barcode in barcode_well_map.keys() {
      barcode_index_map.insert(barcode.clone(), well_index_map[&barcode_well_map[barcode]]);
    }
  }
  else
  if(barcode_type == "ligation_file") {
    for barcode in barcode_well_map.keys() {
      let well_index: usize = barcode_well_map[barcode][3..].parse().unwrap();
      barcode_index_map.insert(barcode.clone(), well_index);
    }
  }

  Ok(barcode_index_map)
}


/// Convert base10 barcode index to a base4 index encoded
/// as a String with 'A'=0, 'C'=1, 'G'=2, and 'T'=3.
///  
/// Arguments: 
///
/// max_index: maximum index value in returned vector of strings.
///
/// Return:
/// 
/// Vector of Strings that encode base4 indices. Vector index 0
/// is 'A'.
/// 
/// Note:
///
/// Convert base4 encoded index to base10 int using
///  let es = String; // base4 encoded index
///  let mut digmap: HashMap<char, u32> = HashMap::new();
///  let _ = digmap.insert('A', 0);
///  let _ = digmap.insert('C', 1);
///  let _ = digmap.insert('G', 2);
///  let _ = digmap.insert('T', 3);
///  y = es.char().rev().enumerate().map(|(ii,k)| 4u32.pow((ii) as u32) as u32 * digmap[&k]).sum();
/// 
fn make_index_encoder(max_index: u64) -> Result<Vec<String>, Error> {
  let base4 = BaseCustom::<char>::new("ACGT".chars().collect());

  let mut convert_vec = Vec::<String>::new();
  for iv in (0..max_index+1) {
    convert_vec.push(base4.r#gen(iv));
  }
  Ok(convert_vec)
}


struct RtSampleMaps {
  sample_index_to_sample_name_vec: Vec<String>,
  rt_index_to_sample_index_vec: Vec<usize>
}


/*
** Make a struct of vectors that map sample indices to sample names and 
** a vector that maps rt indices to sample indices.
*/
fn make_rt_sample_maps(sample_map_vec: Vec<SampleMap>) -> Result<RtSampleMaps, Error> {
  /*
  ** Gather distinct sample names and assign index values
  ** starting at 0.
  ** Use HashMap to map names to sample indices. Later,
  ** copy names to a vector accessed by sample indices.
  */
  let mut sample_name_to_sample_index_map: HashMap<String, usize> = HashMap::new();
  let mut num_sample: usize = 1_usize;

  for sample_map in sample_map_vec.clone() {
    let sample_name = sample_map.sample_id;

    if(!sample_name_to_sample_index_map.contains_key(&sample_name)) {
      sample_name_to_sample_index_map.insert(sample_name, num_sample);
      num_sample += 1;
    }
  }

  /*
  ** Make RT index to sample index vector and initialize to zeros.
  */
  let mut rt_index_to_sample_index_vec: Vec<usize> = vec![0; MAX_NUM_PLATES * 96 + 1];


  for sample_map in sample_map_vec.clone() {
    let sample_name = sample_map.sample_id;
    let sample_index = sample_name_to_sample_index_map[&sample_name];

    /*
    ** Get RT indices from sample_map.
    */
    let ranges_parts: Vec<&str> = sample_map.ranges.split(":").collect();
    let rt_index_vec = make_index_vec(ranges_parts[0]).unwrap();

    for rt_index in rt_index_vec {
      if(rt_index_to_sample_index_vec[rt_index] == 0) {
        rt_index_to_sample_index_vec[rt_index] = sample_index;
      }
      else
      if(rt_index_to_sample_index_vec[rt_index] != sample_index) {
        panic!("rt index used for more than one sample");
      }
    }
  }

  /*
  ** Make a sample index to sample name vector.
  */
  let num_sample: usize = sample_name_to_sample_index_map.len();
  let mut sample_index_to_sample_name_vec: Vec<String> = vec!("".to_string(); num_sample + 1);

  sample_index_to_sample_name_vec[0] = "Undetermined".to_string();
  for sample_name in sample_name_to_sample_index_map.keys() {
    let sample_index = sample_name_to_sample_index_map[sample_name];
    sample_index_to_sample_name_vec[sample_index] = sample_name.clone();
  }

  let rt_sample_maps: RtSampleMaps = RtSampleMaps { sample_index_to_sample_name_vec: sample_index_to_sample_name_vec, rt_index_to_sample_index_vec: rt_index_to_sample_index_vec };

  Ok(rt_sample_maps)
}


fn open_bam_writers(rt_sample_maps: &RtSampleMaps, lane_index: usize, pcr7_index: usize, pcr5_index: usize, thread_pool: ThreadPool) -> Result<Vec<Box<rust_htslib::bam::Writer>>, Error> {
  let mut bam_out_vec: Vec<Box<rust_htslib::bam::Writer>> = Vec::new();
  let date = chrono::offset::Local::now();

  for i in 0..(rt_sample_maps.sample_index_to_sample_name_vec.len()) {
      let filename: String = format!("{}-{:03}_{:03}_{:03}-L{:03}.bam",
                               rt_sample_maps.sample_index_to_sample_name_vec[i],
                               lane_index,
                               pcr7_index,
                               pcr5_index,
                               lane_index);
      let mut header = rust_htslib::bam::Header::new();
      header.push_comment(b"Made by cram2bam");
      header.push_comment(format!("cram2bam run date: {}", date.to_string()).as_bytes());
      bam_out_vec.push(Box::new(rust_htslib::bam::Writer::from_path(filename.clone(), &header, rust_htslib::bam::Format::Bam).expect(&format!("Error: unable to open file {}", filename.clone()))));
  } 


  for i in 0..(rt_sample_maps.sample_index_to_sample_name_vec.len()) {
    let _ = bam_out_vec[i].set_thread_pool(&thread_pool);
  }

  Ok(bam_out_vec)
}


#[inline]
pub fn u8_to_str(in_ru8: &[u8]) -> Result<&str, std::str::Utf8Error> {
  return(std::str::from_utf8(in_ru8));
}


fn process_cram(ucram_filename: String,
                sample_map_vec: Vec<SampleMap>,
                rt_to_index_map: HashMap<String, usize>,
                lig_to_index_map: HashMap<String, usize>,
                lane_index: usize,
                pcr7_index: usize,
                pcr5_index: usize,
                number_threads: u32) -> Result<(), Error> {
  /*
  ** Make RT and ligation 'well' to index maps.
  */
  // let rt_well_to_index_map  = make_rt_well_to_index_map(MAX_NUM_PLATES).expect("unable to make RT well map");
  // let lig_well_to_index_map = make_lig_well_to_index_map(MAX_NUM_PLATES).expect("unable to make ligation well map");

  let rt_index_to_well_map = make_rt_index_to_well_map(MAX_NUM_PLATES).expect("unable to make RT index map");

  if(MAX_NUM_PLATES * 96 + 1 >= 4_usize.pow(7)) {
    eprintln!("Error: the maximum number of wells exceeds the largest well");
    eprintln!("       index encodable as a string of bases. You must");
    eprintln!("       decrease the value of MAX_NUM_PLATES in the");
    eprintln!("       rna_rtlig_demux program source file main.rs." );
    panic!("");
  }
  let index_encoder = make_index_encoder((MAX_NUM_PLATES * 96 + 1) as u64).unwrap();

  /*
  ** Make an  htslib thread pool.
  */
  let thread_pool = ThreadPool::new(number_threads).unwrap();

  /*
  ** Set up cram reader.
  */
  let mut record_in = rust_htslib::bam::Record::new();
  let mut cram_reader = Reader::from_path(ucram_filename.clone()).expect(&format!("Error: unable to open file {}", ucram_filename.clone()));
  cram_reader.set_thread_pool(&thread_pool).unwrap();

  let rt_sample_maps = make_rt_sample_maps(sample_map_vec).unwrap();
  let mut bam_writer_vec = open_bam_writers(&rt_sample_maps, lane_index, pcr7_index, pcr5_index, thread_pool).unwrap();

  let sample_index_to_sample_name_vec = rt_sample_maps.sample_index_to_sample_name_vec;
  let rt_index_to_sample_index_vec = rt_sample_maps.rt_index_to_sample_index_vec;

  let ba_string = "ba".to_string();
  let ba_u8 = ba_string.as_bytes();
  let bb_string = "bb".to_string();
  let bb_u8 = bb_string.as_bytes();
  let umi_string = "UM".to_string();
  let umi_u8 = umi_string.as_bytes();

  #[allow(unused_assignments)]
  let mut p7_well_name: String = String::new();
  if(pcr7_index > 0) {
    let (ipl, p7_well) = index_to_well(pcr7_index, true).unwrap();
    p7_well_name = format!("P{:02}-{}", ipl, p7_well);
  }
  else {
    p7_well_name = "none".to_string();
  }
  let mut p5_well_name: String = String::new();
  if(pcr5_index > 0) {
    let (ipl, p5_well) = index_to_well(pcr5_index, true).unwrap();
    p5_well_name = format!("P{:02}-{}", ipl, p5_well);
  }
  else {
    p5_well_name = "none".to_string();
  }

  let p7_index_encoded = index_encoder[pcr7_index].clone();
  let p5_index_encoded = index_encoder[pcr5_index].clone();

  #[allow(unused_assignments)]
  let mut rt_well_name: String = String::new();
  #[allow(unused_assignments)]
  let mut sample_name: String = String::new();

  let mut read_length_counter: Vec<u64> = vec!(0_u64; MAX_READ_LENGTH+1);
  let mut undetermined_counter: u64 = 0;

  /*
  ** Process cram records.
  */
  let mut nrecord: usize = 0;
  while let Some(result) = cram_reader.read(&mut record_in) {
    result.expect("Error reading CRAM file.");

/*
    if(nrecord >= 2) {
      break;
    }
*/

    nrecord += 1;

    /*
    ** Extract RT, ligation, and UMI sequences as Strings. 
    */
    // tag bb: rt barcode
    let Aux::String(rt_barcode) = record_in.aux(bb_u8).unwrap() else {panic!("")};
    // tag ba: lig barcode
    let Aux::String(lig_barcode) = record_in.aux(ba_u8).unwrap() else {panic!("")};
    // tag UM: UMI sequence
    let Aux::String(umi_seq_string) = record_in.aux(umi_u8).unwrap() else {panic!("")};

    let rt_index  = rt_to_index_map[rt_barcode];
    let lig_index = lig_to_index_map[lig_barcode];
    let rt_index_encoded = &index_encoder[rt_index];
    let lig_index_encoded = &index_encoder[lig_index];

    let lig_well_name = format!("LIG{}", lig_index);

    let sample_index = rt_index_to_sample_index_vec[rt_index];
    sample_name = sample_index_to_sample_name_vec[sample_index].clone();
    rt_well_name = rt_index_to_well_map[rt_index].clone();

    if(sample_index == 0) {
      undetermined_counter += 1;
    }

    let read_name: String = format!("{}-P5{}-P7{}_{}|{}|{}|{}|{}_{}|{}", sample_name,
                                                                         p5_well_name,
                                                                         p7_well_name,
                                                                         nrecord,
                                                                         sample_name,
                                                                         p5_well_name,
                                                                         p7_well_name,
                                                                         rt_well_name,
                                                                         lig_well_name,
                                                                         umi_seq_string);

    let cell_barcode_string: String = format!("{:A>7}{:A>7}{:A>7}{:A>7}",
                                                  rt_index_encoded,
                                                  lig_index_encoded,
                                                  p7_index_encoded,
                                                  p5_index_encoded);

    {
      let mut record_out: rust_htslib::bam::Record = rust_htslib::bam::Record::new();

      let seq = record_in.seq().as_bytes();
      read_length_counter[seq.len()] += 1;

      record_out.set(read_name.as_bytes(),
                 None,
                 &seq,
                 record_in.qual());

      record_out.push_aux("CB".as_bytes(), rust_htslib::bam::record::Aux::String(&cell_barcode_string)).expect("Error: unable to add barcode sequence to BAM record tags.");
      record_out.push_aux("CY".as_bytes(), rust_htslib::bam::record::Aux::String("CCCCCCCCCCCCCCCCCCCCCCCCCCCC")).expect("Error: unable to add barcode sequence to BAM record tags.");
      record_out.push_aux("UB".as_bytes(), rust_htslib::bam::record::Aux::String(&umi_seq_string)).expect("Error: unable to add barcode quality values to BAM record tags.");
      record_out.push_aux("UY".as_bytes(), rust_htslib::bam::record::Aux::String("CCCCCCCC")).expect("Error: unable to add barcode quality values to BAM record tags.");

      bam_writer_vec[sample_index].write(&record_out).expect("Error: unable to write record to BAM file.");
    }
  }

  println!("cram2bam");
  println!("  number of cram records processed: {}", nrecord);
  println!("  number of undetermined reads:     {}", undetermined_counter);
  println!("");
  println!("  read length distribution");
  println!("  length       reads");
  let mut read_total: u64 = 0;
  for i in (0..MAX_READ_LENGTH+1) {
    if(read_length_counter[i] > 0) {
      read_total += read_length_counter[i];
    }
  }

  for i in (0..MAX_READ_LENGTH+1) {
    if(read_length_counter[i] > 0) {
      println!("  {:6}       {}  ({:.4})", i, read_length_counter[i], read_length_counter[i] as f64 / read_total as f64);
    }
  }
  println!("");
  println!("  total: {}", read_total);
  println!("");


  Ok(())
}


fn main() {
  /*
  ** Process command line options.
  */
  let cl_options = set_cl_options().unwrap();
  let cl_arg = cl_options.get_matches();

  let ucram_filename: String       = cl_arg.get_one::<String>("ucram_file").unwrap().to_string();
  let samplesheet_filename: String = cl_arg.get_one::<String>("sample_sheet").unwrap().to_string();
  let default_rt_filename: String  = cl_arg.get_one::<String>("rt_barcode_file").unwrap().to_string();
  let default_lig_filename: String = cl_arg.get_one::<String>("ligation_barcode_file").unwrap().to_string();
  let lane_index: usize            = cl_arg.get_one::<String>("lane_index").unwrap().parse().unwrap();
  let pcr7_index: usize            = cl_arg.get_one::<String>("pcr7_index").unwrap().parse().unwrap();
  let pcr5_index: usize            = cl_arg.get_one::<String>("pcr5_index").unwrap().parse().unwrap();
  let number_threads:u32           = *cl_arg.get_one::<u32>("number_threads").unwrap();

  println!("cram file:        {}", ucram_filename);
  println!("samplesheet file: {}", samplesheet_filename);
  println!("default rt file:  {}", default_rt_filename);
  println!("default lig file: {}", default_lig_filename);
  println!("lane index:       {}", lane_index);
  println!("pcr7 index:       {}", pcr7_index);
  println!("pcr5 index:       {}", pcr5_index);
  println!("ncpu:             {}", number_threads);
  println!("");

  /*
  ** Read samplesheet and get 'lane' samples.
  */
  let lane: usize = 1;
  let sample_map_vec: Vec<SampleMap> = get_input_file_samples(&samplesheet_filename, lane, pcr7_index, pcr5_index).unwrap();

  /*
  ** Make RT and ligation barcode to well maps.
  */
  let rt_to_index_map  = make_barcode_map(&sample_map_vec, "rt_file", &default_rt_filename).expect("unable to read RT barcode file");
  let lig_to_index_map = make_barcode_map(&sample_map_vec, "ligation_file", &default_lig_filename).expect("unable to read ligation barcode file");

  process_cram(ucram_filename, sample_map_vec, rt_to_index_map, lig_to_index_map, lane_index, pcr7_index, pcr5_index, number_threads).expect("bad status: process_cram");

}

