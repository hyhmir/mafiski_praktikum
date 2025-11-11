use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::process;
use pyo3::prelude::*;

/// Read WAV file samples and write them to a text file
#[pyfunction]
pub fn branje(args: Vec<String>) {
    if args.len() != 2 {
        eprintln!("Usage: branje <input.wav> <output.txt>");
        process::exit(1);
    }

    let input_path = &args[0];
    let output_path = &args[1];

    // Open input WAV file
    let file = match File::open(input_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening file {}: {}", input_path, e);
            process::exit(1);
        }
    };
    let mut reader = BufReader::new(file);

    // Create output file
    let mut outfile = match File::create(output_path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error creating file {}: {}", output_path, e);
            process::exit(1);
        }
    };

    // sliding window to find "data"
    let mut window = [0u8; 4];
    let mut idx = 0usize;
    loop {
        let mut byte = [0u8; 1];
        if let Err(_) = reader.read_exact(&mut byte) {
            eprintln!("Reached EOF before finding 'data' chunk");
            process::exit(1);
        }

        if idx < 4 {
            window[idx] = byte[0];
        } else {
            window[0] = window[1];
            window[1] = window[2];
            window[2] = window[3];
            window[3] = byte[0];

            if &window == b"data" {
                break;
            }
        }
        idx += 1;
    }

    // skip 4 bytes (chunk size after "data")
    let mut skip = [0u8; 4];
    if let Err(_) = reader.read_exact(&mut skip) {
        eprintln!("Unexpected EOF while skipping chunk size");
        process::exit(1);
    }

    // read 2-byte little-endian signed samples and write them to the file
    let mut sample_bytes = [0u8; 2];
    loop {
        match reader.read_exact(&mut sample_bytes) {
            Ok(()) => {
                let sample = i16::from_le_bytes(sample_bytes);
                if let Err(e) = writeln!(outfile, "{}", sample) {
                    eprintln!("Error writing to file: {}", e);
                    process::exit(1);
                }
            }
            Err(ref e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
            Err(e) => {
                eprintln!("I/O error while reading samples: {}", e);
                process::exit(1);
            }
        }
    }

    println!("Wrote samples from '{}' to '{}'", input_path, output_path);
}
