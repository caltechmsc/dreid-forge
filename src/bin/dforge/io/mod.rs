mod infer;
mod spec;

pub use infer::{
    bio_input as infer_bio_input_format, chem_input as infer_chem_input_format,
    output as infer_output_format,
};
pub use spec::OutputSpec;

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, IsTerminal, Read, Stdin, StdoutLock, Write};
use std::path::Path;

use anyhow::{Context, Result};

/// Returns `true` if stderr is a terminal (interactive).
pub fn stderr_is_tty() -> bool {
    io::stderr().is_terminal()
}

/// Returns `true` if stdin is a terminal (interactive).
pub fn stdin_is_tty() -> bool {
    io::stdin().is_terminal()
}

/// Returns `true` if stdout is a terminal (interactive).
pub fn stdout_is_tty() -> bool {
    io::stdout().is_terminal()
}

pub enum InputSource {
    File(BufReader<File>),
    Stdin(BufReader<Stdin>),
}

impl Read for InputSource {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            InputSource::File(r) => r.read(buf),
            InputSource::Stdin(r) => r.read(buf),
        }
    }
}

impl BufRead for InputSource {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            InputSource::File(r) => r.fill_buf(),
            InputSource::Stdin(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            InputSource::File(r) => r.consume(amt),
            InputSource::Stdin(r) => r.consume(amt),
        }
    }
}

pub fn open_input(path: Option<&Path>) -> Result<InputSource> {
    match path {
        Some(p) => {
            let file = File::open(p)
                .with_context(|| format!("Failed to open input file: {}", p.display()))?;
            Ok(InputSource::File(BufReader::new(file)))
        }
        None => Ok(InputSource::Stdin(BufReader::new(io::stdin()))),
    }
}

pub enum OutputTarget {
    File(BufWriter<File>),
    Stdout(BufWriter<StdoutLock<'static>>),
}

impl Write for OutputTarget {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutputTarget::File(w) => w.write(buf),
            OutputTarget::Stdout(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputTarget::File(w) => w.flush(),
            OutputTarget::Stdout(w) => w.flush(),
        }
    }
}

pub fn create_output(path: Option<&Path>) -> Result<OutputTarget> {
    match path {
        Some(p) => {
            let file = File::create(p)
                .with_context(|| format!("Failed to create output file: {}", p.display()))?;
            Ok(OutputTarget::File(BufWriter::new(file)))
        }
        None => Ok(OutputTarget::Stdout(BufWriter::new(io::stdout().lock()))),
    }
}
