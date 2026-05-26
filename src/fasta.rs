use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub sequence: Vec<u8>,
}

pub struct FastaReader {
    path: String,
}

impl FastaReader {
    pub fn new(path: impl AsRef<Path>) -> Result<Self> {
        let path = path
            .as_ref()
            .to_str()
            .context("Invalid path")?
            .to_string();
        Ok(FastaReader { path })
    }

    pub fn batch_iterator(&self, batch_size: usize) -> Result<Vec<Vec<Sequence>>> {
        let mut all_batches = Vec::new();
        let mut current_batch = Vec::new();

        let mut reader = parse_fastx_file(&self.path)
            .context("Failed to open FASTA file")?;

        while let Some(record) = reader.next() {
            let record = record.context("Failed to parse FASTA record")?;
            current_batch.push(Sequence {
                id: String::from_utf8_lossy(record.id()).into_owned(),
                sequence: record.seq().to_vec(),
            });

            if current_batch.len() >= batch_size {
                all_batches.push(current_batch);
                current_batch = Vec::new();
            }
        }

        if !current_batch.is_empty() {
            all_batches.push(current_batch);
        }

        Ok(all_batches)
    }

    pub fn read_all(&self) -> Result<Vec<Sequence>> {
        let mut sequences = Vec::new();
        let mut reader = parse_fastx_file(&self.path)
            .context("Failed to open FASTA file")?;

        while let Some(record) = reader.next() {
            let record = record.context("Failed to parse FASTA record")?;
            sequences.push(Sequence {
                id: String::from_utf8_lossy(record.id()).into_owned(),
                sequence: record.seq().to_vec(),
            });
        }

        Ok(sequences)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_fasta_reader() -> Result<()> {
        let mut file = NamedTempFile::new()?;
        writeln!(file, ">seq1\nACGT")?;
        writeln!(file, ">seq2\nTGCA")?;

        let reader = FastaReader::new(file.path())?;
        let seqs = reader.read_all()?;

        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].id, "seq1");
        assert_eq!(seqs[0].sequence, b"ACGT");

        Ok(())
    }
}
