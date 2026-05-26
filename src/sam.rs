use crate::search::{SearchResult, CigarOp};
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub struct SamWriter {
    file: File,
    reference_name: String,
    reference_length: usize,
}

impl SamWriter {
    pub fn new(
        output_path: impl AsRef<Path>,
        reference_name: &str,
        reference_length: usize,
    ) -> Result<Self> {
        let file = File::create(output_path)?;
        Ok(SamWriter {
            file,
            reference_name: reference_name.to_string(),
            reference_length,
        })
    }

    pub fn write_header(&mut self) -> Result<()> {
        // SAM header format
        writeln!(self.file, "@HD\tVN:1.0\tSO:coordinate")?;
        writeln!(
            self.file,
            "@SQ\tSN:{}\tLN:{}",
            self.reference_name, self.reference_length
        )?;
        writeln!(self.file, "@PG\tID:rusty-suffix\tPN:rusty-suffix\tVN:0.1.0")?;
        Ok(())
    }

    pub fn write_result(&mut self, result: &SearchResult) -> Result<()> {
        let flag = 0; // Properly mapped (0)
        let pos = result.reference_position + 1; // Convert to 1-based
        let mapq = calculate_mapq(result.mismatches);

        // Generate CIGAR string from alignment operations or fallback to simple format
        let cigar = if let Some(alignment) = &result.alignment {
            build_cigar_string(&alignment.operations)
        } else {
            format!("{}M", result.match_length)
        };

        // Get the matched portion of the query (use the match length from result)
        let query_seq = &result.query_sequence[..result.match_length.min(result.query_sequence.len())];

        // Build the SAM record
        let mut record = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            result.query_id, flag, self.reference_name, pos, mapq, cigar
        );

        // RNEXT, PNEXT, TLEN (standard unmapped pair fields)
        record.push_str("\t*\t0\t0");

        // SEQ - matched portion of query sequence
        let seq = String::from_utf8_lossy(query_seq);
        record.push('\t');
        record.push_str(&seq);

        // QUAL - unknown quality scores
        record.push_str("\t*");

        // Optional tags
        record.push_str(&format!("\tNM:i:{}", result.mismatches));
        record.push_str(&format!("\tAS:i:{}", calculate_alignment_score(result.mismatches, result.match_length)));

        writeln!(self.file, "{}", record)?;
        Ok(())
    }

    pub fn write_results(&mut self, results: &[SearchResult]) -> Result<()> {
        self.write_header()?;
        for result in results {
            self.write_result(result)?;
        }
        Ok(())
    }
}

fn build_cigar_string(operations: &[CigarOp]) -> String {
    if operations.is_empty() {
        return String::new();
    }

    let mut cigar = String::new();
    let mut current_op = operations[0];
    let mut count = 1;

    for &op in &operations[1..] {
        if op == current_op {
            count += 1;
        } else {
            cigar.push_str(&format!("{}{}", count, op_to_char(current_op)));
            current_op = op;
            count = 1;
        }
    }

    // Add the last operation
    cigar.push_str(&format!("{}{}", count, op_to_char(current_op)));
    cigar
}

fn op_to_char(op: CigarOp) -> char {
    match op {
        CigarOp::Match => '=',
        CigarOp::Mismatch => 'X',
        CigarOp::Insertion => 'I',
        CigarOp::Deletion => 'D',
        CigarOp::SoftClip => 'S',
        CigarOp::HardClip => 'H',
    }
}

fn calculate_mapq(mismatches: usize) -> u8 {
    // MAPQ: mapping quality (0-60)
    // 60 for perfect matches, decrease based on mismatches
    match mismatches {
        0 => 60,
        1 => 50,
        2 => 40,
        3 => 30,
        4 => 20,
        5 => 10,
        _ => 1,
    }
}

fn calculate_alignment_score(mismatches: usize, match_length: usize) -> i32 {
    // Simple alignment score: match_length - (mismatch_penalty * mismatches)
    let mismatch_penalty = 4;
    (match_length as i32) - ((mismatches as i32) * mismatch_penalty)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mapq_calculation() {
        assert_eq!(calculate_mapq(0), 60);
        assert_eq!(calculate_mapq(1), 50);
        assert_eq!(calculate_mapq(2), 40);
    }

    #[test]
    fn test_alignment_score() {
        assert_eq!(calculate_alignment_score(0, 20), 20);
        assert_eq!(calculate_alignment_score(2, 20), 12);
    }

    #[test]
    fn test_cigar_string_generation() {
        let ops = vec![
            CigarOp::Match, CigarOp::Match, CigarOp::Mismatch,
            CigarOp::Match, CigarOp::Match
        ];
        let cigar = build_cigar_string(&ops);
        assert_eq!(cigar, "2=1X2=");
    }

    #[test]
    fn test_cigar_with_indels() {
        let ops = vec![
            CigarOp::Match, CigarOp::Match, CigarOp::Insertion,
            CigarOp::Match, CigarOp::Deletion, CigarOp::Match
        ];
        let cigar = build_cigar_string(&ops);
        assert_eq!(cigar, "2=1I1=1D1=");
    }
}
