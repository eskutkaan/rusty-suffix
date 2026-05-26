use crate::search::{SearchResult, CigarOp};
use anyhow::Result;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub struct TableWriter {
    file: File,
    reference_name: String,
}

impl TableWriter {
    pub fn new(
        output_path: impl AsRef<Path>,
        reference_name: &str,
    ) -> Result<Self> {
        let file = File::create(output_path)?;
        Ok(TableWriter {
            file,
            reference_name: reference_name.to_string(),
        })
    }

    pub fn write_results(&mut self, results: &[SearchResult]) -> Result<()> {
        self.write_header()?;
        for result in results {
            self.write_result(result)?;
        }
        Ok(())
    }

    fn write_header(&mut self) -> Result<()> {
        writeln!(
            self.file,
            "Query_ID\tQuery_Len\tRef_Name\tRef_Pos\tQuery_Start\tQuery_End\tRef_Start\tRef_End\tAlignment_Len\tIdentity%\tMismatches\tGap_Opens\tAlignment_Score\tQuery_Seq\tRef_Seq\tAlignment_Visual"
        )?;
        Ok(())
    }

    fn write_result(&mut self, result: &SearchResult) -> Result<()> {
        let query_len = result.query_sequence.len();
        let ref_pos = result.reference_position + 1; // Convert to 1-based

        // Calculate alignment positions
        let query_start = if let Some(align) = &result.alignment {
            1 + align.query_start_clipped
        } else {
            1
        };

        let query_end = if let Some(align) = &result.alignment {
            query_len - align.query_end_clipped
        } else {
            query_len
        };

        let ref_start = ref_pos;
        let ref_end = ref_pos + result.match_length as usize - 1;

        // Calculate identity percentage
        let identity = if result.match_length > 0 {
            ((result.match_length - result.mismatches) as f64 / result.match_length as f64) * 100.0
        } else {
            0.0
        };

        // Count gap opens
        let gap_opens = count_gap_opens(&result.alignment);

        // Calculate alignment score
        let alignment_score = calculate_alignment_score(result.mismatches, result.match_length);

        // Extract query sequence
        let query_seq = String::from_utf8_lossy(&result.query_sequence[query_start as usize - 1..query_end as usize]);

        // Reference sequence
        let ref_seq = String::from_utf8_lossy(&result.matched_sequence);

        // Generate ASCII alignment visualization
        let alignment_visual = generate_alignment_ascii(&result.alignment, &query_seq, &ref_seq);

        writeln!(
            self.file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}\t{}\t{}",
            result.query_id,
            query_len,
            self.reference_name,
            ref_pos,
            query_start,
            query_end,
            ref_start,
            ref_end,
            result.match_length,
            identity,
            result.mismatches,
            gap_opens,
            alignment_score,
            query_seq,
            ref_seq,
            alignment_visual
        )?;

        Ok(())
    }
}

fn generate_alignment_ascii(
    alignment: &Option<crate::search::AlignmentDetail>,
    query_seq: &str,
    ref_seq: &str,
) -> String {
    match alignment {
        None => String::new(),
        Some(align) => {
            let mut match_line = String::new();
            let query_bytes = query_seq.as_bytes();
            let ref_bytes = ref_seq.as_bytes();

            // Build match line character by character
            for (i, op) in align.operations.iter().enumerate() {
                match op {
                    CigarOp::Match => {
                        if i < query_bytes.len() && i < ref_bytes.len() {
                            if query_bytes[i] == ref_bytes[i] {
                                match_line.push('=');
                            } else {
                                match_line.push('X');
                            }
                        }
                    }
                    CigarOp::Mismatch => match_line.push('X'),
                    CigarOp::Insertion => match_line.push('-'),
                    CigarOp::Deletion => match_line.push('-'),
                    CigarOp::SoftClip => match_line.push('.'),
                    CigarOp::HardClip => match_line.push('.'),
                }
            }

            // Format as 3-line alignment with escaped newlines for TSV
            format!("{}\\n{}\\n{}", query_seq, match_line, ref_seq)
        }
    }
}

fn count_gap_opens(alignment: &Option<crate::search::AlignmentDetail>) -> usize {
    match alignment {
        None => 0,
        Some(align) => {
            let mut gap_opens = 0;
            let mut in_gap = false;

            for op in &align.operations {
                match op {
                    CigarOp::Insertion | CigarOp::Deletion => {
                        if !in_gap {
                            gap_opens += 1;
                            in_gap = true;
                        }
                    }
                    _ => in_gap = false,
                }
            }

            gap_opens
        }
    }
}

fn calculate_alignment_score(mismatches: usize, match_length: usize) -> i32 {
    let mismatch_penalty = 4;
    (match_length as i32) - ((mismatches as i32) * mismatch_penalty)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_count_gap_opens() {
        let align = crate::search::AlignmentDetail {
            operations: vec![
                CigarOp::Match,
                CigarOp::Match,
                CigarOp::Insertion,
                CigarOp::Insertion,
                CigarOp::Match,
                CigarOp::Deletion,
                CigarOp::Match,
            ],
            query_start_clipped: 0,
            query_end_clipped: 0,
        };

        assert_eq!(count_gap_opens(&Some(align)), 2);
    }

    #[test]
    fn test_alignment_score() {
        assert_eq!(calculate_alignment_score(0, 20), 20);
        assert_eq!(calculate_alignment_score(2, 20), 12);
    }
}
