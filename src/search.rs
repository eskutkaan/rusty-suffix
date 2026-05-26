use crate::fasta::Sequence;
use crate::index::SuffixArrayIndex;
use anyhow::Result;
use bio::alignment::distance::levenshtein;
use rayon::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    Match,      // = (match)
    Mismatch,   // X (mismatch)
    Insertion,  // I (query has extra bp)
    Deletion,   // D (reference has extra bp)
    SoftClip,   // S (unaligned query bases)
    HardClip,   // H (hard clipped bases)
}

#[derive(Debug, Clone)]
pub struct AlignmentDetail {
    pub operations: Vec<CigarOp>,
    pub query_start_clipped: usize,
    pub query_end_clipped: usize,
}

#[derive(Debug, Clone)]
pub struct SearchResult {
    pub query_id: String,
    pub query_sequence: Vec<u8>,
    pub reference_position: usize,
    pub matched_sequence: Vec<u8>,
    pub mismatches: usize,
    pub match_length: usize,
    pub alignment: Option<AlignmentDetail>,
}

pub struct ApproximateSearcher<'a> {
    index: &'a SuffixArrayIndex,
    mismatch_tolerance: usize,
    min_seed_length: usize,
}

impl<'a> ApproximateSearcher<'a> {
    pub fn new(
        index: &'a SuffixArrayIndex,
        mismatch_tolerance: usize,
        min_seed_length: usize,
    ) -> Self {
        ApproximateSearcher {
            index,
            mismatch_tolerance,
            min_seed_length,
        }
    }

    pub fn search_batch(&self, queries: &[Sequence]) -> Result<Vec<SearchResult>> {
        let results: Result<Vec<Vec<SearchResult>>> = queries
            .par_iter()
            .map(|query| self.search_single(query))
            .collect();

        Ok(results?.into_iter().flatten().collect())
    }

    pub fn search_single(&self, query: &Sequence) -> Result<Vec<SearchResult>> {
        let query_seq = &query.sequence;

        // Generate seeds from the query
        let seeds = self.generate_seeds(query_seq);

        // Find all seed matches
        let seed_matches: Vec<_> = seeds
            .iter()
            .flat_map(|(seed_pos, seed)| {
                let exact_matches = self.index.find_pattern(seed);
                exact_matches
                    .into_iter()
                    .map(|ref_pos| (ref_pos, *seed_pos, seed.len()))
                    .collect::<Vec<_>>()
            })
            .collect();

        // Expand seeds with fuzzy matching (parallelized)
        let results: Result<Vec<Option<SearchResult>>> = seed_matches
            .par_iter()
            .map(|(ref_pos, query_seed_pos, seed_len)| {
                self.expand_seed(query, *ref_pos, *query_seed_pos, *seed_len)
            })
            .collect();

        Ok(results?
            .into_iter()
            .filter_map(|r| r)
            .collect())
    }

    fn generate_seeds(&self, query: &[u8]) -> Vec<(usize, Vec<u8>)> {
        let mut seeds = Vec::new();
        let step = self.min_seed_length.max(1);

        for i in (0..query.len()).step_by(step) {
            let end = (i + self.min_seed_length).min(query.len());
            if end > i {
                seeds.push((i, query[i..end].to_vec()));
            }
        }

        seeds
    }

    fn expand_seed(
        &self,
        query: &Sequence,
        ref_pos: usize,
        query_seed_pos: usize,
        seed_len: usize,
    ) -> Result<Option<SearchResult>> {
        let query_seq = &query.sequence;
        let ref_seq = self.index.reference();

        // Expand left from seed
        let mut query_start = query_seed_pos;
        let mut ref_start = ref_pos;

        while query_start > 0 && ref_start > 0 && query_seq[query_start - 1] == ref_seq[ref_start - 1]
        {
            query_start -= 1;
            ref_start -= 1;
        }

        // Expand right from seed
        let query_right_bound = query_seq.len();
        let ref_right_bound = ref_seq.len();
        let mut query_end = query_seed_pos + seed_len;
        let mut ref_end = ref_pos + seed_len;

        while query_end < query_right_bound
            && ref_end < ref_right_bound
            && query_seq[query_end] == ref_seq[ref_end]
        {
            query_end += 1;
            ref_end += 1;
        }

        // Try to extend with fuzzy matching
        let (final_query_end, final_ref_end) =
            self.fuzzy_extend_right(query_seq, query_end, ref_seq, ref_end)?;
        let (final_query_start, final_ref_start) =
            self.fuzzy_extend_left(query_seq, query_start, ref_seq, ref_start)?;

        let matched_query = &query_seq[final_query_start..final_query_end];
        let matched_ref = &ref_seq[final_ref_start..final_ref_end];

        // Calculate edit distance
        let distance = levenshtein(matched_query, matched_ref) as usize;

        if distance <= self.mismatch_tolerance {
            let operations = Self::compute_cigar_operations(matched_query, matched_ref);
            let alignment = Some(AlignmentDetail {
                operations,
                query_start_clipped: final_query_start,
                query_end_clipped: query_seq.len() - final_query_end,
            });

            Ok(Some(SearchResult {
                query_id: query.id.clone(),
                query_sequence: query_seq.clone(),
                reference_position: final_ref_start,
                matched_sequence: matched_ref.to_vec(),
                mismatches: distance,
                match_length: final_query_end - final_query_start,
                alignment,
            }))
        } else {
            Ok(None)
        }
    }

    fn fuzzy_extend_right(
        &self,
        query: &[u8],
        mut query_pos: usize,
        reference: &[u8],
        mut ref_pos: usize,
    ) -> Result<(usize, usize)> {
        let mut mismatches = 0;
        let extend_limit = self.mismatch_tolerance / 2;

        while query_pos < query.len()
            && ref_pos < reference.len()
            && mismatches < extend_limit
        {
            if query[query_pos] != reference[ref_pos] {
                mismatches += 1;
            }
            query_pos += 1;
            ref_pos += 1;
        }

        Ok((query_pos, ref_pos))
    }

    fn fuzzy_extend_left(
        &self,
        query: &[u8],
        mut query_pos: usize,
        reference: &[u8],
        mut ref_pos: usize,
    ) -> Result<(usize, usize)> {
        let mut mismatches = 0;
        let extend_limit = self.mismatch_tolerance / 2;

        while query_pos > 0 && ref_pos > 0 && mismatches < extend_limit {
            if query[query_pos - 1] != reference[ref_pos - 1] {
                mismatches += 1;
            }
            query_pos -= 1;
            ref_pos -= 1;
        }

        Ok((query_pos, ref_pos))
    }

    fn compute_cigar_operations(
        query: &[u8],
        reference: &[u8],
    ) -> Vec<CigarOp> {
        let mut operations = Vec::new();

        if query.is_empty() || reference.is_empty() {
            return operations;
        }

        // Use simple character-by-character comparison
        // For full alignment, we would use Smith-Waterman, but for suffix array search
        // this provides fast approximate alignment operations
        let min_len = query.len().min(reference.len());

        for i in 0..min_len {
            if query[i] == reference[i] {
                operations.push(CigarOp::Match);
            } else {
                operations.push(CigarOp::Mismatch);
            }
        }

        // Handle length differences
        if query.len() > reference.len() {
            for _ in reference.len()..query.len() {
                operations.push(CigarOp::Insertion);
            }
        } else if reference.len() > query.len() {
            for _ in query.len()..reference.len() {
                operations.push(CigarOp::Deletion);
            }
        }

        operations
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seed_generation() {
        let index = SuffixArrayIndex::new_for_testing(
            b"ACGTACGT".to_vec(),
            vec![0, 1, 2, 3, 4, 5, 6, 7],
        );
        let searcher = ApproximateSearcher::new(&index, 1, 4);
        let seeds = searcher.generate_seeds(b"ACGTACGT");

        assert!(!seeds.is_empty());
        assert!(seeds[0].1.len() >= 4);
    }
}
