use anyhow::Context;
use anyhow::Result;
use std::fs;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct SuffixArrayIndex {
    reference: Vec<u8>,
    suffix_array: Vec<usize>,
}

impl SuffixArrayIndex {
    pub fn build(reference_path: impl AsRef<Path>) -> Result<Self> {
        log::info!("Building suffix array index from reference genome");
        let reference = fs::read(reference_path)
            .context("Failed to read reference file")?;

        // Convert to uppercase and remove newlines for genomics data
        let reference: Vec<u8> = reference
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .map(|&b| b.to_ascii_uppercase())
            .collect();

        log::info!("Reference size: {} bytes", reference.len());

        // Build suffix array using sorting
        let mut suffix_array: Vec<usize> = (0..reference.len()).collect();
        suffix_array.sort_by(|&a, &b| reference[a..].cmp(&reference[b..]));

        log::info!("Suffix array constructed successfully");

        Ok(SuffixArrayIndex {
            reference,
            suffix_array,
        })
    }

    pub fn reference(&self) -> &[u8] {
        &self.reference
    }

    pub fn suffix_array(&self) -> &[usize] {
        &self.suffix_array
    }

    pub fn len(&self) -> usize {
        self.reference.len()
    }

    pub fn is_empty(&self) -> bool {
        self.reference.is_empty()
    }

    /// Find all occurrences of a pattern using binary search on the suffix array
    pub fn find_pattern(&self, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || pattern.len() > self.reference.len() {
            return Vec::new();
        }

        let pattern_upper: Vec<u8> = pattern.iter().map(|&b| b.to_ascii_uppercase()).collect();
        let mut matches = Vec::new();

        // Binary search for left boundary
        let left = self.binary_search_left(&pattern_upper);
        if left >= self.suffix_array.len() {
            return matches;
        }

        // Collect all matches from left boundary
        for &sa_idx in &self.suffix_array[left..] {
            let remaining = self.reference.len().saturating_sub(sa_idx);
            if remaining < pattern_upper.len() {
                break;
            }

            if &self.reference[sa_idx..sa_idx + pattern_upper.len()] == pattern_upper.as_slice()
            {
                matches.push(sa_idx);
            } else {
                break;
            }
        }

        // Sort matches by position for consistency
        matches.sort();
        matches
    }

    fn binary_search_left(&self, pattern: &[u8]) -> usize {
        let mut left = 0;
        let mut right = self.suffix_array.len();

        while left < right {
            let mid = left + (right - left) / 2;
            let sa_idx = self.suffix_array[mid];
            let remaining = self.reference.len().saturating_sub(sa_idx);

            let cmp_len = pattern.len().min(remaining);
            let cmp = &self.reference[sa_idx..sa_idx + cmp_len].cmp(pattern);

            match cmp {
                std::cmp::Ordering::Less => left = mid + 1,
                _ => right = mid,
            }
        }

        left
    }

    pub fn get_suffix(&self, sa_index: usize) -> &[u8] {
        let start = self.suffix_array[sa_index];
        &self.reference[start..]
    }

    #[cfg(test)]
    pub fn new_for_testing(reference: Vec<u8>, suffix_array: Vec<usize>) -> Self {
        SuffixArrayIndex {
            reference,
            suffix_array,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_suffix_array_build() -> Result<()> {
        let mut file = NamedTempFile::new()?;
        write!(file, "ACGTACGTACGT")?;

        let index = SuffixArrayIndex::build(file.path())?;
        assert_eq!(index.len(), 12);
        assert!(!index.is_empty());

        Ok(())
    }

    #[test]
    fn test_pattern_matching() -> Result<()> {
        let mut file = NamedTempFile::new()?;
        write!(file, "ACGTACGTACGT")?;

        let index = SuffixArrayIndex::build(file.path())?;
        let matches = index.find_pattern(b"ACG");

        assert!(!matches.is_empty());
        assert_eq!(matches[0], 0);

        Ok(())
    }
}
