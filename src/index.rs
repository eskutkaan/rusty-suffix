use anyhow::Context;
use anyhow::Result;
use std::fs;
use std::path::Path;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuffixArrayIndex {
    reference: Vec<u8>,
    suffix_array: Vec<usize>,
    lcp: Vec<usize>,
}

#[derive(Serialize, Deserialize)]
pub struct CacheHeader {
    pub reference_hash: u64,
    pub reference_len: usize,
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

        // Compute LCP array using Kasai's algorithm
        let lcp = Self::compute_lcp(&reference, &suffix_array);
        log::info!("LCP array computed successfully");

        Ok(SuffixArrayIndex {
            reference,
            suffix_array,
            lcp,
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

    fn compute_lcp(reference: &[u8], suffix_array: &[usize]) -> Vec<usize> {
        let n = suffix_array.len();
        let mut lcp = vec![0; n];

        if n == 0 {
            return lcp;
        }

        // Build rank array: rank[i] = position of suffix starting at i in the suffix array
        let mut rank = vec![0; n];
        for (i, &sa_idx) in suffix_array.iter().enumerate() {
            rank[sa_idx] = i;
        }

        // Kasai's algorithm: O(n) LCP computation
        let mut h = 0;
        for i in 0..n {
            if rank[i] > 0 {
                let j = suffix_array[rank[i] - 1];

                while i + h < n && j + h < n && reference[i + h] == reference[j + h] {
                    h += 1;
                }

                lcp[rank[i]] = h;

                if h > 0 {
                    h -= 1;
                }
            }
        }

        lcp
    }

    pub fn lcp(&self) -> &[usize] {
        &self.lcp
    }

    fn compute_hash(data: &[u8]) -> u64 {
        let mut hash: u64 = 5381;
        for &byte in data {
            hash = ((hash << 5).wrapping_add(hash)).wrapping_add(byte as u64);
        }
        hash
    }

    pub fn to_bytes(&self) -> Result<Vec<u8>> {
        let header = CacheHeader {
            reference_hash: Self::compute_hash(&self.reference),
            reference_len: self.reference.len(),
        };

        let header_bytes = bincode::serialize(&header)
            .context("Failed to serialize cache header")?;
        let index_bytes = bincode::serialize(self)
            .context("Failed to serialize suffix array index")?;

        let mut result = Vec::new();
        result.extend_from_slice(&(header_bytes.len() as u32).to_le_bytes());
        result.extend_from_slice(&header_bytes);
        result.extend_from_slice(&index_bytes);

        Ok(result)
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<(Self, CacheHeader)> {
        if bytes.len() < 4 {
            return Err(anyhow::anyhow!("Cache file too small"));
        }

        let header_len = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]) as usize;
        let header_end = 4 + header_len;

        if bytes.len() < header_end {
            return Err(anyhow::anyhow!("Invalid cache file format"));
        }

        let header: CacheHeader = bincode::deserialize(&bytes[4..header_end])
            .context("Failed to deserialize cache header")?;

        let index: SuffixArrayIndex = bincode::deserialize(&bytes[header_end..])
            .context("Failed to deserialize suffix array index")?;

        Ok((index, header))
    }

    pub fn load_or_build(
        reference_path: impl AsRef<Path>,
        cache_read_path: Option<impl AsRef<Path>>,
        cache_write_path: Option<impl AsRef<Path>>,
    ) -> Result<Self> {
        let reference_path = reference_path.as_ref();
        let raw_reference = fs::read(reference_path)
            .context("Failed to read reference file")?;

        // Clean the reference the same way as in build()
        let reference: Vec<u8> = raw_reference
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .map(|&b| b.to_ascii_uppercase())
            .collect();

        let reference_hash = Self::compute_hash(&reference);

        // Try loading from cache
        if let Some(cache_path) = cache_read_path {
            if let Ok(cache_bytes) = fs::read(cache_path) {
                match Self::from_bytes(&cache_bytes) {
                    Ok((index, header)) => {
                        if header.reference_hash == reference_hash && header.reference_len == reference.len() {
                            log::info!("Loaded suffix array index from cache");
                            return Ok(index);
                        } else {
                            log::warn!("Cache validation failed: reference mismatch");
                        }
                    }
                    Err(e) => {
                        log::warn!("Failed to load cache: {}", e);
                    }
                }
            }
        }

        // Build new index
        log::info!("Building suffix array index from reference genome");
        let index = Self::build(reference_path)?;

        // Save to cache if requested
        if let Some(cache_path) = cache_write_path {
            match index.to_bytes() {
                Ok(cache_bytes) => {
                    if let Err(e) = fs::write(cache_path, cache_bytes) {
                        log::warn!("Failed to write cache: {}", e);
                    } else {
                        log::info!("Suffix array index cached to disk");
                    }
                }
                Err(e) => {
                    log::warn!("Failed to serialize index for caching: {}", e);
                }
            }
        }

        Ok(index)
    }

    #[cfg(test)]
    pub fn new_for_testing(reference: Vec<u8>, suffix_array: Vec<usize>) -> Self {
        let lcp = Self::compute_lcp(&reference, &suffix_array);
        SuffixArrayIndex {
            reference,
            suffix_array,
            lcp,
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

    #[test]
    fn test_lcp_array() -> Result<()> {
        let mut file = NamedTempFile::new()?;
        write!(file, "BANANA")?;

        let index = SuffixArrayIndex::build(file.path())?;
        let lcp = index.lcp();

        // LCP array should have same length as reference
        assert_eq!(lcp.len(), 6);

        // First element of LCP should always be 0
        assert_eq!(lcp[0], 0);

        // LCP values should be non-negative and not exceed reference length
        for &val in lcp {
            assert!(val <= 6);
        }

        Ok(())
    }

    #[test]
    fn test_cache_serialization() -> Result<()> {
        let mut file = NamedTempFile::new()?;
        write!(file, "ACGTACGTACGT")?;

        let index = SuffixArrayIndex::build(file.path())?;
        let bytes = index.to_bytes()?;

        let (loaded_index, header) = SuffixArrayIndex::from_bytes(&bytes)?;

        assert_eq!(index.reference(), loaded_index.reference());
        assert_eq!(index.suffix_array(), loaded_index.suffix_array());
        assert_eq!(index.lcp(), loaded_index.lcp());
        assert_eq!(header.reference_len, 12);

        Ok(())
    }

    #[test]
    fn test_load_or_build() -> Result<()> {
        use tempfile::TempDir;

        let mut ref_file = NamedTempFile::new()?;
        write!(ref_file, "ACGTACGTACGT")?;

        let cache_dir = TempDir::new()?;
        let cache_path = cache_dir.path().join("cache.bin");

        // First build - should create cache
        let index1 = SuffixArrayIndex::load_or_build(
            ref_file.path(),
            None::<&std::path::Path>,
            Some(&cache_path),
        )?;

        assert!(cache_path.exists(), "Cache file should be created");

        // Second load - should use cache
        let index2 = SuffixArrayIndex::load_or_build(
            ref_file.path(),
            Some(&cache_path),
            None::<&std::path::Path>,
        )?;

        assert_eq!(index1.reference(), index2.reference());
        assert_eq!(index1.suffix_array(), index2.suffix_array());

        Ok(())
    }
}
