pub mod config;
pub mod fasta;
pub mod index;
pub mod query;
pub mod sam;
pub mod search;

pub use config::Config;
pub use fasta::FastaReader;
pub use index::SuffixArrayIndex;
pub use query::QueryBatch;
pub use sam::SamWriter;
pub use search::ApproximateSearcher;

use anyhow::Result;
use std::path::Path;

pub struct SuffixArraySearcher {
    index: SuffixArrayIndex,
    config: Config,
    reference_name: String,
    reference_length: usize,
}

impl SuffixArraySearcher {
    pub fn new(reference_path: impl AsRef<Path>, config: Config) -> Result<Self> {
        log::info!("Loading reference genome from {:?}", reference_path.as_ref());

        // First, read the reference to get the sequence name
        let reader = FastaReader::new(&reference_path)?;
        let sequences = reader.read_all()?;
        if sequences.is_empty() {
            anyhow::bail!("Reference FASTA file is empty");
        }

        let reference_name = sequences[0].id.clone();
        let reference_length = sequences[0].sequence.len();

        let index = SuffixArrayIndex::build(reference_path)?;

        Ok(Self {
            index,
            config,
            reference_name,
            reference_length,
        })
    }

    pub fn search_queries(&self, query_path: impl AsRef<Path>) -> Result<Vec<search::SearchResult>> {
        let mut all_results = Vec::new();

        let reader = FastaReader::new(query_path)?;
        let batches = reader.batch_iterator(self.config.batch_size)?;

        for batch in batches {
            log::debug!("Processing batch of {} queries", batch.len());

            let searcher = ApproximateSearcher::new(
                &self.index,
                self.config.mismatch_tolerance,
                self.config.min_seed_length,
            );

            let batch_results = searcher.search_batch(&batch)?;
            all_results.extend(batch_results);
        }

        log::info!("Total results: {}", all_results.len());
        Ok(all_results)
    }

    pub fn reference_name(&self) -> &str {
        &self.reference_name
    }

    pub fn reference_length(&self) -> usize {
        self.reference_length
    }
}
