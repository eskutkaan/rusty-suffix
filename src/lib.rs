pub mod config;
pub mod fasta;
pub mod index;
pub mod query;
pub mod search;

pub use config::Config;
pub use fasta::FastaReader;
pub use index::SuffixArrayIndex;
pub use query::QueryBatch;
pub use search::ApproximateSearcher;

use anyhow::Result;
use std::path::Path;

pub struct SuffixArraySearcher {
    index: SuffixArrayIndex,
    config: Config,
}

impl SuffixArraySearcher {
    pub fn new(reference_path: impl AsRef<Path>, config: Config) -> Result<Self> {
        log::info!("Loading reference genome from {:?}", reference_path.as_ref());
        let index = SuffixArrayIndex::build(reference_path)?;

        Ok(Self { index, config })
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
}
