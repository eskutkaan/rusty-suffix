pub mod config;
pub mod fasta;
pub mod index;
pub mod query;
pub mod sam;
pub mod search;
pub mod table;

pub use config::Config;
pub use fasta::FastaReader;
pub use index::SuffixArrayIndex;
pub use query::QueryBatch;
pub use sam::SamWriter;
pub use search::{ApproximateSearcher, CigarOp, AlignmentDetail};
pub use table::TableWriter;

use anyhow::Result;
use std::path::Path;
use rayon::prelude::*;

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

        let index = SuffixArrayIndex::load_or_build(
            &reference_path,
            config.load_index.as_ref(),
            config.cache_index.as_ref(),
        )?;

        Ok(Self {
            index,
            config,
            reference_name,
            reference_length,
        })
    }

    pub fn search_queries(&self, query_path: impl AsRef<Path>) -> Result<Vec<search::SearchResult>> {
        let reader = FastaReader::new(query_path)?;
        let batches = reader.batch_iterator(self.config.batch_size)?;

        let all_results: Result<Vec<Vec<search::SearchResult>>> = batches
            .into_iter()
            .par_bridge()
            .map(|batch| {
                log::debug!("Processing batch of {} queries", batch.len());

                let searcher = ApproximateSearcher::new(
                    &self.index,
                    self.config.mismatch_tolerance,
                    self.config.min_seed_length,
                );

                searcher.search_batch(&batch)
            })
            .collect();

        let results: Vec<search::SearchResult> = all_results?
            .into_iter()
            .flatten()
            .collect();

        log::info!("Total results: {}", results.len());
        Ok(results)
    }

    pub fn reference_name(&self) -> &str {
        &self.reference_name
    }

    pub fn reference_length(&self) -> usize {
        self.reference_length
    }
}
