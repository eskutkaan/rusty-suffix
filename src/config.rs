use clap::Parser;
use serde::{Deserialize, Serialize};
use anyhow::Result;

#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(name = "rusty-suffix")]
#[command(about = "High-throughput suffix array search for genomics with approximate matching")]
pub struct Config {
    #[arg(short, long, help = "Reference genome file (FASTA)")]
    pub reference: String,

    #[arg(short, long, help = "Query sequences file (multiFASTA)")]
    pub queries: String,

    #[arg(short, long, default_value = "results.tsv", help = "Output file")]
    pub output: String,

    #[arg(short = 'm', long, default_value = "2", help = "Maximum mismatches allowed")]
    pub mismatch_tolerance: usize,

    #[arg(short = 's', long, default_value = "20", help = "Minimum seed length for exact matching")]
    pub min_seed_length: usize,

    #[arg(short = 'b', long, default_value = "500", help = "Batch size for query processing")]
    pub batch_size: usize,

    #[arg(short = 't', long, help = "Number of threads for parallel processing")]
    pub threads: Option<usize>,

    #[arg(long, help = "Cache suffix array index to disk")]
    pub cache_index: Option<String>,

    #[arg(long, help = "Load suffix array index from cache")]
    pub load_index: Option<String>,

    #[arg(short = 'v', long, help = "Enable verbose logging")]
    pub verbose: bool,
}

impl Config {
    pub fn validate(&self) -> Result<()> {
        if self.batch_size == 0 {
            anyhow::bail!("batch_size must be greater than 0");
        }

        if self.min_seed_length == 0 {
            anyhow::bail!("min_seed_length must be greater than 0");
        }

        if self.mismatch_tolerance == 0 {
            log::warn!("mismatch_tolerance is 0 - only perfect matches will be found");
        }

        Ok(())
    }
}

impl Default for Config {
    fn default() -> Self {
        Self {
            reference: String::new(),
            queries: String::new(),
            output: "results.tsv".to_string(),
            mismatch_tolerance: 2,
            min_seed_length: 20,
            batch_size: 500,
            threads: None,
            cache_index: None,
            load_index: None,
            verbose: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = Config::default();
        assert_eq!(config.mismatch_tolerance, 2);
        assert_eq!(config.min_seed_length, 20);
        assert_eq!(config.batch_size, 500);
    }
}
