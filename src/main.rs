use anyhow::Result;
use clap::Parser;
use rusty_suffix::{Config, SuffixArraySearcher, SamWriter, TableWriter, filter_top_alignments};
use std::time::Instant;

fn main() -> Result<()> {
    let config = Config::parse();

    // Validate configuration parameters
    config.validate()?;

    if config.verbose {
        env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Debug)
            .init();
    } else {
        env_logger::Builder::from_default_env()
            .filter_level(log::LevelFilter::Info)
            .init();
    }

    log::info!("rusty-suffix: High-throughput suffix array search for genomics");
    log::info!("Reference: {}", config.reference);
    log::info!("Queries: {}", config.queries);
    log::info!("Mismatch tolerance: {}", config.mismatch_tolerance);
    log::info!("Min seed length: {}", config.min_seed_length);

    let start = Instant::now();

    log::info!("Initializing suffix array searcher...");
    let searcher = SuffixArraySearcher::new(&config.reference, config.clone())?;

    let index_time = start.elapsed();
    log::info!("Index construction time: {:?}", index_time);

    let search_start = Instant::now();
    log::info!("Searching queries...");
    let results = searcher.search_queries(&config.queries)?;
    let search_time = search_start.elapsed();

    log::info!("Search completed in {:?}", search_time);
    log::info!("Found {} matches", results.len());

    // Filter to top N alignments per query if requested
    let results = filter_top_alignments(results, config.max_alignments);

    if let Some(max_n) = config.max_alignments {
        log::info!("Filtered to top {} alignments per query: {} matches remain", max_n, results.len());
    }

    // Write results to SAM file
    let mut writer = SamWriter::new(
        &config.output,
        searcher.reference_name(),
        searcher.reference_length(),
    )?;
    writer.write_results(&results)?;

    // Write results to table file if requested
    if let Some(table_path) = &config.table_output {
        let mut table_writer = TableWriter::new(
            table_path,
            searcher.reference_name(),
        )?;
        table_writer.write_results(&results)?;
        log::info!("Table output written to: {}", table_path);
    }

    // Print performance metrics
    let total_time = start.elapsed();
    let throughput = results.len() as f64 / total_time.as_secs_f64();

    println!("\n=== Performance Metrics ===");
    println!("Index construction: {:?}", index_time);
    println!("Search time: {:?}", search_time);
    println!("Total time: {:?}", total_time);
    println!("Throughput: {:.2} matches/sec", throughput);
    println!("Output file: {} (SAM format)", config.output);

    Ok(())
}
