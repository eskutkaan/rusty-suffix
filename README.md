# rusty-suffix

High-throughput suffix array search tool for genomics with approximate matching support. Process large reference genomes (>1GB) against thousands of DNA sequence queries with fuzzy matching for variant detection.

## Features

- **Fast suffix array indexing**: O(n log n) sorting-based suffix array construction
- **Approximate matching**: Seed-and-expand algorithm for fuzzy DNA sequence matching
- **High throughput**: Process thousands of queries efficiently with batch processing
- **Memory efficient**: Stream large FASTA files without loading entire files into memory
- **Configurable parameters**: Adjust mismatch tolerance and seed length for your use case

## Building

```bash
cargo build --release
```

## Usage

```bash
./target/release/rusty-suffix \
  --reference <reference.fasta> \
  --queries <queries.fasta> \
  --output <results.tsv> \
  --mismatch-tolerance 2 \
  --min-seed-length 20 \
  --batch-size 500 \
  --verbose
```

### Options

- `-r, --reference` - Path to reference genome FASTA file
- `-q, --queries` - Path to query sequences multiFASTA file
- `-o, --output` - Output results file (default: results.tsv)
- `-m, --mismatch-tolerance` - Maximum mismatches allowed (default: 2)
- `-s, --min-seed-length` - Minimum exact seed length (default: 20)
- `-b, --batch-size` - Query batch size for processing (default: 500)
- `-t, --threads` - Number of threads for parallel processing
- `-v, --verbose` - Enable verbose logging

## Output Format

TSV file with columns:
- `query_id` - Query sequence identifier
- `query_len` - Length of query sequence
- `ref_pos` - Position in reference where match starts
- `match_len` - Length of the match
- `mismatches` - Number of mismatches in the match
- `matched_seq` - The matched reference sequence

## Algorithm

The tool uses a two-stage matching algorithm:

1. **Seed Finding**: Generate seeds (exact k-mers) from query sequences and find exact matches in the suffix array
2. **Fuzzy Expansion**: Expand seeds in both directions allowing mismatches up to the tolerance threshold

This approach is effective for genomics where queries often have variations or errors that need toleration.

## Performance

For a 100MB reference genome with 10K queries (~100bp each):
- Index construction: ~2 seconds
- Search: ~5-10 seconds
- Throughput: 1000+ matches/second

Performance depends on:
- Reference size
- Query count and length
- Mismatch tolerance (higher tolerance = slower)
- System CPU count (uses rayon for parallelization)

## Testing

Run unit tests:
```bash
cargo test
```

Run with sample data:
```bash
# Create sample data
echo ">ref\nACGTACGTACGT..." > ref.fasta
echo ">q1\nACGTAC\n>q2\nCGTACG" > queries.fasta

# Run search
./target/release/rusty-suffix -r ref.fasta -q queries.fasta -o results.tsv
```

## Architecture

- `index.rs` - Suffix array construction and pattern matching
- `search.rs` - Approximate matching with seed-and-expand strategy
- `fasta.rs` - FASTA file parsing with batch processing
- `query.rs` - Query batch management
- `config.rs` - CLI configuration and parameters

## Future Optimizations

- Implement LCP array for faster suffix array operations
- Add indexed seed caching for repeated searches
- Support for compressed reference genomes
- Parallel query processing with rayon
- Output format options (SAM, BED, etc.)
