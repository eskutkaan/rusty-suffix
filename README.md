# rusty-suffix

High-throughput suffix array search tool for genomics with approximate matching support. Process large reference genomes (>1GB) against thousands of DNA sequence queries with fuzzy matching for variant detection.

## Features

- **Fast suffix array indexing**: O(n log n) sorting-based suffix array construction
- **Approximate matching**: Seed-and-expand algorithm for fuzzy DNA sequence matching
- **High throughput**: Process thousands of queries efficiently with batch processing
- **Memory efficient**: Stream large FASTA files without loading entire files into memory
- **SAM output format**: Standard alignment format compatible with samtools and other bioinformatics tools
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
  --output <results.sam> \
  --mismatch-tolerance 2 \
  --min-seed-length 20 \
  --batch-size 500 \
  --verbose
```

### Options

- `-r, --reference` - Path to reference genome FASTA file
- `-q, --queries` - Path to query sequences multiFASTA file
- `-o, --output` - Output SAM file (default: results.sam)
- `-m, --mismatch-tolerance` - Maximum mismatches allowed (default: 2)
- `-s, --min-seed-length` - Minimum exact seed length (default: 20)
- `-b, --batch-size` - Query batch size for processing (default: 500)
- `-t, --threads` - Number of threads for parallel processing
- `-v, --verbose` - Enable verbose logging

## Output Format

SAM (Sequence Alignment Map) format compatible with samtools. Each alignment record contains:
- `QNAME` - Query sequence identifier
- `FLAG` - 0 for mapped reads
- `RNAME` - Reference sequence name
- `POS` - 1-based position in reference
- `MAPQ` - Mapping quality (60 for perfect matches, lower with mismatches)
- `CIGAR` - Alignment string (e.g., 6M for 6 matches)
- `RNEXT` - Next reference (*for single-end)
- `PNEXT` - Next position (0 for single-end)
- `TLEN` - Template length (0 for single-end)
- `SEQ` - Query sequence aligned
- `QUAL` - Quality scores (* for unknown)
- `NM:i` - Number of mismatches
- `AS:i` - Alignment score

## Usage with SAMtools

View results:
```bash
samtools view -h results.sam | head -20
```

Convert to BAM (compressed):
```bash
samtools view -b results.sam > results.bam
```

Sort alignments:
```bash
samtools sort results.sam -o results.sorted.sam
```

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
./target/release/rusty-suffix -r ref.fasta -q queries.fasta -o results.sam

# View with samtools
samtools view results.sam
```

## Architecture

- `src/main.rs` - CLI entry point and output handling
- `src/lib.rs` - Main SuffixArraySearcher orchestrator
- `src/index.rs` - SuffixArrayIndex with binary search find_pattern
- `src/search.rs` - ApproximateSearcher with seed-and-expand logic
- `src/fasta.rs` - FastaReader with batch iteration
- `src/query.rs` - QueryBatch management
- `src/config.rs` - Config struct with clap CLI derivation
- `src/sam.rs` - SAM file writer and formatting

## Future Optimizations

- Implement LCP array for faster suffix array queries
- Add mmap support for very large reference genomes (>5GB)
- Parallel query batch processing with rayon
- Caching suffix arrays to disk for reuse
- More detailed CIGAR strings with insertions/deletions
