# Test Files for rusty-suffix

This directory contains example reference and query files for testing the rusty-suffix tool.

## Files

### reference.fa
A 1008 bp reference sequence containing:
- A start codon (ATG)
- Highly repetitive AGCTAGCT sequences throughout
- Suitable for testing exact and approximate matching

### queries.fa
10 test queries with different characteristics:

| Query | Type | Expected Behavior |
|-------|------|-------------------|
| `exact_match_short` | 9 bp exact match | Multiple hits in repetitive region |
| `exact_match_medium` | 22 bp exact match | Single or few hits |
| `exact_match_long` | 63 bp exact match | Single hit |
| `query_with_1_mismatch` | Contains 1 mismatch | Tests tolerance with mismatch_tolerance ≥ 1 |
| `query_with_2_mismatches` | Contains 2 mismatches | Tests tolerance with mismatch_tolerance ≥ 2 |
| `query_in_repetitive_region` | 28 bp from repetitive region | Many hits expected |
| `query_at_start` | 8 bp from start | Tests boundary conditions |
| `query_at_middle` | From middle repetitive section | Multiple hits |
| `query_partial_match` | 18 bp without mismatches | Single or few hits |
| `query_reverse_complement_like` | From repetitive region | Multiple hits |

## Quick Start

### Test 1: Exact Matching (No Mismatches)
```bash
./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 0 \
  -o results_exact.sam
```

### Test 2: Approximate Matching (Allow Mismatches)
```bash
./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  -o results_approx.sam
```

### Test 3: With Caching (Build and Cache)
```bash
./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  --cache-index cache.bin \
  -o results_cached_build.sam
```

### Test 4: With Caching (Load from Cache)
```bash
# This will be much faster if cache.bin exists from Test 3
./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  --load-index cache.bin \
  -o results_cached_load.sam
```

### Test 5: Parallel Processing
```bash
./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  -t 4 \
  -o results_parallel.sam
```

### Test 6: Verbose Output
```bash
./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  -v \
  -o results_verbose.sam
```

## Expected Results

Running the tool on these files should produce:
- **Exact queries**: 1-many hits depending on repetitive content
- **Mismatches**: Additional hits when tolerating edits
- **Short queries**: Higher hit count due to more matching opportunities
- **Repetitive queries**: Very high hit count (AGCTAGCT repeats often)

## Inspecting Results

View the SAM output:
```bash
head -20 results_exact.sam
```

Check CIGAR strings (should be detailed operations like `8=`, `20=2X1=`):
```bash
grep -v "^@" results_exact.sam | cut -f6 | head -10
```

## Performance Testing

Time the full pipeline:
```bash
time ./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  -o results_perf.sam
```

Verify cache speedup:
```bash
time ./target/release/rusty-suffix \
  -r examples/reference.fa \
  -q examples/queries.fa \
  -m 2 \
  --load-index cache.bin \
  -o results_cached.sam
```

The second run should be measurably faster if loading from cache.
