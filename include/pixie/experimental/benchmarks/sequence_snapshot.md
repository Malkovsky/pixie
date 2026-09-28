# Sequence and permutation comparison snapshot

Current comparison run: **2026-09-27**. Experimental implementations remain
experimental; production defaults have not changed. This snapshot retains the
complete registered comparison matrix, including native controls, Immer,
primitive probes, and sorted wrappers. Raw JSON and logs are not retained in
the repository.

## Coverage and validation

| Benchmark target               | Registered | Measured | RAM skips |
| ------------------------------ | ---------- | -------- | --------- |
| permutable_sequence_benchmarks | 1038       | 1022     | 16        |
| permutation_benchmarks         | 124        | 124      | 0         |
| bit_sequence_benchmarks        | 795        | 795      | 0         |
| sorted_sequence_benchmarks     | 27         | 27       | 0         |
| Total                          | 1984       | 1968     | 16        |

Every registered case was attempted. The 16 skipped cases are listed with their
RAM-preflight reason in the detailed tables. The default 4096 MiB configured
budget was retained; the harness applies its additional planning headroom.
Both navigation-beyond-L3 probes were rerun successfully with
`PIXIE_SEQUENCE_L3_MIB=16`. They check navigation-storage size against the
declared cache size; they do not measure hardware cache misses.

Final Release and ASan validation each passed **507 tests**: 221 bit-sequence,
228 permutable-sequence, 49 permutation, and 9 sorted-wrapper tests. The latter
validation includes failed-insertion rollback and distinct-sentinel input
ranges. Benchmark coverage was checked against the exact registered names and
the final executable hashes.

## Reading the results

Timings are medians of three repetitions, rounded to three significant digits.
The runner requested a 0.05 s minimum and randomized repetition order within
each operation/size group. Explicit harness iteration counts still apply.
CPU time is divided by `operations_per_iteration` when that counter exists;
otherwise the tables report one complete benchmark iteration. Construction
uses elements as operations; mixed sorted workloads use insert/query pairs.
Fixture construction, resets, and destruction follow each registered harness's
timing boundaries. Batched probes are not standalone-call latency measurements.

This is a broad screening snapshot, not a precise ranking of close alternatives:
**864 of 1968 measured rows have CPU-time CV above 10%**. The detailed tables
retain CV so noisy rows remain visible. In particular, fresh independent reads
at 32 MiB varied substantially. Rerun a small relevant subset before promoting
a layout or relying on a small percentage difference.

Owned-byte counters describe requested live storage, not RSS or allocator
bookkeeping. Immer's accounting excludes its heap free lists and temporary
persistent snapshots. A missing memory counter is shown as a dash. Packed,
indirect, and large-object reads have different payload semantics; scalar/first
word probes must not be interpreted as full-object reads. Permutation merges
rebase donor indices, whereas permutable-sequence and raw controls preserve
values. Their merge timings are not interchangeable contracts.

## Representative uint64_t results

N = 4,194,304 values (32 MiB useful payload). `Default` is
`Packed64_B256_F8_Prefix`; `F32/R96` is the experimental
`PowerOfTwo64_L256_F32_R96`; `Immer` is `Immer64_F32_L32`.

| Operation                    | Unit       | Default | F32/R96 | Immer |
| ---------------------------- | ---------- | ------- | ------- | ----- |
| Build                        | ns/value   | 5.14    | 2.83    | 3.82  |
| Fresh independent read       | ns/read    | 360     | 59.6    | 35.0  |
| Fresh dependent read         | ns/read    | 320     | 204     | 221   |
| Edited independent read      | ns/read    | 181     | 70.4    | 72.3  |
| Edited dependent read        | ns/read    | 372     | 193     | 214   |
| Random insertion             | us/insert  | 2.22    | 2.53    | 8.62  |
| General rotation             | us/rotate  | 16.9    | 12.5    | 23.5  |
| Whole-sequence rotation      | us/rotate  | 7.94    | 4.51    | 9.37  |
| Local 65-value rotation      | us/rotate  | 2.37    | 1.13    | 12.5  |
| Equal consuming merge        | us/merge   | 6.13    | 10.4    | 11.2  |
| Unequal consuming merge      | us/merge   | 6.61    | 3.37    | 8.35  |

These results support retaining multiple experiments, not choosing a universal
winner. The F32/R96 layout remains a useful mutation-oriented candidate, while
the default also wins some mutation cases. Immer remains a strong read
reference. Small read differences and the noisy fresh-read row need a focused
confirmation before drawing stronger conclusions.

## Why the alternatives remain experimental

- **Power-of-two leaves and F32/R96:** regular-prefix traversal and a bounded
  reserve can reduce mutation overhead. The reserve activates at 1 MiB of
  logical payload and keeps at most 96 internal nodes, never payload leaves.
  Direct F32 nodes occupy 512 bytes here, so the maximum retained-node payload
  is 48 KiB. Smaller owners retain no reserve. This is not the production
  transaction-local allocation policy.
- **Root cache:** the cached F32/R96 variant measured 35.6 ns for fresh
  independent reads versus 59.6 ns without the cache, but 85.2 ns versus
  70.4 ns after edits. The fresh rows have CV above 40%. This run does not
  establish a consistent cache benefit, and the cache adds owner state.
- **Byte child order:** on the F32-aligned rotation probe, F32/R96 byte order
  measured 306 ns versus 401 ns with direct pointers. General rotations were
  17.5 us versus 12.5 us, and whole rotations 5.81 us versus 4.51 us. The map
  is lazy and confined to leaf parents; upper nodes remain direct. It adds an
  indexed load and enlarges an F32 node from 512 to 560 bytes. The current
  split/join paths still materialize child order, so a fast local shuffle is
  not a general mutation win. F64 variants are retained in the full matrix.
- **16-slot packed order and 16x16 grouping:** these remain registered,
  tested alternatives. The grouped implementation uses a uniform fanout of
  256, not separate medium/high fanouts. Its order metadata occupies 176 bytes;
  the integrated node grows from 4096 to 4272 bytes. Whole-group rotations can
  change order metadata, while spill cases still move pointers. The 16-slot
  dependent primitive uses a legacy short cycle; the 32/64 probes use
  full-period address cycles. Do not compare those probes as equal coverage
  of different capacities.
- **Split-bit 32/64 order:** retained as a primitive experiment, not an
  integrated tree policy. Encodings occupy 20/48 bytes and complete objects
  24/56 bytes, versus 33/65 bytes for byte-order objects. Its rotation is
  scalar; byte-order rotation has locally dispatched SIMD and scalar paths.
  Primitive timings exclude tree traversal, measures, and allocation.

The shared production engine retains its promoted traversal, construction,
seam-repair, and rotation improvements. The experimental policy header owns
retained reserves and mapped child order. Production catalogs do not include
experimental implementations.

## Sorted wrappers versus std::set

N = 1,048,576 unique uint64_t keys. Insert timings average construction from
empty to N. Lower-bound timing repeats 128 fixed probes, half hits, with fixture
construction excluded; it is a warm repeated-query workload, not a stream of
fresh random queries. Both wrappers preserve duplicates, but these benchmark
inputs are unique so the comparison remains compatible with `std::set`.

| Workload                 | Unit     | Sequence wrapper | Vector/permutation | std::set |
| ------------------------ | -------- | ---------------- | ------------------ | -------- |
| Random insertion         | ns/key   | 647              | 1540               | 717      |
| Lower bound              | ns/query | 121              | 269                | 33.8     |
| Insert + lower bound     | ns/pair  | 1460             | 2860               | 1440     |

The sequence wrapper is competitive on this insertion workload but remains
behind `std::set` on the repeated lower-bound probe. Mixed work is close here;
this does not establish one best sorted representation.

## Complete current tables

- Permutable sequence: [part 1](permutable_sequence-1.md),
  [part 2](permutable_sequence-2.md), [part 3](permutable_sequence-3.md),
  [part 4](permutable_sequence-4.md).
- Permutation: [all cases](permutation-1.md).
- Bit sequence and primitives: [part 1](bit_sequence-1.md),
  [part 2](bit_sequence-2.md), [part 3](bit_sequence-3.md).
- Sorted wrappers: [all cases](sorted_sequence-1.md).

## Reproduction

Host: AMD Ryzen 7 8845HS, 16 logical CPUs, 16 MiB shared L3, Linux under WSL.
GCC 13.3, C++20, Release `-O3 -DNDEBUG`, `MARCH=native`, AVX-512 enabled
(including the host's VBMI support), Google Benchmark 1.9.4. Timings were pinned
to logical CPU 4 and serialized under one lock; builds and tests did not run
concurrently with measurements. Hardware performance counters were unavailable.
Immer is pinned to `bd4fc749b97dfa2b66a8f3de00bbf234db4856ef` and remains optional.

The measured ad hoc build was `build/immer-review-release`, with tests,
benchmarks, and third-party backends enabled. The supported comparison preset
provides the equivalent benchmark configuration:

```sh
cmake --preset benchmark-all-backends
cmake --build --preset benchmark-all-backends --target \
  permutable_sequence_benchmarks permutation_benchmarks \
  bit_sequence_benchmarks sorted_sequence_benchmarks -j 2
python3 scripts/benchmark_sequences.py \
  --build-dir build/benchmark-all-backends \
  --output-dir /tmp/pixie-sequence-results --cpu 4 \
  --repetitions 3 --min-time 0.05s --resume
env PIXIE_SEQUENCE_L3_MIB=16 python3 scripts/benchmark_sequences.py \
  --build-dir build/benchmark-all-backends \
  --output-dir /tmp/pixie-sequence-results --cpu 4 \
  --binary bit_sequence_benchmarks --filter NavigationBeyondL3 \
  --repetitions 3 --min-time 0.05s --resume
```

Use the actual cache size on another host. The complete run returns a nonzero
status when a RAM guard skips a case; inspect its logs and summary. A filtered
run writes its own invocation summary, while the per-group results remain in
the output directory. For a combined snapshot, select the latest result for
each registered case with the matching executable hash. The runner's
`--list-only`, `--binary`, and `--filter` options allow a small inventory before
an experiment; the full matrix is intentionally a finalization run.

Final executable SHA-256 values:

```text
permutable_sequence_benchmarks 116ee550408bbcaf514ae66c38108659175b89f17622f7b09d4c2059b25945b3
permutation_benchmarks         5952a434df724a040a4ac01e975be1a9a1253d73b06e7f265ffd81c21bc75c9f
bit_sequence_benchmarks        5f282951a113bc2011faa38adbc267579af8b0d2998d2d22f478db3e6fe2a2d5
sorted_sequence_benchmarks     6482189eda3c1bd6384120f1f1a96acf5c5decac628653cea24ce81e83a80168
```
