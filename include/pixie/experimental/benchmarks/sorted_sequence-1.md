# sorted_sequence_benchmarks, part 1

[Methodology and conclusions](sequence_snapshot.md). Medians of three repetitions; rounded to three significant digits. An operation is one unit of the benchmark’s `operations_per_iteration` counter (for example, one read, rotation, merge, insertion, insert/query pair, or constructed element); otherwise the unit is one whole timed iteration. CV is the CPU-time coefficient of variation across repetitions. Owned bytes are the harness counter, not process RSS; a dash means unavailable.

| Registered case                                                             | CPU ns/unit | Unit      | CV %  | Owned bytes |
| --------------------------------------------------------------------------- | ----------- | --------- | ----- | ----------- |
| `SortedLowerBound/PermutableSequence64_B256_F8_Prefix/N:1048576`            | 121         | operation | 9.29  | 14138520    |
| `SortedLowerBound/PermutableSequence64_B256_F8_Prefix/N:4096`               | 20.9        | operation | 0.909 | 55832       |
| `SortedLowerBound/PermutableSequence64_B256_F8_Prefix/N:65536`              | 53.5        | operation | 2.89  | 883096      |
| `SortedLowerBound/StdSet64/N:1048576`                                       | 33.8        | operation | 5.34  | —           |
| `SortedLowerBound/StdSet64/N:4096`                                          | 8.84        | operation | 0.836 | —           |
| `SortedLowerBound/StdSet64/N:65536`                                         | 19.5        | operation | 3.56  | —           |
| `SortedLowerBound/VectorPermutation64_B256_F8_Prefix/N:1048576`             | 269         | operation | 9.49  | 26467888    |
| `SortedLowerBound/VectorPermutation64_B256_F8_Prefix/N:4096`                | 29.9        | operation | 0.884 | 103792      |
| `SortedLowerBound/VectorPermutation64_B256_F8_Prefix/N:65536`               | 79.3        | operation | 1.89  | 1648432     |
| `SortedMixedInsertLowerBound/PermutableSequence64_B256_F8_Prefix/N:1048576` | 1.46e+03    | operation | 15.8  | 14138520    |
| `SortedMixedInsertLowerBound/PermutableSequence64_B256_F8_Prefix/N:4096`    | 246         | operation | 1.82  | 55832       |
| `SortedMixedInsertLowerBound/PermutableSequence64_B256_F8_Prefix/N:65536`   | 446         | operation | 1.22  | 883096      |
| `SortedMixedInsertLowerBound/StdSet64/N:1048576`                            | 1.44e+03    | operation | 12.6  | —           |
| `SortedMixedInsertLowerBound/StdSet64/N:4096`                               | 98.1        | operation | 5.68  | —           |
| `SortedMixedInsertLowerBound/StdSet64/N:65536`                              | 248         | operation | 7.29  | —           |
| `SortedMixedInsertLowerBound/VectorPermutation64_B256_F8_Prefix/N:1048576`  | 2.86e+03    | operation | 6.16  | 26467888    |
| `SortedMixedInsertLowerBound/VectorPermutation64_B256_F8_Prefix/N:4096`     | 433         | operation | 5.54  | 103792      |
| `SortedMixedInsertLowerBound/VectorPermutation64_B256_F8_Prefix/N:65536`    | 827         | operation | 7.23  | 1648432     |
| `SortedRandomInsert/PermutableSequence64_B256_F8_Prefix/N:1048576`          | 647         | operation | 1.52  | 14138520    |
| `SortedRandomInsert/PermutableSequence64_B256_F8_Prefix/N:4096`             | 150         | operation | 1.16  | 55832       |
| `SortedRandomInsert/PermutableSequence64_B256_F8_Prefix/N:65536`            | 268         | operation | 3.38  | 883096      |
| `SortedRandomInsert/StdSet64/N:1048576`                                     | 717         | operation | 5.94  | —           |
| `SortedRandomInsert/StdSet64/N:4096`                                        | 52.4        | operation | 2.02  | —           |
| `SortedRandomInsert/StdSet64/N:65536`                                       | 112         | operation | 9.37  | —           |
| `SortedRandomInsert/VectorPermutation64_B256_F8_Prefix/N:1048576`           | 1.54e+03    | operation | 6.25  | 26467888    |
| `SortedRandomInsert/VectorPermutation64_B256_F8_Prefix/N:4096`              | 300         | operation | 1.16  | 103792      |
| `SortedRandomInsert/VectorPermutation64_B256_F8_Prefix/N:65536`             | 538         | operation | 0.919 | 1648432     |
