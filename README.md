# Bam-to-Contacts

`Bam_to_Contacts_grisha_fix_final.py` converts BAM alignments into tab-separated raw contact tables for downstream [nf-rnachrom](https://github.com/ilnitsky/nf-rnachrom) processing.

The script groups alignments by read identifier, excludes supplementary alignment records while retaining the remaining non-supplementary alignments of the read, separates primary and secondary alignments, extracts alignment metadata and selected SAM tags, assigns mapping status and pairtype, and splits contacts into the main output and the `Other` output.

## Requirements

Python dependency:

```text
pysam
```

Input files must be readable BAM files.

## Supported experiment types

| Experiment type | R1 | R2 | Main output pairtypes |
|---|---|---|---|
| `ATA` | RNA | DNA | `UU`, `UM` |
| `OTA_SE` | DNA | not used | `U` |
| `OTA_PE` | DNA mate 1 | DNA mate 2 | `UU` |
| `RNA_SEQ_SE` | RNA | not used | `U` |
| `RNA_SEQ_PE` | RNA mate 1 | RNA mate 2 | `UU` |

For `OTA_SE` and `RNA_SEQ_SE`, only R1 is processed.

## Command line

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 <R1.bam> \
[-r2 <R2.bam>] \
-e <experiment_type> \
[-t TAG [TAG ...]] \
-p <output_prefix>
```

| Argument | Required | Description |
|---|---|---|
| `-r1`, `--rna_bam` | yes | Input R1 BAM file |
| `-r2`, `--dna_bam` | paired modes only | Input R2 BAM file |
| `-e`, `--exp_type` | yes | `ATA`, `OTA_PE`, `OTA_SE`, `RNA_SEQ_PE`, or `RNA_SEQ_SE` |
| `-t`, `--other_tags` | no | Additional SAM tags to extract; default: `NH` |
| `-p`, `--prefix` | yes | Output file prefix |

`-r2` is required for `ATA`, `OTA_PE`, and `RNA_SEQ_PE`.

The argument names `--rna_bam` and `--dna_bam` originate from the original ATA-oriented implementation, where R1 corresponds to the RNA part of the contact and R2 to the DNA part. Therefore, for `ATA`, `-r1` must be the RNA BAM and `-r2` the DNA BAM. For other experiment types, the biological meaning of R1 and R2 follows the experiment-specific table above.

## Input requirements

For paired experiment types (`ATA`, `OTA_PE`, and `RNA_SEQ_PE`):

- both BAM files must be sorted by query name;
- both BAM files must contain the same number of unique read identifiers.

The script checks query-name sorting from the BAM header (`HD:SO=queryname`) and compares the number of unique `query_name` values in the two files.

For `OTA_PE` and `RNA_SEQ_PE`, proper-pair selection must be performed upstream, immediately after mapping during BAM preparation. Bam-to-Contacts does not repeat proper-pair validation.

For `OTA_SE` and `RNA_SEQ_SE`, only the R1 BAM is required and processed.

Alternative mappings must be emitted as separate secondary-alignment records in the BAM file. Bam-to-Contacts determines multimapping from BAM secondary records (`is_secondary == True`) and does not use mapper-specific tags containing alternative hits for this purpose.

## Read grouping

Reads are grouped using the suffix after the final dot in `query_name`:

```python
query_name.split('.')[-1]
```

A BAM record is assigned to the current read group when its full `query_name` either matches this suffix or ends with:

```text
.<suffix>
```

Therefore, input BAM files are expected to use read identifiers compatible with this naming convention.

## Alignment handling

Supplementary alignment records are excluded from contact construction, while the read itself is retained and processed using its remaining non-supplementary alignments.

For each read component:

- a mapped non-secondary alignment is stored as the primary alignment;
- secondary alignments are retained separately and used to determine multimapping status;
- unmapped primary-alignment fields are written as `*`.

Primary alignments store chromosome, 1-based start and end coordinates, strand, CIGAR, NM, MAPQ, and selected SAM tags.

Secondary alignments retain chromosome, start position, strand, CIGAR, and NM.

Primary `*_start` coordinates are written as 1-based values, whereas the start position stored inside the compact secondary-alignment representation currently corresponds to the 0-based `pysam reference_start`.

## Mapping status and pairtype

Each read component is assigned one mapping status:

```text
U = uniquely mapped
M = multimapped
N = not mapped
```

The current rules are:

```text
mapped primary alignment, no secondary alignments -> U
mapped primary alignment + secondary alignment(s) -> M
no mapped primary alignment -> N
```

`U/M/N` status is determined from primary/secondary alignment structure and is not inferred from MAPQ or `NH`.

For paired experiments, R1 and R2 statuses are concatenated to form the pairtype:

```text
UU
UM
MU
MM
UN
NU
NN
...
```

For single-end experiments, the pairtype is the single R1 status:

```text
U
M
N
```

## Output files

Each run produces two tab-separated contact tables:

```text
<prefix>_Unique_RNA.tab.rc
<prefix>_Other.tab.rc
```

Both files use the same experiment-specific schema.

### `<prefix>_Unique_RNA.tab.rc`

Contains contacts selected for downstream pipeline processing:

```text
ATA        -> UU, UM
OTA_SE     -> U
OTA_PE     -> UU
RNA_SEQ_SE -> U
RNA_SEQ_PE -> UU
```

For `ATA`, `UM` is retained because the RNA component is uniquely mapped even though the DNA component is multimapped.

### `<prefix>_Other.tab.rc`

Contains all remaining mapping-status combinations that are not selected for downstream analysis.

The `Other` output therefore contains valid reconstructed contact records rather than technically invalid input.

The main output is used as the input contact table for the downstream [EditDistance-CIGAR filter](https://github.com/ryabykh2018/nf-rnachrom-EditDistance_CIGAR_filter).

## Output schema

Each output row contains:

```text
read_id
pairtype
7 primary-alignment fields for R1
7 primary-alignment fields for R2
R1 secondary alignments
R2 secondary alignments
R1 selected SAM tags
R2 selected SAM tags
```

Missing values are written as `*`.

For single-end experiments, all fields corresponding to the unused R2 component are `*`.

### ATA

```text
read_id
ATA_pairtype
rna_chr
rna_start
rna_end
rna_strand
rna_cigar
rna_NM
rna_mapq
dna_chr
dna_start
dna_end
dna_strand
dna_cigar
dna_NM
dna_mapq
rna_secondary_alignments
dna_secondary_alignments
rna_other_tags
dna_other_tags
```

### OTA_SE

```text
read_id
OTA_SE_pairtype
dna1_chr
dna1_start
dna1_end
dna1_strand
dna1_cigar
dna1_NM
dna1_mapq
dna2_chr
dna2_start
dna2_end
dna2_strand
dna2_cigar
dna2_NM
dna2_mapq
dna1_secondary_alignments
dna2_secondary_alignments
dna1_other_tags
dna2_other_tags
```

For `OTA_SE`, all `dna2_*` fields are `*`.

### OTA_PE

```text
read_id
OTA_PE_pairtype
dna1_chr
dna1_start
dna1_end
dna1_strand
dna1_cigar
dna1_NM
dna1_mapq
dna2_chr
dna2_start
dna2_end
dna2_strand
dna2_cigar
dna2_NM
dna2_mapq
dna1_secondary_alignments
dna2_secondary_alignments
dna1_other_tags
dna2_other_tags
```

### RNA_SEQ_SE

```text
read_id
RNAseq_SE_pairtype
rna1_chr
rna1_start
rna1_end
rna1_strand
rna1_cigar
rna1_NM
rna1_mapq
rna2_chr
rna2_start
rna2_end
rna2_strand
rna2_cigar
rna2_NM
rna2_mapq
rna1_secondary_alignments
rna2_secondary_alignments
rna1_other_tags
rna2_other_tags
```

For `RNA_SEQ_SE`, all `rna2_*` fields are `*`.

### RNA_SEQ_PE

```text
read_id
RNAseq_PE_pairtype
rna1_chr
rna1_start
rna1_end
rna1_strand
rna1_cigar
rna1_NM
rna1_mapq
rna2_chr
rna2_start
rna2_end
rna2_strand
rna2_cigar
rna2_NM
rna2_mapq
rna1_secondary_alignments
rna2_secondary_alignments
rna1_other_tags
rna2_other_tags
```

## Usage examples

### ATA

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 RNA.queryname.bam \
-r2 DNA.queryname.bam \
-e ATA \
-t NH \
-p sample
```

### OTA_SE

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 DNA.queryname.bam \
-e OTA_SE \
-t NH \
-p sample
```

### OTA_PE

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 DNA_R1.queryname.bam \
-r2 DNA_R2.queryname.bam \
-e OTA_PE \
-t NH \
-p sample
```

### RNA_SEQ_SE

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 RNA.queryname.bam \
-e RNA_SEQ_SE \
-t NH \
-p sample
```

### RNA_SEQ_PE

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 RNA_R1.queryname.bam \
-r2 RNA_R2.queryname.bam \
-e RNA_SEQ_PE \
-t NH \
-p sample
```

## Assumptions and limitations

- Supplementary alignment records are excluded from contact construction, while the read itself is retained and processed using its remaining non-supplementary alignments.
- Secondary alignments are retained and define multimapping status.
- `U/M/N` status is determined from primary/secondary alignment structure and is not inferred from MAPQ or `NH`.
- `NM` is taken directly from the BAM record and is not recalculated.
- For paired inputs, the script checks query-name sorting and equality of the number of unique read identifiers, but does not verify that the complete R1 and R2 read-ID sets are identical.
- Proper-pair selection for `OTA_PE` and `RNA_SEQ_PE` must be performed upstream.
- Alternative mappings must be emitted as separate secondary-alignment records. Mapper-specific tags containing alternative hits are not used to determine multimapping.
- Primary `*_start` coordinates are 1-based, whereas start positions stored inside `*_secondary_alignments` correspond to the 0-based `pysam reference_start`.

For example, BWA-MEM should be run with `-a` when alternative mappings need to be emitted explicitly as separate secondary-alignment records. Without this option, alternative hits may instead be reported through mapper-specific tags and will not be used by Bam-to-Contacts for multimapping classification.

## Performance

Output files are opened with:

```text
buffering = 262144
```

Records are accumulated separately for the main and `Other` outputs and written in batches of:

```text
100000
```

For `OTA_SE` and `RNA_SEQ_SE`, the second BAM file is not opened.
