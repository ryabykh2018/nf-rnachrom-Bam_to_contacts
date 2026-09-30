# Bam-to-Contacts

`Bam_to_Contacts_grisha_fix_final.py` converts BAM alignments into tab-separated raw contact tables for downstream `nf-rnachrom` processing.

The script groups BAM records belonging to the same read identifier, removes supplementary alignments, separates primary and secondary alignments, extracts alignment metadata and selected SAM tags, assigns mapping status / pairtype, and splits contacts into a main `Unique_RNA` output and an `Other` output.

The main output is intended for downstream processing by the [EditDistance-CIGAR filter](https://github.com/ryabykh2018/nf-rnachrom-EditDistance_CIGAR_filter).

## Requirements

Python dependency:

```text
pysam
```

Input files must be readable BAM files.

## Supported experiment types

```text
ATA
OTA_PE
OTA_SE
RNA_SEQ_PE
RNA_SEQ_SE
```

## Command line

```bash
python Bam_to_Contacts_grisha_fix_final.py \
-r1 <R1.bam> \
[-r2 <R2.bam>] \
-e <experiment_type> \
[-t TAG [TAG ...]] \
-p <output_prefix>
```

Arguments:

| Argument | Required | Description |
|---|---|---|
| `-r1`, `--rna_bam` | yes | Input R1 BAM file |
| `-r2`, `--dna_bam` | paired modes only | Input R2 BAM file |
| `-e`, `--exp_type` | yes | `ATA`, `OTA_PE`, `OTA_SE`, `RNA_SEQ_PE`, or `RNA_SEQ_SE` |
| `-t`, `--other_tags` | no | Additional SAM tags to extract; default: `NH` |
| `-p`, `--prefix` | yes | Prefix for output files |

Despite the historical argument names `rna_bam` / `dna_bam`, the biological meaning of R1 and R2 depends on the selected experiment type.

## Experiment-specific BAM inputs

| Experiment type | R1 | R2 | R2 required |
|---|---|---|---|
| `ATA` | RNA | DNA | yes |
| `OTA_SE` | DNA | not used | no |
| `OTA_PE` | DNA mate 1 | DNA mate 2 | yes |
| `RNA_SEQ_SE` | RNA | not used | no |
| `RNA_SEQ_PE` | RNA mate 1 | RNA mate 2 | yes |

For `OTA_SE` and `RNA_SEQ_SE`, the second BAM is not opened.

## Paired-BAM requirements

For:

```text
ATA
OTA_PE
RNA_SEQ_PE
```

both BAM files must report queryname sorting in the BAM header:

```text
HD:SO = queryname
```

The script also counts unique `query_name` values in both BAM files and requires the counts to be equal.

Important: this validation compares the **number** of unique read identifiers. It does not verify that the complete R1 and R2 read-ID sets are identical.

## Read grouping

Supplementary alignments are removed before grouping.

Reads are grouped using the suffix returned by:

```python
query_name.split('.')[-1]
```

A BAM record belongs to the current group when its full query name either equals this suffix or ends with:

```text
.<suffix>
```

This behavior should be considered when preparing read names. Query-name sorting and naming consistency are required for correct grouping.

## Alignment handling

### Supplementary alignments

Records with:

```text
is_supplementary = True
```

are ignored completely. They are not written as primary or secondary alignments.

### Primary alignment

A mapped record that is not secondary is used as the primary alignment.

The following seven fields are written:

```text
chr
start
end
strand
cigar
NM
mapq
```

For an unmapped read, all seven fields are:

```text
*
```

If more than one non-secondary mapped record is encountered within the same grouped read component, the later record replaces the previously stored primary alignment. Correct BAM primary-alignment structure is therefore assumed.

### Secondary alignments

Records with:

```text
is_secondary = True
```

are stored in:

```text
*_secondary_alignments
```

One secondary alignment is serialized as:

```text
(chr,start,strand+cigar,NM)
```

Example:

```text
(chr2,12456,-98M2S,1)
```

Multiple secondary alignments are joined with:

```text
;
```

If no secondary alignments are present:

```text
*
```

## Coordinate conventions

### Primary alignments

Primary alignment coordinates are written as 1-based values:

```text
start = read.reference_start + 1
end   = read.reference_end
```

Example:

```text
pysam reference_start = 99
pysam reference_end   = 150

output:
start = 100
end   = 150
```

### Secondary alignments

The current compact secondary-alignment representation uses:

```python
read.reference_start
```

without adding 1.

Therefore, the `start` value inside:

```text
(chr,start,strand+cigar,NM)
```

is currently the 0-based `pysam` `reference_start`, whereas the primary `*_start` fields are 1-based.

This difference reflects the current implementation and should be taken into account by downstream users of `*_secondary_alignments`.

## Strand

The strand field is derived from `is_reverse`:

```text
is_reverse = False -> +
is_reverse = True  -> -
```

## NM and additional SAM tags

`NM` is read from the BAM `NM` tag.

If the tag is absent:

```text
*
```

No replacement NM value is calculated.

Additional SAM tags are controlled by:

```text
-t / --other_tags
```

Default:

```text
NH
```

Example:

```bash
-t NH AS XS
```

Available tags are serialized as:

```text
TAG:"value"
```

and multiple tags are joined with `;`.

Example:

```text
NH:"1";AS:"98";XS:"42"
```

If none of the requested tags are present:

```text
*
```

## Mapping status

Each read component is assigned one mapping status:

```text
U = uniquely mapped
M = multimapped
N = not mapped
```

The current implementation uses the following rules:

```text
no primary alignment / "*" -> N
mapped primary + secondary alignments -> M
mapped primary + no secondary alignments -> U
```

`U/M/N` is therefore based on primary/secondary alignment structure.

It is **not** determined from:

```text
MAPQ
NH
```

## Pairtype

For paired experiments, pairtype is constructed by concatenating R1 and R2 mapping statuses.

Examples:

```text
UU
UM
MU
MM
UN
NU
NN
```

For single-end experiments, pairtype consists only of the R1 status:

```text
U
M
N
```

## Output files

Each run creates:

```text
<prefix>_Unique_RNA.tab.rc
<prefix>_Other.tab.rc
```

Both files contain the same experiment-specific header.

### `<prefix>_Unique_RNA.tab.rc`

Contains the mapping classes selected for downstream processing:

| Experiment type | Pairtypes written to `Unique_RNA` |
|---|---|
| `ATA` | `UU`, `UM` |
| `OTA_PE` | `UU` |
| `OTA_SE` | `U` |
| `RNA_SEQ_PE` | `UU` |
| `RNA_SEQ_SE` | `U` |

For ATA, `UM` is intentionally retained:

```text
U = uniquely mapped RNA
M = multimapped DNA
```

This allows downstream processing of contacts with a uniquely mapped RNA component even when the DNA component is multimapped.

### `<prefix>_Other.tab.rc`

Contains all mapping-status combinations not selected for the experiment-specific `Unique_RNA` output.

Examples:

```text
ATA: MU, MM, UN, NU, ...
OTA_PE: UM, MU, MM, UN, ...
OTA_SE: M, N
RNA_SEQ_PE: UM, MU, MM, ...
RNA_SEQ_SE: M, N
```

## Output schema

Every row contains:

```text
read_id
pairtype
7 primary-alignment fields for R1
7 primary-alignment fields for R2
R1 secondary alignments
R2 secondary alignments
R1 other tags
R2 other tags
```

Missing data are represented as:

```text
*
```

For SE experiments, the complete unused R2 side is filled with `*`.

### ATA header

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

### OTA_PE header

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

### OTA_SE header

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

### RNA_SEQ_PE header

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

### RNA_SEQ_SE header

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

Outputs:

```text
sample_Unique_RNA.tab.rc
sample_Other.tab.rc
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

## Representative pairtype examples

### ATA

```text
RNA: primary only
DNA: primary only
-> UU
-> Unique_RNA
```

```text
RNA: primary only
DNA: primary + secondary alignment(s)
-> UM
-> Unique_RNA
```

```text
RNA: primary + secondary alignment(s)
DNA: primary only
-> MU
-> Other
```

### RNA_SEQ_PE

```text
rna1: primary only
rna2: primary only
-> UU
-> Unique_RNA
```

```text
rna1: primary only
rna2: primary + secondary alignment(s)
-> UM
-> Other
```

## Error conditions

The script terminates with an error when required paired input is missing.

Examples include:

```text
BAM file (-r2) is required for PE RNA-seq experiments
DNA BAM file (-r2) is required for RNA-DNA experiments
```

For `ATA`, `OTA_PE`, and `RNA_SEQ_PE`, the script also terminates if:

- either BAM does not report `HD:SO=queryname`;
- the counts of unique read identifiers in the two BAM files differ;
- a BAM cannot be opened during these checks.

## Performance

Output files use:

```text
buffering = 262144
```

Rows are accumulated separately for `Unique_RNA` and `Other` and flushed in batches of:

```text
100000
```

For SE modes, R2 is not opened, avoiding redundant BAM processing.

## Assumptions and limitations

- BAM records are expected to be organized consistently with queryname-based grouping.
- Supplementary alignments are intentionally ignored.
- Secondary alignments are retained and define `M` status.
- `U/M/N` status is not inferred from MAPQ or `NH`.
- Missing `NM` is reported as `*`.
- The paired-BAM pre-check validates equal counts of unique read IDs, not equality of the complete read-ID sets.
- Primary `*_start` fields are 1-based, but starts embedded in `*_secondary_alignments` are currently 0-based.
- Legacy mapper-specific BWA/STAR/HISAT code paths are commented out and are not active.
- The script assumes the BAM contains a conventional primary-alignment structure; if multiple non-secondary mapped records occur in one grouped component, the last one encountered becomes the stored primary alignment.

## Testing and validation status

No dedicated automated test suite for `Bam_to_Contacts_grisha_fix_final.py` was present in the project materials used to prepare this README.

At minimum, repository-level verification should include:

```bash
python -m py_compile Bam_to_Contacts_grisha_fix_final.py
```

and representative integration runs for all five supported experiment types, checking:

- expected headers;
- U/M/N and pairtype assignment;
- `Unique_RNA` / `Other` routing;
- primary and secondary alignment serialization;
- unmapped records;
- missing NM / requested tags;
- paired-BAM sorting/count validation;
- SE behavior without R2.

