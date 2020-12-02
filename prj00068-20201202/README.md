- Definition of MUT column names:

|Column Name                 |Definition|
|----                        |--        |
|contig                      |The reference sequence name
|start                       |0-based start position of the feature in contig
|end                         |Half-open end position of the feature in contig
|sample                      |The sample name
|variation_type              |The category this variant is assigned too
|ref                         |The reference allele at this position
|alt                         |The left-aligned, normalized, alternate allele at this position
|raw_alt_depth               |The raw alternate allele depth as computed by the variant caller before post-processing
|depth                       |The total read depth at this position (including No-calls)
|subtype                     |The substitution type for variants of class “snv”
|context                     |The 3-mer context about this variant’s left-aligned position
|informative_total_depth     |The informative total depth at this position (excluding N-calls)
|vaf                         |The variant allele frequency to 10 decimals
|is_snp                      |Whether we think this variant is a germline polymorphism or not
|is_hom_alt_with_somatic_ref |If this variant is a homozygous alternate, did we observe somatic reference calls at this position?
|final_somatic_alt_depth     |The final post-processing adjusted alternate allele depth (for mutation frequency and spectral analysis)
