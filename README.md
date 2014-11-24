Reference Processing
====================

AnnotationSnakefile
-------------------

`AnnotationSnakefile` is a
[Snakemake](https://bitbucket.org/johanneskoester/snakemake) file containing
rules for generating bacterial annotation files from NCBI GFF files.

`extract_ids_from_gff.py` requires [gffutils](https://pypi.python.org/pypi/gffutils).

NCBI GFF downloads: ftp://ftp.ncbi.nlm.nih.gov/genomes/

Input:
#   `{organism}.gff` - GFF downloaded from NCBI
#   `{organism}.fasta` - Fasta sequence
#   `ORGANISMS` - List of organism names

Output:
#   `{organism}.gtf` - GTF file generated by 'gffread' from tophat package
#   `{organism}.bed` - BED12 file generated by converting the GTF file
#   `{organism}_parsed.gff` - GFF file generated by parsing NCBI GFF through 'gffread'
#       Note: GFF parsing sometimes fails with segmentation fault
#   `{organism}_gene_ids_unique_coverage.txt` - Gene ids along with amount of the
#       gene that is unambiguous with respect to other genes on the same strand.
#   `{organism}_gene_annotations.txt` - Gene annotations for each gene id


*Template*

    ORGANISMS = ["escherichia_coli_k12_nc_000913_3",
                 "pseudomonas_aeruginosa_pao1_nc_002516_2",
                 "staphylococcus_aureus_subsp__aureus_str__newman_nc_009641_1"]
   
    rule all:
        input: expand("{organism}.bed", organism=ORGANISMS),
               # GFF Parsing sometimes fails with segmentation fault
               #expand("{organism}_parsed.gff", organism=ORGANISMS),
               # Don't always need gene names file
               # expand("{organism}_gene_names.txt", organism=ORGANISMS),
               expand("{organism}_gene_ids_unique_coverage.txt", organism=ORGANISMS),
               expand("{organism}_gene_annotations.txt", organism=ORGANISMS)
   
   
    REFPROCESSING_DIR = "/Users/lparsons/Documents/projects/reference_processing"
    include: "%s/AnnotationSnakefile" % REFPROCESSING_DIR



Utilities
---------

### `add_seqid_to_gff_id.py`

Add the sequence id to the ID (and Parent) attributes in a GFF file. This is
useful when there are multiple GFF files provided by NCBI that have duplicate
IDs (such as when they include a plasmid with a bacterial sequence).


### `change_gff_id.py`

Change the ID (and Parent) attributes in a GFF file. This is is useful when
there is another suitable identifier to use (such as `locus_tag`). Note that
not all attributes may have this id, so it might be necessary to use *both*
`add_seqid_to_gff_id.py` and `change_gff_id.py` to ensure no duplicates.
