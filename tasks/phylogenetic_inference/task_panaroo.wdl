version 1.0

task panaroo {

meta {
    version: "0.1.0"
    description: "Panaroo pangenome construction task for Core Gene SNP pipeline"
  }

  input {

    Array[File] gff3
    String cluster_name

    # Alignment options
    # align = true is required to produce core/pan alignment FASTAs
    Boolean align = true
    String alignment_type = "core"                      # Options: "core", "pan" — passed to --alignment flag
                                                        # Note: "pan" also produces core_gene_alignment.aln
                                                        # To get both alignments, run twice with different alignment_type values

    # Core genome definition
    Float core_threshold = 0.98                         # Proportion of isolates a gene must be in to be considered core
    String? core_subset                                 # Optional: define core relative to a subset of isolates

    # Pangenome construction behavior
    String clean_mode = "strict"                        # Options: strict, moderate, sensitive
    Float seq_id = 0.95                                 # Sequence identity threshold for initial clustering [-c]
    Float family_threshold = 0.7                        # Family-level sequence identity threshold [-f]
    Float len_dif_percent = 0.98                        # Maximum length difference to cluster sequences
    Float family_len_dif_percent = 0.98                 # Family-level length difference threshold 

    # Graph and edge parameters
    Boolean clean_edges = true                          # Clean edges in the pangenome graph [--no_clean_edges if false]
    Float edge_support_threshold = 0.0                  # Minimum edge support threshold
    Int min_edge_support_sv = 5                         # Minimum number of isolates sharing an edge for structural variant calls
    Int min_trailing_support = 2                        # Minimum cluster size to keep a gene called at the end of a contig
    Int trailing_recursive = 2                          # Number of recursive trailing removals
    Float length_outlier_support_proportion = 0.01      # Proportion of isolates required to retain length outliers

    # Gene refinding
    String refind_mode = "default"                      # Options: default, strict, off
    Int search_radius = 5000                            # Search radius for gene refinding in bp
    Float refind_prop_match = 0.2                       # Minimum proportion of gene that must match to be refound

    # Gene family behavior
    Boolean remove_invalid_genes = true                 # Remove genes that fail basic validity checks
    Boolean merge_paralogs = false                      # Merge paralogous gene families
    Boolean remove_by_consensus = true                  # Use consensus-based removal during graph cleaning
    Boolean all_seq_in_graph = false                    # Retain all sequences in graph even if filtered [set to true if debugging]

    # Alignment tool — only used when align = true
    String aligner = "mafft"                            # Options: prank, clustal, mafft, none

    # High variation flagging
    Int high_var_flag = 100                             # Flag gene families with more than this many graph cycles
                                                        # High-cycle families often indicate annotation errors or highly variable loci

    # Core entropy filter (optional, applied to core alignment)
    String? core_entropy_filter                         # Filters high-entropy positions from core alignment if set

    # Runtime and infrastructure
    # Plenty of resources for ~1,500 genomes
    Int memory = 128                                  
    Int cpu = 32        
    Int disk_size = 750
    String docker_image = "staphb/panaroo:1.6.0"
  }

  command <<<
    # Fail loudly and immediately on any error, unset variable, or failed pipe
    # -e: exit on error
    # -u: error on unset variables
    # -x: print each command to logs before executing (invaluable for debugging long runs)
    # -o pipefail: if any command in a pipe fails, the whole pipe fails
    set -euxo pipefail

    # Date/version capture
    date | tee DATE
    panaroo --version 2>&1 | tee PANAROO_VERSION_RAW

    # Strip the "panaroo " prefix so VERSION contains just the number e.g. "1.6.0"
    sed 's/panaroo //' PANAROO_VERSION_RAW | tee VERSION

    # Panaroo requires a text file listing absolute GFF paths, one per line
    # We symlink all input GFFs into a local staging directory to normalize paths,
    # then generate the file list using absolute paths for robustness
    # Using a bash array + quoted loop to safely handle any filenames with spaces
    mkdir INPUT_DIR
    gff_files=(~{sep=' ' gff3})
    for gff in "${gff_files[@]}"; do
      ln -s "$gff" INPUT_DIR/
    done
    ls -1 "$(pwd)"/INPUT_DIR/* > local_gffs.txt

    # Create Panaroo output directory
    mkdir panaroo_output

    # Build alignment flags conditionally using a bash array.
    # When align = false, the array stays empty and expands to nothing.
    ALIGN_FLAGS_ARR=()
    if [[ "~{align}" == "true" ]]; then
      ALIGN_FLAGS_ARR+=(--alignment "~{alignment_type}")
      ALIGN_FLAGS_ARR+=(--aligner "~{aligner}")
    fi

    # Run Panaroo
    panaroo \
      -i local_gffs.txt \
      -o panaroo_output/ \
      --clean-mode ~{clean_mode} \
      -c ~{seq_id} \
      -f ~{family_threshold} \
      --len_dif_percent ~{len_dif_percent} \
      --family_len_dif_percent ~{family_len_dif_percent} \
      --core_threshold ~{core_threshold} \
      ~{if defined(core_subset) then "--core_subset " + select_first([core_subset]) else ""} \
      ~{true="--remove-invalid-genes" false="" remove_invalid_genes} \
      ~{true="" false="--no_clean_edges" clean_edges} \
      --edge_support_threshold ~{edge_support_threshold} \
      --min_edge_support_sv ~{min_edge_support_sv} \
      --min_trailing_support ~{min_trailing_support} \
      --trailing_recursive ~{trailing_recursive} \
      --length_outlier_support_proportion ~{length_outlier_support_proportion} \
      --refind-mode ~{refind_mode} \
      --search_radius ~{search_radius} \
      --refind_prop_match ~{refind_prop_match} \
      ~{true="--merge_paralogs" false="" merge_paralogs} \
      ~{true="--remove_by_consensus True" false="--remove_by_consensus False" remove_by_consensus} \
      ~{true="--all_seq_in_graph" false="" all_seq_in_graph} \
      --high_var_flag ~{high_var_flag} \
      ~{if defined(core_entropy_filter) then "--core_entropy_filter " + select_first([core_entropy_filter]) else ""} \
      "${ALIGN_FLAGS_ARR[@]}" \
      -t ~{cpu}

    # Rename all outputs with cluster_name prefix for traceability (mirrors PIRATE convention)
    mv panaroo_output/summary_statistics.txt        panaroo_output/~{cluster_name}_summary_statistics.txt
    mv panaroo_output/final_graph.gml               panaroo_output/~{cluster_name}_final_graph.gml
    mv panaroo_output/pre_filt_graph.gml            panaroo_output/~{cluster_name}_pre_filt_graph.gml
    mv panaroo_output/pan_genome_reference.fa       panaroo_output/~{cluster_name}_pan_genome_reference.fa
    mv panaroo_output/gene_data.csv                 panaroo_output/~{cluster_name}_gene_data.csv
    mv panaroo_output/gene_presence_absence.Rtab    panaroo_output/~{cluster_name}_gene_presence_absence.Rtab
    mv panaroo_output/gene_presence_absence.csv     panaroo_output/~{cluster_name}_gene_presence_absence.csv
    mv panaroo_output/gene_presence_absence_roary.csv \
                                                    panaroo_output/~{cluster_name}_gene_presence_absence_roary.csv
    mv panaroo_output/combined_DNA_CDS.fasta        panaroo_output/~{cluster_name}_combined_DNA_CDS.fasta
    mv panaroo_output/combined_protein_CDS.fasta    panaroo_output/~{cluster_name}_combined_protein_CDS.fasta
    mv panaroo_output/combined_protein_cdhit_out.txt \
                                                    panaroo_output/~{cluster_name}_combined_protein_cdhit_out.txt
    mv panaroo_output/combined_protein_cdhit_out.txt.clstr \
                                                    panaroo_output/~{cluster_name}_combined_protein_cdhit_out.txt.clstr
    mv panaroo_output/struct_presence_absence.Rtab  panaroo_output/~{cluster_name}_struct_presence_absence.Rtab

    # Rename alignment outputs only if alignment was requested
    # alignment_type "pan"  produces both pan_genome_alignment.aln AND core_gene_alignment.aln; "core" only produces core_gene_alignment.aln
    if [[ "~{align}" == "true" ]]; then
      if [[ -f panaroo_output/core_gene_alignment.aln ]]; then
        mv panaroo_output/core_gene_alignment.aln panaroo_output/~{cluster_name}_core_gene_alignment.aln
      fi
      if [[ -f panaroo_output/pan_genome_alignment.aln ]]; then
        mv panaroo_output/pan_genome_alignment.aln panaroo_output/~{cluster_name}_pan_genome_alignment.aln
      fi
    fi
  >>>

  output {

    # Version and date traceability
    String date                 = read_string("DATE")
    String panaroo_version      = read_string("VERSION")
    String panaroo_docker_image = docker_image

    # Pangenome summary and graph outputs
    File panaroo_summary_stats        = "panaroo_output/" + cluster_name + "_summary_statistics.txt"
    File panaroo_final_graph          = "panaroo_output/" + cluster_name + "_final_graph.gml"
    File panaroo_pre_filt_graph       = "panaroo_output/" + cluster_name + "_pre_filt_graph.gml"
    File panaroo_pan_genome_reference = "panaroo_output/" + cluster_name + "_pan_genome_reference.fa"

    # Gene presence/absence outputs
    File panaroo_gene_data                   = "panaroo_output/" + cluster_name + "_gene_data.csv"
    File panaroo_gene_presence_absence_rtab  = "panaroo_output/" + cluster_name + "_gene_presence_absence.Rtab"
    File panaroo_gene_presence_absence_csv   = "panaroo_output/" + cluster_name + "_gene_presence_absence.csv"
    File panaroo_gene_presence_absence_roary = "panaroo_output/" + cluster_name + "_gene_presence_absence_roary.csv"

    # Structural variation presence/absence
    File panaroo_struct_presence_absence = "panaroo_output/" + cluster_name + "_struct_presence_absence.Rtab"

    # Combined sequence outputs
    File panaroo_combined_DNA_CDS                   = "panaroo_output/" + cluster_name + "_combined_DNA_CDS.fasta"
    File panaroo_combined_protein_CDS               = "panaroo_output/" + cluster_name + "_combined_protein_CDS.fasta"
    File panaroo_combined_protein_cdhit_out         = "panaroo_output/" + cluster_name + "_combined_protein_cdhit_out.txt"
    File panaroo_combined_protein_cdhit_out_cluster = "panaroo_output/" + cluster_name + "_combined_protein_cdhit_out.txt.clstr"

    # Alignment outputs — only produced when align = true
    File? panaroo_core_alignment_fasta = "panaroo_output/" + cluster_name + "_core_gene_alignment.aln"
    File? panaroo_pan_alignment_fasta  = "panaroo_output/" + cluster_name + "_pan_genome_alignment.aln"
  }

  runtime {

    docker:      docker_image
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB" # TES
    preemptible: 0
    maxRetries:  3
  }
}
