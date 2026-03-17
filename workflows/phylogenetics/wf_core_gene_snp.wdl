version 1.0

import "../../tasks/phylogenetic_inference/task_iqtree.wdl" as iqtree
import "../../tasks/phylogenetic_inference/task_panaroo.wdl" as panaroo_task
import "../../tasks/phylogenetic_inference/utilities/task_reorder_matrix.wdl" as reorder_matrix
import "../../tasks/phylogenetic_inference/utilities/task_snp_dists.wdl" as snp_dists
import "../../tasks/phylogenetic_inference/utilities/task_snp_sites.wdl" as snp_sites
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_summarize_data.wdl" as data_summary


workflow core_gene_snp_workflow {
  input {
    Array[File] gff3
    String cluster_name

    # If align = true, Panaroo will produce core and/or pangenome alignments for the sample set.
    # align = false will produce pangenome summary and presence/absence outputs only,
    # with no downstream SNP or phylogenetic analysis.
    Boolean align = true

    # alignment_type controls which alignment Panaroo produces:
    #   "core" — core genome alignment only (feeds core SNP pipeline)
    #   "pan"  — pangenome alignment + core genome alignment
    # Note: if pan_tree = true, alignment_type should be set to "pan" so the
    # pan alignment FASTA is available for the pan_tree branch below.
    String alignment_type = "core"

    # use core_tree = true to produce a phylogenetic tree and SNP distance matrix
    # from the core genome alignment
    Boolean core_tree = true

    # use pan_tree = true to produce a phylogenetic tree and SNP distance matrix
    # from the pangenome alignment. Requires alignment_type = "pan".
    Boolean pan_tree = false

    # -------------------------------------------------------------------------
    # Panaroo pangenome construction parameters
    # See task_panaroo.wdl for full parameter documentation
    # -------------------------------------------------------------------------

    # Core genome definition
    Float panaroo_core_threshold = 0.98   # Proportion of isolates a gene must be in to be core
    String? panaroo_core_subset           # Optional: define core relative to a subset of isolates

    # Pangenome construction behavior
    String panaroo_clean_mode = "sensitive"     # Options: strict, moderate, sensitive
    Float panaroo_seq_id = 0.95                 # Sequence identity threshold for initial clustering
    Float panaroo_family_threshold = 0.7        # Family-level sequence identity threshold
    Float panaroo_len_dif_percent = 0.98        # Maximum length difference to cluster sequences
    Float panaroo_family_len_dif_percent = 0.98 # Family-level length difference threshold

    # Graph and edge parameters
    Boolean panaroo_clean_edges = true                      # Clean edges in the pangenome graph
    Float panaroo_edge_support_threshold = 0.0              # Minimum edge support threshold
    Int panaroo_min_edge_support_sv = 5                     # Minimum isolates sharing edge for SV calls
    Int panaroo_min_trailing_support = 2                    # Minimum cluster size at contig ends
    Int panaroo_trailing_recursive = 2                      # Number of recursive trailing removals
    Float panaroo_length_outlier_support_proportion = 0.01  # Proportion to retain length outliers

    # Gene refinding
    String panaroo_refind_mode = "default"  # Options: default, strict, off
    Int panaroo_search_radius = 5000        # Search radius for gene refinding in bp
    Float panaroo_refind_prop_match = 0.2   # Minimum proportion of gene to be refound

    # Gene family behavior
    Boolean panaroo_remove_invalid_genes = true   # Remove genes failing basic validity checks
    Boolean panaroo_merge_paralogs = false         # Merge paralogous gene families
    Boolean panaroo_remove_by_consensus = true     # Use consensus-based removal during graph cleaning
    Boolean panaroo_all_seq_in_graph = false       # Retain all sequences in graph even if filtered

    # Alignment tool
    String panaroo_aligner = "mafft"  # Options: prank, clustal, mafft, none

    # High variation flagging
    Int panaroo_high_var_flag = 100   # Flag gene families with more graph cycles than this threshold

    # Core entropy filter — optional
    String? panaroo_core_entropy_filter

    # Panaroo runtime
    Int panaroo_memory = 128
    Int panaroo_cpu = 32
    Int panaroo_disk_size = 750
    String panaroo_docker = "staphb/panaroo:1.6.0"

    # -------------------------------------------------------------------------
    # Data summary optional inputs
    # -------------------------------------------------------------------------
    Array[String]? sample_names
    String? data_summary_terra_project
    String? data_summary_terra_workspace
    String? data_summary_terra_table
    String? data_summary_column_names
    Boolean phandango_coloring = false

    # Tree options
    Boolean midpoint_root_tree = true
  }

  # Replace spaces in cluster_name with underscores for safe use in filenames
  String cluster_name_updated = sub(cluster_name, " ", "_")

  # -------------------------------------------------------------------------
  # Pangenome construction with Panaroo
  # -------------------------------------------------------------------------
  call panaroo_task.panaroo {
    input:
      gff3                              = gff3,
      cluster_name                      = cluster_name_updated,
      align                             = align,
      alignment_type                    = alignment_type,
      core_threshold                    = panaroo_core_threshold,
      core_subset                       = panaroo_core_subset,
      clean_mode                        = panaroo_clean_mode,
      seq_id                            = panaroo_seq_id,
      family_threshold                  = panaroo_family_threshold,
      len_dif_percent                   = panaroo_len_dif_percent,
      family_len_dif_percent            = panaroo_family_len_dif_percent,
      clean_edges                       = panaroo_clean_edges,
      edge_support_threshold            = panaroo_edge_support_threshold,
      min_edge_support_sv               = panaroo_min_edge_support_sv,
      min_trailing_support              = panaroo_min_trailing_support,
      trailing_recursive                = panaroo_trailing_recursive,
      length_outlier_support_proportion = panaroo_length_outlier_support_proportion,
      refind_mode                       = panaroo_refind_mode,
      search_radius                     = panaroo_search_radius,
      refind_prop_match                 = panaroo_refind_prop_match,
      remove_invalid_genes              = panaroo_remove_invalid_genes,
      merge_paralogs                    = panaroo_merge_paralogs,
      remove_by_consensus               = panaroo_remove_by_consensus,
      all_seq_in_graph                  = panaroo_all_seq_in_graph,
      aligner                           = panaroo_aligner,
      high_var_flag                     = panaroo_high_var_flag,
      core_entropy_filter               = panaroo_core_entropy_filter,
      memory                            = panaroo_memory,
      cpu                               = panaroo_cpu,
      disk_size                         = panaroo_disk_size,
      docker_image                      = panaroo_docker
  }

  # -------------------------------------------------------------------------
  # Core genome SNP analysis
  # Runs when align = true AND core_tree = true
  # -------------------------------------------------------------------------
  if (align) {
    if (core_tree) {
      # Extract SNP sites from core genome alignment
      call snp_sites.snp_sites as core_snp_sites {
        input:
          msa_fasta             = select_first([panaroo.panaroo_core_alignment_fasta]),
          output_name           = cluster_name_updated + "_core",
          allow_wildcard_bases  = true,
          output_vcf            = false,
          output_phylip         = false,
          output_multifasta     = true,
          output_pseudo_ref     = false,
          output_monomorphic    = false
      }
      # Build maximum likelihood phylogenetic tree from core SNP alignment
      call iqtree.iqtree as core_iqtree {
        input:
          alignment    = select_first([core_snp_sites.snp_sites_multifasta]),
          cluster_name = cluster_name_updated
      }
      # Compute pairwise SNP distance matrix from core SNP alignment
      call snp_dists.snp_dists as core_snp_dists {
        input:
          alignment    = select_first([core_snp_sites.snp_sites_multifasta]),
          cluster_name = cluster_name_updated
      }
      # Reorder SNP matrix to match tree leaf order for visualization
      call reorder_matrix.reorder_matrix as core_reorder_matrix {
        input:
          input_tree        = core_iqtree.ml_tree,
          matrix            = core_snp_dists.snp_matrix,
          cluster_name      = cluster_name_updated + "_core",
          midpoint_root_tree = midpoint_root_tree,
          phandango_coloring = phandango_coloring
      }
    }

    # -------------------------------------------------------------------------
    # Pangenome SNP analysis
    # Runs when align = true AND pan_tree = true
    # Requires alignment_type = "pan" — pan_genome_alignment.aln must be present
    # -------------------------------------------------------------------------
    if (pan_tree) {
      # Build maximum likelihood phylogenetic tree from pangenome alignment
      call iqtree.iqtree as pan_iqtree {
        input:
          alignment    = select_first([panaroo.panaroo_pan_alignment_fasta]),
          cluster_name = cluster_name_updated
      }
      # Compute pairwise SNP distance matrix from pangenome alignment
      call snp_dists.snp_dists as pan_snp_dists {
        input:
          alignment    = select_first([panaroo.panaroo_pan_alignment_fasta]),
          cluster_name = cluster_name_updated
      }
      # Reorder SNP matrix to match tree leaf order for visualization
      call reorder_matrix.reorder_matrix as pan_reorder_matrix {
        input:
          input_tree         = pan_iqtree.ml_tree,
          matrix             = pan_snp_dists.snp_matrix,
          cluster_name       = cluster_name_updated + "_pan",
          midpoint_root_tree = midpoint_root_tree,
          phandango_coloring = phandango_coloring
      }
    }
  }

  # -------------------------------------------------------------------------
  # Optional data summary task
  # Runs only when data_summary_column_names is provided
  # -------------------------------------------------------------------------
  if (defined(data_summary_column_names)) {
    call data_summary.summarize_data {
      input:
        sample_names    = sample_names,
        terra_project   = data_summary_terra_project,
        terra_workspace = data_summary_terra_workspace,
        terra_table     = data_summary_terra_table,
        column_names    = data_summary_column_names,
        output_prefix   = cluster_name_updated,
        phandango_coloring = phandango_coloring
    }
  }

  # -------------------------------------------------------------------------
  # Workflow version capture
  # -------------------------------------------------------------------------
  call versioning.version_capture {
    input:
  }

  # -------------------------------------------------------------------------
  # Workflow outputs
  # -------------------------------------------------------------------------
  output {
    # Workflow version
    String core_gene_snp_wf_version       = version_capture.phb_version
    String core_gene_snp_wf_analysis_date = version_capture.date

    # -----------------------------------------------------------------------
    # Panaroo outputs
    # -----------------------------------------------------------------------

    # Version traceability
    String panaroo_version      = panaroo.panaroo_version
    String panaroo_docker_image = panaroo.panaroo_docker_image

    # Pangenome summary and graph outputs
    File panaroo_summary_stats        = panaroo.panaroo_summary_stats
    File panaroo_final_graph          = panaroo.panaroo_final_graph
    File panaroo_pre_filt_graph       = panaroo.panaroo_pre_filt_graph
    File panaroo_pan_genome_reference = panaroo.panaroo_pan_genome_reference

    # Gene presence/absence outputs
    File panaroo_gene_data                   = panaroo.panaroo_gene_data
    File panaroo_gene_presence_absence_rtab  = panaroo.panaroo_gene_presence_absence_rtab
    File panaroo_gene_presence_absence_csv   = panaroo.panaroo_gene_presence_absence_csv
    File panaroo_gene_presence_absence_roary = panaroo.panaroo_gene_presence_absence_roary

    # Structural variation
    File panaroo_struct_presence_absence = panaroo.panaroo_struct_presence_absence

    # Combined sequence outputs
    File panaroo_combined_DNA_CDS                   = panaroo.panaroo_combined_DNA_CDS
    File panaroo_combined_protein_CDS               = panaroo.panaroo_combined_protein_CDS
    File panaroo_combined_protein_cdhit_out         = panaroo.panaroo_combined_protein_cdhit_out
    File panaroo_combined_protein_cdhit_out_cluster = panaroo.panaroo_combined_protein_cdhit_out_cluster

    # Alignment outputs — present only when align = true
    File? panaroo_core_alignment_fasta = panaroo.panaroo_core_alignment_fasta
    File? panaroo_pan_alignment_fasta  = panaroo.panaroo_pan_alignment_fasta

    # -----------------------------------------------------------------------
    # snp_sites outputs — present only when align = true and core_tree = true
    # -----------------------------------------------------------------------
    File?   snp_sites_multifasta = core_snp_sites.snp_sites_multifasta
    String? snp_sites_version    = core_snp_sites.snp_sites_version
    String? snp_sites_docker     = core_snp_sites.snp_sites_docker

    # -----------------------------------------------------------------------
    # snp_dists outputs
    # -----------------------------------------------------------------------
    String? snp_dists_version = select_first([core_snp_dists.snp_dists_version, pan_snp_dists.snp_dists_version, ""])

    # -----------------------------------------------------------------------
    # IQTree outputs
    # -----------------------------------------------------------------------
    String? iqtree_version = select_first([core_iqtree.version, pan_iqtree.version, ""])

    # -----------------------------------------------------------------------
    # Reorder matrix outputs
    # -----------------------------------------------------------------------
    File?   core_snp_matrix  = core_reorder_matrix.ordered_matrix
    File?   core_tree_output = core_reorder_matrix.tree
    File?   pan_snp_matrix   = pan_reorder_matrix.ordered_matrix
    File?   pan_tree_output  = pan_reorder_matrix.tree

    # -----------------------------------------------------------------------
    # Data summary outputs — present only when data_summary_column_names is set
    # -----------------------------------------------------------------------
    File? summarized_data    = summarize_data.summarized_data
    File? filtered_metadata  = summarize_data.filtered_metadata
  }
}
