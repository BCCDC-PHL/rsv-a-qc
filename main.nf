#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { nextclade_dataset } from './modules/rsv-a-qc.nf'
include { identify_complete_genomes } from './modules/rsv-a-qc.nf'
include { prepare_multi_fasta } from './modules/rsv-a-qc.nf'
include { nextclade } from './modules/rsv-a-qc.nf'
include { detect_resistance_mutations } from './modules/rsv-a-qc.nf'
include { augur_align } from './modules/rsv-a-qc.nf'
include { augur_tree } from './modules/rsv-a-qc.nf'
include { make_alleles } from './modules/rsv-a-qc.nf'
include { plot_tree_snps } from './modules/rsv-a-qc.nf'
include { build_snpeff_db } from './modules/rsv-a-qc.nf'
include { snpeff } from './modules/rsv-a-qc.nf'
include { make_aa_table } from './modules/rsv-a-qc.nf'
include { primer_bed_to_amplicon_bed } from './modules/rsv-a-qc.nf'
include { calc_amplicon_depth } from './modules/rsv-a-qc.nf'
include { calc_per_base_depth } from './modules/rsv-a-qc.nf'
include { create_primer_snp_bed } from './modules/rsv-a-qc.nf'
include { make_genome_bed } from './modules/rsv-a-qc.nf'
include { make_sample_qc_summary } from './modules/rsv-a-qc.nf'
include { write_qc_summary } from './modules/rsv-a-qc.nf'


workflow {

  ch_run_dir            = Channel.fromPath(params.run_dir)

  ch_run_id             = Channel.of(file(params.run_dir).getName())

  ch_artic_analysis_dir = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir)

  ch_consensus = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir + "/" + params.artic_consensus_subdir + "/*${params.artic_consensus_filename_suffix}").map{ it -> [it.baseName.split("\\.")[0], it] }

  ch_variants = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir + "/" + params.artic_variants_subdir + "/*${params.artic_variants_filename_suffix}").map{ it -> [it.baseName.split("\\.")[0], it] }

  ch_alignments = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir + "/" + params.artic_alignment_subdir + "/*${params.artic_alignment_filename_suffix}").map{ it -> [it.baseName.split("\\.")[0], it] }

  ch_alignment_indexes = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir + "/" + params.artic_alignment_subdir + "/*${params.artic_alignment_filename_suffix}.bai").map{ it -> [it.baseName.split("\\.")[0], it] }

  ch_alignments_with_index = ch_alignments.join(ch_alignment_indexes)

  ch_primer_bed = Channel.fromPath(params.bed)

  ch_primer_pairs_tsv = Channel.fromPath(params.primer_pairs_tsv)

  ch_resistance_mutations = Channel.fromPath(params.resistance_mutations)

  ch_snpEff_config = Channel.fromPath(params.snpEff_config)

  ch_snpEff_dataDir = Channel.fromPath(params.snpEff_dataDir)


  main:

    identify_complete_genomes(ch_run_id.combine(ch_artic_analysis_dir))

    prepare_multi_fasta(identify_complete_genomes.out.combine(ch_artic_analysis_dir))

    nextclade_dataset(ch_run_id)

    nextclade(prepare_multi_fasta.out.join(nextclade_dataset.out.dataset))

    detect_resistance_mutations(nextclade.out.fcds.combine(ch_resistance_mutations))

    augur_align(prepare_multi_fasta.out.join(nextclade_dataset.out.ref))

    augur_tree(augur_align.out.join(nextclade_dataset.out.ref))

    make_alleles(augur_align.out.join(nextclade_dataset.out.ref))

    plot_tree_snps(augur_tree.out.join(make_alleles.out).join(nextclade.out.qc))

    if (params.build_snpeff_db) {
      build_snpeff_db(nextclade_dataset.out.ref)
    }

    snpeff(ch_variants.combine(nextclade_dataset.out.ref.map{ it -> it[1] }).combine(ch_snpEff_config).combine(ch_snpEff_dataDir))

    make_aa_table(snpeff.out)

    primer_bed_to_amplicon_bed(ch_primer_bed.combine(ch_primer_pairs_tsv))

    calc_amplicon_depth(ch_alignments_with_index.combine(primer_bed_to_amplicon_bed.out))

    create_primer_snp_bed(ch_variants.combine(ch_primer_bed))

    make_genome_bed(nextclade_dataset.out.ref)

    calc_per_base_depth(ch_alignments_with_index.combine(make_genome_bed.out))

    make_sample_qc_summary(ch_consensus.join(ch_variants).join(calc_per_base_depth.out).join(make_aa_table.out).combine(make_alleles.out))

    write_qc_summary(ch_run_id.combine(make_sample_qc_summary.out.map{ it -> it[1] }.splitCsv(sep: '\t').unique().toList()).map{ it -> [it[0], it[1..-1]] })
}
