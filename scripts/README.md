After creating distance files between pdbId-ligandId pairs,

'mindist' below can be substituted with any of 9 possible distance metrics.
    - mindist
    - meandist
    - fracin4
    - maxstd
    - meanstd
    - sumstd
    - maxvdw
    - meanvdw
    - sumvdw

0. Use hmmr to identify domains

    Resulting file at /domains/BioLiP_2018-09-12-domains-pfam_v31.tsv.gz

    tab-delimited file with columns PDB ID-PDB Chain, Domain Name (unique), and
    comma-delimited list of (1-indexed domain match state : 0-indexed sequence position : amino acid value)

1. Use the file containing all domain hits to create per-domain, per-ligand alignments

    eval_uniqueness --create_alignments --distance mindist

    a) creates alignment files with names:
    /domains/alignments/mindist/<domain>_<super_ligand_type>_<distance_metric>.tmpfa
    e.g.
    /domains/alignments/mindist/PF00036_EF-hand_1_ALL__mindist.tmpfa

    i.e. each file has to do with one domain and 1 ligand 'type'
    each line of this file is of the form:
        2lqcA_10_38  1:10:E,2:11:F, ..
        (i.e. pdb id_start_end -> hmmer scores)

    b) convert each file created in (a) above to file with name:
    /domains/alignments/mindist/<domain>_<super_ligand_type>_<distance_metric>.aln.fa

2. eval_uniqueness --distance mindist

    creates file:
    /domains/uniqueness_scores_mindist.txt.gz

3. generate_domain_scores

    creates files:
    /domains/binding_scores/mindist/PF0009_GTD_EFTU_binding_scores.mindist.txt.gz
    ...

4. cross_validate_scores

    creates files:
    /domains/cross_validation/mindist/precision_recall
        PF0009_GTD_EFTU_pr.mindist.txt.gz
        ...
    /domains/cross_validation/mindist/precision_threshold
    /domains/cross_validation/mindist/receiver_operator

5. interacdome_webserver --webserver

    creates files:
    /interacdome-webserver/pfms/*.pfm
    /interacdome-webserver/interacdome_allresults.tsv
    /nteracdome-webserver/interacdome_fordownload.tsv


