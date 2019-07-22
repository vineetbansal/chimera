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

1. eval_uniqueness --create-alignments --distance mindist

    creates files:
    /domains/alignments/mindist/*.tmp.fa

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


