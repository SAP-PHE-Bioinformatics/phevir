process PANGO_COLLAPSE {
    tag "collapsing samples"
    label 'process_medium'
    container "quay.io/uphl/pango-collapse:0.8.2-2024-09-10"


    input:
    path(nextclade)
    path collapse_file
    val runID

    output:
    path("*.tsv"), emit: collapsed
    path  "versions.yml"            , emit: versions
    script:
    def args = task.ext.args ?: ''
    def prefix = "${runID}"
    """
    pango-collapse $nextclade --output ${prefix}_collapsed.tsv --collapse-file $collapse_file --collapse-column VOC_Lineage -l Nextclade_pango

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pango-collapse: \$(pango-collapse -v 2>&1 | head -n1 | sed 's/^.*pango-collapse //; s/ .*\$//')
    END_VERSIONS
    """
}