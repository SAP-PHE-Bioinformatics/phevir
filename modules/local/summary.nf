process SUMMARY {
  tag        "Creating summary files"
  label      "process_single"
  publishDir path: params.outdir, mode: 'copy'
  container  'quay.io/biocontainers/pandas:1.5.2'

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple file(files), val(versions), file(multiqc)

  output:
  path "mpx_results.{csv,txt}", emit: summary_file

  shell:
  def multiqc_files = multiqc.join(" ")
  """
    echo "${versions}" | cut -f 1,3,5,7,9,11  -d ',' | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/, /,/g' >  versions.csv
    echo "${versions}" | cut -f 2,4,6,8,10,12 -d ',' | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/, /,/g' | awk '{(\$1=\$1); print \$0}' >> versions.csv

    echo "Summary files are ${files}"

    mkdir multiqc_data
    for file in ${multiqc_files}
    do
      if [ -f "\$file" ]; then mv \$file multiqc_data/. ; fi
    done

    if [ -n "\$(find . -iname *ampliconstats.txt | head -n 1)" ] 
    then
      cat *_ampliconstats.txt | grep -h ^FREADS > ampliconstats.summary
    else
      touch ampliconstats.summary
    fi

    if [ -s "vadr.vadr.sqa" ] ; then tail -n +2 "vadr.vadr.sqa" | grep -v "#-" | tr -s '[:blank:]' ',' > vadr.csv ; fi

    python combine.py ${params.minimum_depth}
  """
}