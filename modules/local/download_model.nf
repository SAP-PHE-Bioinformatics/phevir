process DOWNLOAD_MODEL {

    container "quay.io/biocontainers/wget:1.20.1"

    input:
    val(model_file)

    output:
    path("medaka_model.tar.gz"), emit: model

    script:
    """
    wget -O medaka_model.tar.gz $model_file
    """
}