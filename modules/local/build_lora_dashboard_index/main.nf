process BUILD_LORA_DASHBOARD_INDEX {
    label "process_medium"
    label "error_retry"

    input:
        val dashboards

    output:
        path("index.html"), emit: dashboard_index

    when:
        task.ext.when == null || task.ext.when

    script:
        def links = dashboards.collect { genome_id ->
            "<li><a href=\"${genome_id}/index.html\">${genome_id}</a></li>"
        }.join("\n")
        """
        printf '%s\\n' \\
            '<!doctype html>' \\
            '<html lang="en">' \\
            '<head>' \\
            '<meta charset="utf-8">' \\
            '<meta name="viewport" content="width=device-width, initial-scale=1">' \\
            '<title>LoRA Dashboard Index</title>' \\
            '<style>' \\
            'body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; color: #172033; background: #f6f8fb; }' \\
            'header { background: #172033; color: white; padding: 24px 30px; }' \\
            'main { max-width: 900px; margin: 0 auto; padding: 26px 30px; }' \\
            'a { color: #2563eb; font-weight: 700; }' \\
            'li { margin: 10px 0; }' \\
            '.card { background: white; border: 1px solid #d8e1ea; border-radius: 8px; padding: 18px; }' \\
            '</style>' \\
            '</head>' \\
            '<body>' \\
            '<header><h1>LoRA Dashboard</h1><div>${params.runID}</div></header>' \\
            '<main><div class="card"><h2>Reference Genomes</h2><ul>' \\
            '${links}' \\
            '</ul></div></main>' \\
            '</body>' \\
            '</html>' \\
            > index.html
        """
}
