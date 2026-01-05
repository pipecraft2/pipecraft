process create_final_table {
    input:
        val(json_files)  // This is a list of file paths (strings or Path objects)
        val(run_id)
        val(use_itsx)
        val(rel_abu_threshold)

    output:
        path("${run_id}.results.xlsx"), emit: results_excel

    publishDir "${run_id}_results/", mode: 'move'

    script:
        // Join all file paths into a single string
        def files_arg = json_files.collect { it.toString() }.join(" ")
        """
        #!/usr/bin/env python3
        import sys
        import pandas as pd
        import json

        # Parse arguments from Nextflow variables
        json_file_paths = "${files_arg}".split()
        rel_abu_threshold = float("${rel_abu_threshold}")
        output_file = "${run_id}.results.xlsx"
        itsx_used = int("${use_itsx}")

        empty_json_dfs = []
        json_dfs = []
        for json_path in json_file_paths:
            with open(json_path, 'r') as f:
                json_data = json.load(f)
                print(json_path)
                # manage nonprocessed barcodes
                if "cluster_data" not in json_data:
                    empty_data = [{
                        "Sample": json_data["barcode_id"],
                        "Number of clusters": json_data["message"]
                    }]
                    empty_df = pd.DataFrame(empty_data)
                    empty_json_dfs.append(empty_df)
                    continue
                
                clusters_data = []
                for cluster in json_data["cluster_data"]:
                    clusters_data.append({
                        "Cluster ID": cluster["cluster_id"],
                        "Cluster size": cluster["cluster_size"],
                        "Cluster relative abundance": cluster["relative_abundance"] * 100,
                        "Cluster sequence": cluster["cluster_sequence"],
                        "Cluster sequence untrimmed": cluster.get("cluster_sequence_untrimmed", ""),
                        "BLASTn taxonomy assignment": cluster.get("blastn_tax_name", "No BLAST hit"),
                        "BLASTn perc. ident.": cluster.get("blastn_pident", "N/A"),
                        "BLASTn query coverage": cluster.get("blastn_query_coverage", "N/A"),
                        "BLASTn query length": cluster.get("blastn_query_length", "N/A"),
                        "BLASTn subject length": cluster.get("blastn_subject_length", "N/A"),
                        "BLASTn evalue": cluster.get("blastn_evalue", "N/A"),
                        "BLASTn subject SH": cluster.get("blastn_sh_id", "N/A"),
                        "BLASTn full taxonomy": cluster.get("blastn_full_taxonomy", "N/A")
                    })
                    
                clusters_df = pd.DataFrame(clusters_data)
                clusters_df["Sample"] = json_data["barcode_id"]
                clusters_df["Number of clusters"] = json_data["number_of_clusters"]
                clusters_df["Total reads after filtering"] = json_data["total_reads_after_filtering"]
                # Reorder columns
                cols = ["Sample", "Number of clusters", "Total reads after filtering"] + \\
                       [col for col in clusters_df.columns if col not in ["Sample", "Number of clusters", "Total reads after filtering"]]
                clusters_df = clusters_df[cols]
                json_dfs.append(clusters_df)

        if len(json_dfs) == 0:
            df = pd.DataFrame()
        else:
            df = pd.concat(json_dfs, ignore_index=True)

            # apply rel abu threshold
            df = df[df["Cluster relative abundance"] >= rel_abu_threshold]
            
            # Count passed clusters per sample
            passed_clusters = df.groupby("Sample").size().reset_index(name="Number of passed clusters")
            df = df.merge(passed_clusters, on="Sample", how="left")
            
            # Reorder to place "Number of passed clusters" after "Number of clusters"
            cols = df.columns.tolist()
            if "Number of passed clusters" in cols and "Number of clusters" in cols:
                cols.remove("Number of passed clusters")
                idx = cols.index("Number of clusters") + 1
                cols.insert(idx, "Number of passed clusters")
                df = df[cols]
            
            df = df.sort_values(by=["Sample", "Cluster relative abundance"], ascending=[True, False])

        if empty_json_dfs:
            # put nonprocessed barcodes on top
            to_concat = empty_json_dfs + [df]
            df = pd.concat(to_concat, ignore_index=True)

        # mark no itx extraction
        if itsx_used != 1:
            df["Cluster sequence"] = "No ITSx sequence extraction. Taxonomy was assigned based on full, untrimmed sequence."

        df.to_excel(output_file, engine="openpyxl", index=False)
        """
}