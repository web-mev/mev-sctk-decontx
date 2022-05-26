workflow SctkDecontx {
    
    # An integer matrix of counts
    File raw_counts

    call runDecontx {
        input:
            raw_counts = raw_counts
    }

    output {
        File decontaminate_output = runDecontx.fout
    }
}

task runDecontx {

    File raw_counts
    
    String output_name_prefix = "sctk_doublet_cells_reduced_counts"
    Int disk_size = 20

    command {
        Rscript /opt/software/decontX.R \
        -f ${raw_counts} \
        -o ${output_name_prefix}
    }

    output {
        File fout = glob("${output_name_prefix}*")[0]
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-sctk-decontx"
        cpu: 2
        memory: "16 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}
