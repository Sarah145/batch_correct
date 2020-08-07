#!/usr/bin/env nextflow




Channel.fromPath(params.dataset_list)
  .flatMap{ it.readLines() }
  .view()
  .set { DATALIST }

// fetch data
process get_datasets {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
 	tag "$datasetname"

    	input:
    	val datasetname from DATALIST
    	output:
    	set val(datasetname), file('*.Rds') into RUN_SEU
    	
    
		
    	shell:
    	'''
    	Rfile="!{params.data_dir}/!{datasetname}.Rds"
    	
    	if [[ ! -e $Rfile ]]; then
    	  echo "Please check existence of $Rfile"
    	  false
    	fi
    	ln -s $Rfile .
    	'''
}

process Seurat_3{
 	publishDir "./data/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Seurat 3 $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from RUN_SEU
    	
	output:
    	set val(datasetname), val('Seurat3'), file('Seurat3.*.rds')	
	
    	"""
    	seurat_integrate.R\
		--input_object ${datain}\
		--batch_key ${params.batch_key}\
		--output_object Seurat3.${datasetname}.rds
    	"""
}
