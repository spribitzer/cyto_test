#### Compensate_fcs.r
# Description: Takes a users' query from Advanced search; 
#    downloads .fcs files; finds the appropriate compensation matrix for each .fcs file; 
#    uses flowCore::compensate() for each .fcs file; and exports to '/home/jupyter/compensate_fcs' directory
# Arguments: 
#.  query_id : ID from Advanced Search 
# Output: .fcs files prefixed with "compensated_<original_fcs_fileName>.fcs" 
# NOTE: only works with queryID! 
# Contributors: James Harvey 
####


### 
# GLOBALS 
###
FILETYPE <- 'FlowCytometry-supervised-comp'
OUTDIR <- '/home/jupyter/compensate_fcs'


###
# Libraries 
### 
library(hise) 
library(flowCore) 
library(data.table) 
library(stringr) 


library(unix)
rlimit_as(1e12)


#### 
# Description: downloads .fcs files onto IDE and returns descriptor object
# Arguments: qid: queryId from Advanced search
# output: cacheFiles() object
####
download_fcs <- function(qid) { 
    print('downloading .fcs files based on QueryID from user...') 
    cache_out <- cacheFiles(list(qid))
    return(cache_out)
}


###
# Description: creates list of ids from a cacheResult() object 
# Arguments: 
#   cache_result (list) : output from cacheFiles() SDK method 
#   field_name (str) : possible values: [id,name]
# output: list of ids 
create_id_list <- function(cache_result, field_name) { 
    if (field_name == 'id') { 
        id_list <- lapply(cache_result, function(x) { 
            x$descriptors$file$id})
    } 
    else { 
        id_list <- lapply(cache_result, function(x) { 
            x$descriptors$file$name})
    } 
    return(id_list)
} 


###
# Description: Goes through each .fcs file downloaded and downloads the appropriate compensation matrix file 
#
# Note: There's an assumption here that should be carefully considered. For any unique (panel, batchID), the same compensation file is used. 
# There exists individual files for each subject within a given (panel, batchID), we'll just take the first entry when you loop through each item 
# in this list of files from cacheFiles(). 
###
download_comp_mats <- function(cache_result) { 
    print('downloading compensation matrices...') 
    
    # a lambda function that assumes a object from cacheFiles() as input.
    # this essentially looks at the result you get from the users' queryID input;
    # loops through each entry and call getDescriptors() to pull the compensation 
    # matrices by querying on panel and batchID
    create_filter <- function(x) { 
    filter_list <- list(file.panel=x$descriptors$file$panel, 
                file.batchID=x$descriptors$file$batchID)
    desc <- getFileDescriptors(fileType=FILETYPE, 
                               filter=filter_list)[[1]][[1]]
    }
    desc_list <- lapply(cache_result, create_filter)
    
    # now we get all the fileIDs to download all the compensation matrices onto the IDE 
    compensation_fileId_list <- lapply(desc_list, function(x) { 
        x$file$id
        })
    
    # and now download
    cache_out_list <- sapply(compensation_fileId_list, function(x) { 
        cacheFiles(list(x))
    }) 
    
    return(cache_out_list)
}


###
# Description: Goes through each .fcs file downloaded and download all control files for all batches listed 
# NOTE: assumes all batch_control files has subjectGuid == 'HMN169517'
###
download_batch_control_files <- function(cache_result) { 
    # id to determine whether a file is a batch control file 
    batch_control_id <- 'HMN169517'
    unique_batches <- unique(lapply(cache_result, function(x) { 
        x$descriptors$file$batchID}))
    unique_panels <- unique(lapply(cache_result, function(x) { 
        x$descriptors$file$panel}))
    
    # pull descriptors/fileIds using query(subjectGuid='HMN169517, batchID=this_batchID) 
    q_filter <- list(subject.subjectGuid=batch_control_id, file.batchID=unique_batches, file.panel=unique_panels) 
    batch_control_desc <- getFileDescriptors(fileType='FlowCytometry', filter=q_filter)
    
    # loop through all our query results, and download the files.
    # also save the cache_result for each entry 
    cache_batch_control_list <- c()
    for (i in 1:length(batch_control_desc)) { 
        for (j in 1:length(batch_control_desc[[i]])) { 
            tmp_cache <- cacheFiles(list(batch_control_desc[[i]][[j]]$file$id))
            cache_batch_control_list <- c(cache_batch_control_list, tmp_cache) 
        }
    }
    
    
    # download those files 
    return(cache_batch_control_list) 
    }

###
# Description: For a given panel and batch, apply the compensation matrix to the bridging control and export 
###
apply_compensation_to_batch_control <- function(batch, panel) { 
    
    
    outdir <- paste0(OUTDIR, '/', 'final_batch_control')
    dir.create(outdir, showWarnings = FALSE)
    comp_matrix_dir <- paste0(OUTDIR, '/comp_matrix/')
    batch_control_dir <- paste0(OUTDIR, '/raw_batch_control/')
    
    batch_control_fileN <- list.files(batch_control_dir, pattern=paste0(batch, '_', panel))
    print(batch_control_fileN) 
    comp_matrix_fileN <- list.files(comp_matrix_dir, pattern=paste0(batch, '-', panel))

    comp_mat <- fread(paste0(comp_matrix_dir, comp_matrix_fileN))
    comp_mat[,V1:=NULL] # original export attached this V1 column 
    
    batch_control <- readRDS(paste0(batch_control_dir, batch_control_fileN))

    compensated_bridging_control <- flowCore::compensate(x=batch_control, spillover=comp_mat)
    batch_control_fileN <- str_replace(batch_control_fileN, '.fcs', '.rds') # save as RDS instead of 
    saveRDS(compensated_bridging_control, file=paste0(outdir, '/comp_', batch_control_fileN))
                        
}


###
# Description: Takes the list of ids and names to load in the .fcs file and corresponding matrix.
# then uses the compensate function and exports a new .fcs file 
# Arguments: 
#   fcs_fileId_list (list) : string of fileIds for fcs files 
#   compmat_fileId_list (list) : string of fileIds for compensation matrices
# Note:
#.   - orders of lists matter. the first entry in the compmat_fileId_list will be applied to the first entry of the  .fcs fileId list
#    - length(fcs_fileId_list) == length(compmat_fileId_list) must be true
#        (Why this has to be the case, please see note in download_comp_mats()). 
##
apply_export_compensation <- function(fcs_fileId_list, fcs_fileName_list,
                               compMat_fileId_list, compMat_fileName_list, batch_control_cache) { 
    print('applying compsenation...')
    # find where this script is saved and assume that's where the files of interest are downloaded 
    indir <- getwd() 
    
    # create output directory 
    dir.create(OUTDIR, showWarnings = FALSE) 
    
    # create batch_control dir; save batch control files 
    dir.create(paste0(OUTDIR, '/', 'raw_batch_control'), showWarnings = FALSE)
    for (i in 1:length(batch_control_cache)) { 
        # load the file 
        f_tmp <- read.FCS(paste0(indir,
                                 '/',
                                 batch_control_cache[[i]]$filePath))
        control_file_name <- stringr::word(stringr::str_replace_all(batch_control_cache[[i]]$filePath, '/', ' '), 3,3)
        control_file_name <- stringr::str_replace(control_file_name, '.fcs', '.rds') 
        saveRDS(f_tmp, paste0(OUTDIR,
                                '/',
                                'raw_batch_control/', 
                                control_file_name))
    }
    
    for (i in 1:length(fcs_fileId_list)) { 
        # read the fcs_file
        fcsName <- word(str_replace(fcs_fileName_list[i], '/', ' '), 2,2)
        fcs_f <- read.FCS(paste0(indir, 
                                 '/cache/', 
                                 as.character(fcs_fileId_list[i]), 
                                 '/', 
                                 fcsName))
        
        # read the compensation file 
        compName <- word(str_replace_all(compMat_fileName_list[i], '/', ' '), 6,6)
        comp_f <- fread(paste0(indir, 
                               '/cache/', 
                               as.character(compMat_fileId_list[i]),
                               '/',
                               compName))
        
        
        
        # utilize flowCore compensate function 
        compensated_fcs <- flowCore::compensate(x=fcs_f, spillover=comp_f)
        
        # export as .fcs suffixed with "_comp"
        print(paste0('export file number ',i, ' out of', length(fcs_fileId_list)))
        comp_fileN <- word(str_replace(fcs_fileName_list[i], '/', ' '), 2,2)
        comp_fileN <- stringr::str_replace(comp_fileN, '.fcs', '.rds')
        
        # get panel and batch; create a directory 
        thisBatch <- word(str_replace_all(fcsName, '_', ' '), 1, 1) 
        thisPanel <- word(str_replace_all(fcsName, '_', ' '), 2, 2)
        
        
        # create directories for raw_fcs; compensation_matrix; and compensated_fcs 
        rawDir <- 'raw_fcs'
        matDir <- 'comp_matrix' 
        newFCSDir <- 'final_fcs'
        for (i in c(rawDir, matDir, newFCSDir)) {     
            #thisDir <- (paste0(OUTDIR, '/', thisBatch, '_', thisPanel))
            dir.create(paste0(OUTDIR, '/', i), showWarnings = FALSE) 
        
            # export the raw data, as well as the compensated file 
            if (i == rawDir) { 
                #saveRDS(fcs_f, file = paste0(OUTDIR, '/', i, '/', fcsName))
                saveRDS(fcs_f, file = paste0(OUTDIR, '/', i, '/', stringr::str_replace(fcsName, '.fcs', '.rds'))) 
                }
            else if (i == matDir) { 
                write.csv(comp_f, paste0(OUTDIR, '/',i,'/', compName)) 
            }
            else if (i == newFCSDir) { 
                saveRDS(compensated_fcs, file=paste0(OUTDIR, '/', i, 
                                              '/', 
                                              word(str_replace(comp_fileN, '\\.', ' '), 1,1), '_comp.rds'))
            }
        }
        
        # apply compensation to briding controls and export; 
        apply_compensation_to_batch_control(thisBatch, thisPanel) 
    }

}




# wrapper function. Makes it easy to just insert this and source it in another script. 
# Run_compensate() after source('.../compensate_fcs.r') 
Run_compensate <- function(query_id, return_all_fcs) { 
    fcs_cache_result <- download_fcs(query_id) 
    
    # NOTE: dependency from download_fcs here 
    compMat_cache_result <- download_comp_mats(fcs_cache_result) 
    
    # download batch control file 
    batch_control_cache_result <- download_batch_control_files(fcs_cache_result)
    
    # get lists of fileNames and fileIds for both fcs and compensation matrices 
    fcs_file_ids <- create_id_list(fcs_cache_result, 'id')
    fcs_file_names <- create_id_list(fcs_cache_result, 'name') 
    comp_file_ids <- create_id_list(compMat_cache_result, 'id') 
    comp_file_names <- create_id_list(compMat_cache_result, 'name') 
    
    apply_export_compensation(fcs_fileId_list=fcs_file_ids, 
                              fcs_fileName_list=fcs_file_names, 
                              compMat_fileId_list=comp_file_ids, 
                              compMat_fileName_list=comp_file_names,
                              batch_control_cache_result)
    

    print('ALL FILES EXPORTED SUCCESSFULLY!!') 
    
} 
                                                                                                           

         
## look at users' input
# NOTE: if sourcing in an interactive session, run "options(run.compensate=FALSE)" first. Otherwise, the script will try running 
if (getOption('run.compensate', default=TRUE)) { 
    args <- commandArgs(trailing=TRUE)
    query_id <- as.character(args[1])
    if (length(args) > 1) {
        return_obj <- as.logical(args[2]) 
    } else {
        return_obj <- FALSE
    }
    Run_compensate(query_id, return_obj) 
}



