
vcf_file_name="/cbscratch/franco/datasets/gtex_v8/genotypes/vcfs/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.01_chr22.vcf.gz"
ds_file_name="/home/franco/cbscratch/datasets/gtex_v8/genotypes/vcfs_0.05/GTEX_v8_2019-07-29_WGS_838Indiv_Freeze_NoMissingGT_SNPfilter_MAF0.05_chr22.dosage.gz"
donor_file="/home/franco/cbscratch/datasets/gtex_v8/genotypes/gtex_v8.sample"

parse_gt<-function(gt, gt_index){
  dosage = NA
  GTstr = strsplit(gt, ":")[[1]][gt_index]
  if (GTstr == "0/0") { dosage = 0.0 }
  if (GTstr == "1/0" || GTstr == "0/1" ) { dosage = 1.0 }
  if (GTstr == "1/1") { dosage = 2.0 }
  return(dosage)
}

parse_ds<-function(gt, ds_index){
  dosage = as.numeric(strsplit(gt, ":")[[1]][ds_index])
  return(dosage)
}

read_vcf_new<-function(vcf_file_name) {
  
  con = file(vcf_file_name, "r")
  gt_flag_start = F
  use_GT = F
  use_DS = F
  lines = 0
  dosage_matrix = c()
  snpspos = c()
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      colnames(snpspos) = c("snpid","chr","pos")
      snpspos = as.data.frame(snpspos)
      snpspos[,3] = as.numeric(as.character(snpspos[,3]))
      snpspos[,2] = as.character(snpspos[,2])
      snpspos[,1] = as.character(snpspos[,1])
      colnames(dosage_matrix)<-donors
      break
    }
    # if ( lines > 100 ) {
    #   colnames(snpspos) = c("snpid","chr","pos")
    #   snpspos = as.data.frame(snpspos)
    #   snpspos[,3] = as.numeric(as.character(snpspos[,3]))
    #   snpspos[,2] = as.character(snpspos[,2])
    #   snpspos[,1] = as.character(snpspos[,1])
    #   colnames(dosage_matrix)<-donors
    #   break
    # }
    # lines = lines + 1

    if (gt_flag_start) {
      
      gt_arr = strsplit(line, "\t")[[1]]
      chrom = gt_arr[1]
      if (! grepl("chr", chrom)) {
        chrom = paste("chr", chrom, sep="")
      }
      pos   = gt_arr[2]
      snpid = gt_arr[3]
      snpspos = rbind(snpspos, c(snpid, chrom, pos))
      
      if (! use_GT && ! use_DS) {
        format = strsplit(gt_arr[9], ":")[[1]]
        gt_index = grep("GT", format)
        if (length(gt_index) == 0 ) {
          ds_index = grep("DS", format)
          if ( length(ds_index) == 0) {
            message("Error, no GT or DS fields found in VCF file")
            break
          } else { use_DS = T}
        } else { use_GT = T }
      }
      
      gts = gt_arr[10:length(gt_arr)]
      if (use_GT){
        dosages = sapply(gts, FUN=parse_gt, gt_index = gt_index)
      } 
      if (use_DS){
        dosages = sapply(gts, FUN=parse_gt, ds_index = ds_index)
      }
      
      dosage_matrix = rbind(dosage_matrix, as.numeric(dosages))
    } else { 
      if ( grepl("#CHROM", line) ) {
        header = strsplit(line, "\t")[[1]]
        donors = header[10:length(header)]  
        print(paste("Number of samples: ",length(donors)))
        gt_flag_start = T
      }
    }
  }
  return(list(dosage_matrix, snpspos))
}



read_vcf <- function (vcf_file_name, maf_filter=T, maf_thres=0.05) {
  message("Reading Genotype ...")
  message(vcf_file_name)
  
  # read two times the vcf file, first for the columns names, second for the data
  # Read the vcf file to a string
  tmp_vcf_strings = readLines(vcf_file_name)
  header_colnum = grep("#CHROM", tmp_vcf_strings)
  total_colnum  = length(tmp_vcf_strings)
  data_colnum   = total_colnum - header_colnum
  
  header_string = tmp_vcf_strings[header_colnum] #tmp_vcf_strings[-(header_colnum + 1):-total_colnum]
  colnames<-unlist(strsplit(header_string,"\t"))
  
  snps_data = read.table(text="", col.names = colnames)
  colnames(snps_data) = colnames
  data_strings  = tmp_vcf_strings[ (header_colnum + 1): total_colnum]
  
  for (i in 1 : data_colnum) {
    tmpstrvals = unlist(strsplit(data_strings[i], "\t"))
    for (j in 1:9) {
      snps_data[i, j] = tmpstrvals[j]
    }
    for (j in 10:length(colnames)) {
      GT = unlist(strsplit(tmpstrvals[j], ":"))[1]
      dosage = 0.0
      if (GT == "./.") { dosage = as.numeric(unlist(strsplit(tmpstrvals[j], ":")[3]))}
      if (GT == "0/0") { dosage = 0.0 }
      if (GT == "1/0" || GT == "0/1" ) { dosage = 1.0 }
      if (GT == "1/1") { dosage = 2.0 }
      snps_data[i, j] = dosage
    }
  }
  
  # get SNPs positions for cis and trans analysis (before cropping the snp matrix)
  row_names = snps_data[,3]
  snpspos = snps_data[,c(3,1,2)]
  if ( !grepl("chr", snps_data[1, 1]) )
  { snpspos[,2] = paste("chr", snpspos[,2], sep="") }
  colnames(snpspos) = c("snpid","chr","pos")
  snpspos[,3] = as.numeric(snpspos[,3])
  
  snps_data = snps_data[, 10:ncol(snps_data)]
  rownames(snps_data) = row_names
  
  snps_mat = as.matrix(snps_data)
  
  return (list(snps_mat, snpspos))
}



read_genotype <- function(SNP_file_name, donors_file_name, maf_filter=T, maf_thres=0.01) {
  message("Reading Genotype ...")
  message(SNP_file_name)
  snps_mat = read.csv(file=SNP_file_name, sep=" ", stringsAsFactors=F, header=F)
  row_names = snps_mat[,2]
  donor_ids = read.csv(file=donors_file_name, sep=" ", stringsAsFactors=F, header=F, skip=2)[,1]
  
  # get SNPs positions for cis and trans analysis (before cropping the snp matrix)
  snpspos = snps_mat[,c(2,1,3)]
  snpspos[,2] = paste("chr", snpspos[,2], sep="")
  colnames(snpspos) = c("snpid","chr","pos")
  snpspos[,3] = as.numeric(snpspos[,3])
  
  mafs = snps_mat[,6]
  
  snps_mat = snps_mat[,7:ncol(snps_mat)]
  rownames(snps_mat) = row_names
  colnames(snps_mat) = donor_ids
  
  if (maf_filter) {
    mafs2remove = mafs < maf_thres | mafs > (1-maf_thres)
    mafs      = mafs[!mafs2remove]
    snps_mat  = snps_mat[!mafs2remove, ]
    snpspos   = snpspos[!mafs2remove,]
  }
  snps_mat = as.matrix(snps_mat)
  return(list(snps_mat, snpspos))
}


# start_time <- Sys.time()
# res1 = read_vcf_new(vcf_file_name)
# end_time <- Sys.time()
# diff = end_time - start_time
# print(diff)
# 
# 
# start_time <- Sys.time()
# res2 = read_vcf(vcf_file_name)
# end_time <- Sys.time()
# diff = end_time - start_time
# print(diff)
