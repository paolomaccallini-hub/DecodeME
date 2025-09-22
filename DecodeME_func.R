# file name: DecodeME_func
#
library(httr) # for API query
library(R.utils) # for .bgz files
library(data.table) # for data tables
library(gtexr) # GTEx 
library(otargen) # Open Target
library(reticulate) # for LD format conversion
library(Matrix)
library(LDlinkR)
library(susieR) # for fine mapping
library(yaml) # for settings reading
library(MungeSumstats)
#
#-------------------------------------------------------------------------------
# Read configurations from YAML file
#-------------------------------------------------------------------------------
#
config<-read_yaml("DecodeME_config.yml")
#
#-------------------------------------------------------------------------------
# Set Python environment
#-------------------------------------------------------------------------------
#
reticulate::use_python(config$python$path,required=T)
py_config()
#
#-------------------------------------------------------------------------------
# API keys
#-------------------------------------------------------------------------------
#
OpenGWAS_key<-config$api_keys$OpenGWAS
LDlinkR_key<-config$api_keys$LDlinkR
#
#-------------------------------------------------------------------------------
# Statistical settings and filters
#-------------------------------------------------------------------------------
#
pco_c<-as.numeric(config$filters$p_value_common) # p value cut-off for common myvariants
pco_uc<-as.numeric(config$filters$p_value_uncommon) # p value cut-off for uncommon myvariants
MAFco_c<-as.numeric(config$filters$maf_common) # lower cut-off for minor allele frequency of common myvariants
MAFco_uc<-as.numeric(config$filters$maf_uncommon) # lower cut-off for minor allele frequency of uncommon myvariants
PIPco<-as.numeric(config$filters$pip_cutoff) # remove variants with PIP below PIPco after finemapping
ABCco<-as.numeric(config$filters$abc_score_cutoff) # remove variants with ABC.score below ABCco after finemapping
exonic_variants <- c(
  "missense_variant",
  "synonymous_variant",
  "frameshift_variant",
  "stop_gained",
  "inframe_deletion",
  "start_lost",
  "inframe_insertion",
  "stop_lost",
  "stop_retained_variant",
  "protein_altering_variant",
  "coding_sequence_variant",
  "incomplete_terminal_codon_variant"
)
#
#-------------------------------------------------------------------------------
# Phenotypes (binary traits)
#-------------------------------------------------------------------------------
#
phenotypes<-config$phenotypes
#
#-------------------------------------------------------------------------------
# Radius of a locus
#-------------------------------------------------------------------------------
#
Radius<-config$locus$Radius
#
#-------------------------------------------------------------------------------
# GTEx version
#-------------------------------------------------------------------------------
#
gtex_version<-config$GTEx$version
#
#-------------------------------------------------------------------------------
# Sex-specific tissues to be removed (GTEx)
#-------------------------------------------------------------------------------
#
male_only_tissues<-c("Testis","Prostate") # Male-only tissues
female_only_tissues<-c("Ovary","Uterus","Vagina","Fallopian_Tube") # Female-only tissues
#
#-------------------------------------------------------------------------------
# Create output folder, if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"DecodeME_output")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
#
#-------------------------------------------------------------------------------
# Create input folders, if absent
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/DecodeME")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
current_dir<-getwd()
folder_path<-file.path(current_dir,"Data/LD")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
} 
#
#-------------------------------------------------------------------------------
# Build data base from DecodeME
#-------------------------------------------------------------------------------
#
# Download summary statistics, if not present
#
url<-"https://files.de-1.osf.io/v1/resources/rgqs3/providers/osfstorage/688881b3c8c15886a1e7cc93/?zip="
file_path<-file.path(current_dir,"Data/DecodeME/DecodeME_summary.zip")
if(!file.exists(file_path)) {
  print("Downloading the summary statistics...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404,403) # don't retry on these errors
  )
}   
#
#-------------------------------------------------------------------------------
# Build ABC database
#-------------------------------------------------------------------------------
#
url<-"https://mitra.stanford.edu/engreitz/oak/public/Nasser2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
destfile<-"ABC.tsv.bgz"
file_path<-file.path(current_dir,"Data",destfile)
if(!file.exists(file_path)) {
  print("Downloading and editing the full set ABC predictions. This may take a while...")
  RETRY(
    verb = "GET",
    url = url,
    write_disk(file_path, overwrite = TRUE),
    times = 5,           # up to 5 attempts
    pause_min = 5,       # wait 5s between attempts
    terminate_on = c(404, 403) # don't retry on these errors
  ) 
} 
print("Reading ABC database. This may take a while...")
ABC<-fread(file_path)
ABC<-ABC[ABC.Score>=ABCco] # https://www.nature.com/articles/s41586-021-03446-x
#
#-------------------------------------------------------------------------------
# Load RegulomeDB
#-------------------------------------------------------------------------------
#
print("Loading RegulomeDB database...")
file.name<-"ENCFF250UJY.tsv.gz"
#
if(file.exists(file.name)) {
  regulome<-fread(file.name)
  file.remove(file.name)
  regulome<-regulome[,c("rsid","ranking")]
  colnames(regulome)<-c("rsid","RegulomeDB")
  file.name<-file.path("Data","ENCFF250UJY.tsv.gz")
  fwrite(regulome,file.name) # write a lighter version of the file and gzip it
  file.remove("ENCFF250UJY.tsv.gz")  
} else if (file.exists(file.path("Data","ENCFF250UJY.tsv.gz"))) {
  regulome<-fread(file.path("Data","ENCFF250UJY.tsv.gz"))
} else {
  print("You have not downloaded the RegulomeDB score: Do it now!")
  print("Go here: https://www.regulomedb.org/regulome-search/")  
}
#
#-------------------------------------------------------------------------------
# Read the proper NPZ LD matrix
#-------------------------------------------------------------------------------
#
download_UKBB_LD<-function(chrom,position) {
  #
  # Round to 3Mb window (used by UKBB LD panels)
  #
  window_start<-floor((position-Radius)/3e6)*3e6+1
  window_end<-window_start+3e6
  #
  # Check
  #
  while (window_end<=position+Radius) {
    window_start<-window_start+1e6
    window_end<-window_start+3e6
  }
  #
  # Search for the proper LD matrix, if present at all
  #
  search.test<-0 # used to end search
  while(search.test==0) {
    #
    # File base name
    #
    base_name<-paste0("chr",chrom,"_",window_start,"_",window_end)
    url_base<-"https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/"
    #
    # Full URLs
    #
    urls<-paste0(url_base,base_name,c(".gz",".npz"))  
    #
    # Check if the .gz exists remotely before proceeding
    #
    check_url<-urls[2]
    http_status<-tryCatch(
      {
        h<-httr::HEAD(check_url)
        httr::status_code(h)
      },
      error=function(e) return(404)  # assume not found if request fails
    )
    if (http_status==404) {
      window_start<-window_start+1e6 # add 1e6 bases to check next LD
      window_end<-window_start+3e6
      if (!(position-Radius>=window_start&position+Radius<=window_end)) {
        print("I could not find a LD matrix for your lead SNP!")
        search.test<-1
      }
    } else {
      print("I found the proper LD matrix for your lead SNP!")
      search.test<-1
    }
  }
  #
  # If search succeeded, proceed with download
  #
  if (http_status!=404) {
    #
    # Destination paths
    #
    dest_file<-file.path("Data/LD",basename(urls))
    #
    # Download
    #
    for (i in seq_along(urls)) {
      if (!file.exists(dest_file[i])) {
        message("Downloading: ",urls[i])
        tryCatch(
          curl::curl_download(urls[i],dest_file[i]),
          error=function(e) stop("Failed to download: ",urls)
        )
      } else {
        message("Already exists: ",dest_file[i])
      }
    }
  }
  return(c(base_name,http_status))
}
#
#-------------------------------------------------------------------------------
# Convert the LD matrix from NPZ to RDS 
#-------------------------------------------------------------------------------
#
npz2rsd<-function(chrom,position,base_name) {
  #
  # Import numpy
  #
  np<-reticulate::import("numpy")
  scipy_sparse<-reticulate::import("scipy.sparse")
  #
  # Build paths to files
  #
  npz_path<-paste0("Data/LD/",base_name,".npz")
  gz_path<-paste0("Data/LD/",base_name,".gz")
  #
  # Load .npz file
  #
  print("Reading NPZ LD matrix...")
  ld<-np$load(npz_path,allow_pickle=T)
  #
  # Extract sparse matrix components
  #
  data<-reticulate::py_to_r(ld$f[["data"]])
  row<-reticulate::py_to_r(ld$f[["row"]])
  col<-reticulate::py_to_r(ld$f[["col"]])
  shape<-reticulate::py_to_r(ld$f[["shape"]])
  shape<-as.integer(shape)
  remove(ld)
  #
  # Construct a COO sparse matrix (becomes dgTMatrix in R)
  #
  print("Converting NPZ LD matrix to RDS...")
  ld_matrix<-scipy_sparse$coo_matrix(tuple(data,tuple(row,col)),shape=shape)
  remove(row,col,data)
  #
  # Convert dgTMatrix to CsparseMatrix
  #
  ld_matrix<-as(ld_matrix,"CsparseMatrix")
  diag(ld_matrix)<-1
  #
  # Read SNP names
  #
  snps<-fread(gz_path)
  if (nrow(snps)!=nrow(ld_matrix)) {
    warning("Number of SNPs and matrix dimensions do not match.")
  } else {
    dimnames(ld_matrix)<-list(snps$rsid,snps$rsid)
  }
  #
  # Create a variant column (we cannot select by rsid, because different variants have same rsid)
  #
  snps<-snps[,variant:=do.call(paste,c(.SD,sep=":")),.SDcols=c("chromosome","position","allele1","allele2")]
  #
  # Keep only variants included in our locus and in our GWAS
  #
  pos<-position
  snps_temp<-subset.data.frame(snps,snps$position>=pos-Radius&snps$position<=pos+Radius)
  snps_temp<-subset.data.frame(snps_temp,variant%in%mydata$variant)
  index<-which(snps$variant%in%snps_temp$variant)
  ld_matrix<-ld_matrix[index,index]
  #
  # Fill the upper triangle
  #
  ld_matrix<-forceSymmetric(ld_matrix,uplo="L")
  #
  # Save 
  #
  print("Saving converted matrix...")
  rds_out<-gsub(".npz","",npz_path)
  rds_out<-paste0(rds_out,"_",position,".RDS")
  saveRDS(ld_matrix,file=rds_out)
  message("Saved LD matrix as RDS to: ",rds_out)
}
#
#-------------------------------------------------------------------------------
# Fine Mapping with SusieR and LD matrix from UKBiobank
#-------------------------------------------------------------------------------
#
Fine_map_LD<-function(myvariants) {
  #
  # Search path to LD matrix
  #
  leadSNP<-which.min(myvariants$pval)
  position<-myvariants$pos[leadSNP]
  base_name<-myvariants$LD.name[1]
  LD_path<-paste0("Data/LD/",base_name,"_",position,".RDS")
  #
  # Load LD matrix (this is cut to the risk locus)
  #
  R<-readRDS(LD_path)
  R<-data.matrix(R)
  #
  # Retrieve all the variants included in both R and mydata
  #
  index<-which(mydata$rsid%in%rownames(R))
  mydata.2<-mydata[index,]
  index<-which(rownames(R)%in%mydata.2$rsid)
  R<-R[index,index]
  #
  # Calculate Z scores and correlation matrix for alternative alleles
  #
  z_scores<-mydata.2$beta/mydata.2$se # we calculate the z score
  #
  # Set the max number admitted of causal myvariants (L) and the sample size (n) 
  #
  n<-n_cases+n_controls # sample size
  phi<-n_cases/n_controls
  ne<-n*phi*(1-phi) # effective sample size
  L<-10 # max number of causal myvariants
  #
  # Use the function fitted_rss of the susieR to calculate 
  # the posterior inclusion probability (PIP) of each variant
  #
  fitted_rss=susie_rss(z_scores,R,n=ne,L=L,z_ld_weight=0)
  #
  # Summary of the results
  #
  SusieResult<-summary(fitted_rss)
  #
  # Collect relevant results
  #
  PIP<-SusieResult[[1]][order(SusieResult[[1]]$variable),]
  mydata.2$CS<-PIP$cs
  mydata.2$PIP<-PIP$variable_prob
  mydata.2<-as.data.table(mydata.2)
  mydata.2[,c("chr","pos","ref","alt"):=tstrsplit(variant,":",fixed=TRUE)] # split
  mydata.2$locus<-rep(myvariants$locus[1],nrow(mydata.2))
  #
  # plot and save it
  #
  RL<-myvariants$locus[1]
  file.name<-paste0(paste0("DecodeME_output/",pheno,"_",RL,".png"))
  png(file.name,res=100,width=1000,height=500)
  susie_plot(fitted_rss,xlab="",y="PIP",add_legend="topright",xaxt="n")  # suppress x-axis
  keep_idx<-seq(1,length(mydata.2$pos),by=50)  # or any spacing you want
  axis(
    side = 1,
    at = keep_idx,
    labels = mydata.2$pos[keep_idx],
    las = 2,         # Rotate labels vertically
    cex.axis = 0.6   # Shrink text size
  )
  title(main=paste0(pheno,", ","chr ",myvariants$chr[1]),cex.main=0.8)
  dev.off()
  #
  # Add lead SNP and save
  #
  mydata.2$LeadSNP<-rep(F,nrow(mydata.2))
  index<-which.min(mydata.2$pval)
  mydata.2$LeadSNP[index]<-T
  return(mydata.2)
}
#
#-------------------------------------------------------------------------------
# This function add RegulomeDB ranking
#-------------------------------------------------------------------------------
#
RegulomeDB<-function(myvariants) {
  #
  print("Adding RegulomeDB annotation, this may take a while...")
  myvariants<-merge(myvariants,regulome,by="rsid",all.x=T)
  #
  return(myvariants)
}
