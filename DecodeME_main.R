# file name: DecodeME_main
#
#-------------------------------------------------------------------------------
# This script performs GWAS analysis on DecodeME GWAS-1 summary statistics
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
# Link to module with functions, packages, and settings
#-------------------------------------------------------------------------------
#
BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
source("DecodeME_func.R",echo=F)
#
#-------------------------------------------------------------------------------
# Filter summary statistics
#-------------------------------------------------------------------------------
#
file_name<-"My_genes_DecodeME.csv"
if (file.exists(file_name)) file.remove("My_genes_DEcodeME.csv")
#
# Read variants that passed quality filter
#
input_phenotype<-"gwas_qced.var.gz"
zip_path<-file.path(current_dir,"Data/DecodeME","DecodeME_summary.zip")
tmpdir<-tempdir()
unzip(zip_path,files=input_phenotype,exdir=tmpdir)
myvariants<-fread(file.path(tmpdir,input_phenotype))
#
for (pheno in phenotypes) {
  file_name<-paste0("DecodeME_output/",pheno,"_finemapped.csv")
  if(!dir.exists(file_name)) {
    #
    #-------------------------------------------------------------------------
    # Read summary statistics 
    #-------------------------------------------------------------------------
    #
    print(paste0("Filtering phenotype ",pheno,"..."))
    #
    # Read summary statistics for selected sex and phenotype
    #
    input_phenotype<-paste0(pheno,".regenie.gz")
    zip_path<-file.path(current_dir,"Data/DecodeME","DecodeME_summary.zip")
    tmpdir<-tempdir()
    unzip(zip_path,files=input_phenotype,exdir=tmpdir)
    mydata<-fread(file.path(tmpdir,input_phenotype))
    #
    # Keep only variants that passed quality filter
    #
    mydata<-merge(myvariants,mydata,by="ID",all.x=T)
    #
    # Read sample size for selected sex and phenotype
    #
    n_cases<-mydata$N_CASES[1]
    n_controls<-mydata$N_CONTROLS[1]
    print(paste0("I found ",n_cases," cases and ",n_controls," controls."))
    #
    # Standardization and lift-over
    #
    munge_path<-file.path(current_dir,"Data/DecodeME",paste0(pheno,".tsv.gz"))
    #
    if (!file.exists(munge_path)) {
      format_sumstats(path=file.path(tmpdir,input_phenotype),ref_genome="GRCh38",
                      convert_ref_genome="GRCh37",save_path=munge_path)  
      mydata<-fread(munge_path)
    } else {
      mydata<-fread(munge_path)
    }
    #
    # Rewrite the ID column with GRCh37 coordinates
    #
    mydata[,ID:=paste(CHR,BP,A1,A2,sep=":")]
    #
    # Calculate MAF
    #
    mydata[,minor_AF:=ifelse(FRQ>0.5,1-FRQ,FRQ)]
    #
    # Build a proxy for INFO as 1/(SE^2 * Neff_cc * 2*MAF*(1-MAF)),
    # where Neff_cc = N * pi_hat * (1 - pi_hat)
    #
    p<-mydata$minor_AF
    pi_hat<-mydata$N_CAS/mydata$N
    Neff_cc<-mydata$N*pi_hat*(1-pi_hat)              
    den<-mydata$SE^2*Neff_cc*2*p*(1-p)
    INFO_proxy<-1/den
    #
    INFO_proxy[!is.finite(INFO_proxy)]<-NA
    INFO_proxy[INFO_proxy>1]<-1
    INFO_proxy[INFO_proxy<0]<-0
    #
    mydata$INFO_proxy<-INFO_proxy
    #
    # Filter by INFO_PROXY
    #
    mydata<-mydata[INFO_proxy>INFOco]
    #
    # Edit the names of the columns to match Neale Lab standard
    #
    mydata<-subset.data.frame(mydata,select=-CHR)
    mydata<-subset.data.frame(mydata,select=-BP)
    mydata<-subset.data.frame(mydata,select=-A1)
    mydata<-subset.data.frame(mydata,select=-A2)
    mydata<-subset.data.frame(mydata,select=-LOG10P)
    #
    for (j in 1:ncol(mydata)) {
      if (colnames(mydata)[j]=="ID") colnames(mydata)[j]<-"variant"
      if (colnames(mydata)[j]=="SNP") colnames(mydata)[j]<-"rsid"
      if (colnames(mydata)[j]=="N") colnames(mydata)[j]<-"n_complete_samples"
      if (colnames(mydata)[j]=="BETA") colnames(mydata)[j]<-"beta"
      if (colnames(mydata)[j]=="SE") colnames(mydata)[j]<-"se"
      if (colnames(mydata)[j]=="P") colnames(mydata)[j]<-"pval"
      if (colnames(mydata)[j]=="N_CAS") colnames(mydata)[j]<-"n_cases"
      if (colnames(mydata)[j]=="N_CON") colnames(mydata)[j]<-"n_controls"
    }
    #
    #-------------------------------------------------------------------------
    # Search for variants significantly different between patients and controls
    #-------------------------------------------------------------------------
    #
    # Filter variants by p value and MAF: keep variants with p value below pco_c
    # and MAF above MAFco_uc (uncommon variants will be filtered later using pco_uc,
    # now we keep them for fine mapping)
    #
    print("I am searching for variants significantly associated with phenotype...")
    myvariants<-subset.data.frame(mydata,minor_AF>=MAFco_c&pval<=pco_c)
    #
    # Exit if no variants remain
    #
    if (nrow(myvariants)==0) next
    #
    # Remove uncommon variants with p value above pco_uc
    #
    temp<-subset.data.frame(myvariants,minor_AF<MAFco_c) # all uncommon variants
    if(nrow(temp)>0) {
      temp<-subset.data.frame(temp,pval>pco_uc) # all non-significant uncommon variants 
      if (nrow(temp)>0) {
        for (vi in 1:nrow(temp)) {                 
          myvariants<-subset.data.frame(myvariants,variant!=temp$variant[vi])  
        }
      } 
    }
    #
    # Exit if no variant remains
    #
    if (nrow(myvariants)==0) next
    #
    #-------------------------------------------------------------------------
    # Edit 
    #-------------------------------------------------------------------------
    #
    setDT(myvariants)
    myvariants[,c("chr","pos","ref","alt"):=tstrsplit(variant,":",fixed=TRUE)] # split
    #
    #-------------------------------------------------------------------------
    # Find loci
    #-------------------------------------------------------------------------
    #
    myvariants$locus<-rep(NA_integer_,nrow(myvariants))
    myvariants$pos<-as.numeric(myvariants$pos)
    uniqueCHR<-unique(myvariants$chr)
    locus<-1
    for(chri in uniqueCHR){
      temp<-myvariants[chr==chri]
      temp<-temp[order(temp$pval,decreasing=F),]
      if(nrow(temp)>0){
        assigned<-rep(FALSE,nrow(temp))
        for(i in 1:nrow(temp)){
          if(!assigned[i]){
            window<-(abs(temp$pos-temp$pos[i])<=Radius)&!assigned
            temp$locus[window]<-locus
            assigned[window]<-TRUE
            locus<-locus+1
          }
        }
        myvariants[chr==chri,locus:=temp$locus]
      }
    }
    #
    #-------------------------------------------------------------------------
    # Prepare LD matrices for lead variants
    #-------------------------------------------------------------------------
    #
    myvariants$LD.name<-rep(NA,nrow(myvariants))
    loci<-unique(myvariants$locus)
    for (lc in loci) {
      temp<-myvariants[locus==lc]
      leadSNP<-which.min(temp$pval)
      chrom<-temp$chr[leadSNP]
      position<-temp$pos[leadSNP]
      #
      # Download NPZ LD matrix (if not already there)
      #
      result<-download_UKBB_LD(chrom,position)
      #
      # Convert to RDS format and subset around leadSNP
      #
      if (result[2]==404) {
        index<-(which(myvariants$variant%in%temp$variant))
        myvariants$LD.name[index]<-"none"
        break
      } else {
        index<-(which(myvariants$variant%in%temp$variant))
        myvariants$LD.name[index]<-result[1]
        base_name<-result[1]
        npz2rsd(chrom,position,base_name) 
      }
    }
    #
    # Remove variants that cannot be fine mapped for lack of LD matrix
    #
    myvariants<-subset.data.frame(myvariants,LD.name!="none")
    #
    # Exit if no variant remains
    #
    if (nrow(myvariants)==0) next
    #
    #-------------------------------------------------------------------------
    # Finemapping with SusieR
    #-------------------------------------------------------------------------
    #
    mylist<-list()
    loci<-unique(myvariants$locus)
    for (lc in loci) {
      temp<-subset.data.frame(myvariants,locus==lc)
      print(paste0("Finemapping locus ",lc))
      mylist[[lc]]<-Fine_map_LD(temp)
    }
    myvariants<-do.call(rbind,mylist)
    myvariants<-subset.data.frame(myvariants,PIP>PIPco&CS!="-1")
    #
    # Exit if no variant remains
    #
    if (nrow(myvariants)==0) next
    #
    #-------------------------------------------------------------------------
    # Add variants in high LD with top SNPs and not in the original GWAS
    #-------------------------------------------------------------------------
    #
    print("Adding variants in high LD with lead SNPs...")
    temp<-myvariants # a temporary copy
    for (snp in 1:nrow(temp)) {
      if (temp$LeadSNP[snp]) {
        LDblock<-LDproxy(temp$rsid[snp],pop="GBR",r2d="r2",token=LDlinkR_key,file=F,
                         genome_build="grch37",win_size="500000",
                         api_root="https://ldlink.nih.gov/LDlinkRest")
        #
        if (is.null(LDblock)|ncol(LDblock)==1) next
        #
        # Remove variants included in our GWAS
        #
        index<-which(LDblock$RS_Number%in%mydata$rsid)
        LDblock<-LDblock[-index,] 
        if (nrow(LDblock)==0) next
        #
        # Remove variants already added
        #
        index<-which(LDblock$RS_Number%in%myvariants$rsid)
        LDblock<-LDblock[-index,] 
        if (nrow(LDblock)==0) next
        #
        # Select by R2
        #
        LDblock<-subset.data.frame(LDblock,R2>=0.8)
        if (nrow(LDblock)==0) next
        #
        # Edit LDblock 
        #
        for (p in 1:nrow(LDblock)) {
          LDblock$Coord[p]<-gsub("chr","",LDblock$Coord[p])
          alleles<-gsub("\\(","",LDblock$Alleles[p])
          alleles<-gsub("\\)","",alleles)
          alleles<-strsplit(alleles,"/")[[1]]
          LDblock$Coord[p]<-paste0(c(LDblock$Coord[p],alleles[1],alleles[2]),collapse=":")
        }
        LDblock<-subset.data.frame(LDblock,select=c("Coord","RS_Number","MAF","R2"))
        colnames(LDblock)<-c("variant","rsid","minor_AF","R2")
        LDblock$locus<-rep(temp$locus[snp],nrow(LDblock))
        LDblock$CS<-rep(temp$CS[snp],nrow(LDblock))
        #
        # Add these new variants
        #
        myvariants<-merge(myvariants,LDblock,all.x=T,all.y=T)
      }
    }
    #
    #-------------------------------------------------------------------------
    # Add coordinates of variants on GRCh38
    #-------------------------------------------------------------------------
    #
    myvariants$var.GRCh38<-rep(NA,nrow(myvariants))
    for (vi in 1:nrow(myvariants)) {
      variant<-get_variant(snpId=myvariants$rsid[vi],.verbose=F) 
      if(nrow(variant)>0) {
        myvariants$var.GRCh38[vi]<-variant$variantId # GTEx variant identifier  
      }
    }
    #
    #-------------------------------------------------------------------------
    # Search for significant variant-gene association using eQTL 
    #-------------------------------------------------------------------------
    #
    print("Collecting significant variant-gene associations by eQTL")
    myvariants$eQTL.gene<-rep(NA,nrow(myvariants))
    myvariants$eQTL.tissue<-rep(NA,nrow(myvariants))
    myvariants$eQTL.nes<-rep(NA,nrow(myvariants))
    for (vi in 1:nrow(myvariants)) {
      if (!is.na(myvariants$var.GRCh38[vi])) {
        result<-get_significant_single_tissue_eqtls(variantIds=myvariants$var.GRCh38[vi],
                                                    datasetId=gtex_version,.verbose=F)
      } else {
        result<-get_significant_single_tissue_eqtls(variantIds=myvariants$rsid[vi],.verbose=F)
      }
      #
      # Remove sex specific tissues
      #
      if(nrow(result)>0) {
        index<-which(result$tissueSiteDetailId%in%c(male_only_tissues,
                                                    female_only_tissues))
        if (length(index)>0) result<-result[-index,]
      }
      #
      # Add gene description and keep only protein coding genes
      #
      if(nrow(result)>0) {
        a<-get_genes(result$geneSymbol,.verbose=F)
        if (nrow(a)>0) {
          a<-subset.data.frame(a,select=c("geneSymbol","geneType"))
          result<-merge(result,a,by="geneSymbol",all.x=T)
          result<-subset.data.frame(result,geneType=="protein coding")
        }
      }
      #
      # Collaspe results
      #
      if(nrow(result)>0) {
        myvariants$eQTL.gene[vi]<-paste0(result$geneSymbolUpper,collapse="/")
        myvariants$eQTL.tissue[vi]<-paste0(result$tissueSiteDetailId,collapse="/")
        myvariants$eQTL.nes[vi]<-paste0(result$nes,collapse="/")
      }
    }
    #
    #-------------------------------------------------------------------------
    # Small edit 
    #-------------------------------------------------------------------------
    #
    myvariants[myvariants==""]<-NA
    #
    #-------------------------------------------------------------------------
    # Add missing VEP annotation (Open Target) 
    #-------------------------------------------------------------------------
    #
    myvariants$consequence<-rep(NA,nrow(myvariants))
    for (vi in 1:nrow(myvariants)) {
      if (is.na(myvariants$consequence[vi])) {
        if (!is.na(myvariants$var.GRCh38[vi])) {
          variant<-gsub("chr","",myvariants$var.GRCh38[vi])
          variant<-gsub("_b38","",variant)
        } else {
          variant<-myvariants$rsid[vi]
        } 
        result<-try(variantEffectQuery(variant),silent=T)
        if (!is.null(result)) {
          result<-subset.data.frame(result,method=="VEP")
          myvariants$consequence[vi]<-result$assessment[1]
        }
      }
    }
    #
    #-------------------------------------------------------------------------
    # Add closest gene (Open Target) 
    #-------------------------------------------------------------------------
    #
    myvariants$closest.gene<-rep(NA,nrow(myvariants))
    for (vi in 1:nrow(myvariants)) {
      if (!is.na(myvariants$var.GRCh38[vi])) {
        variant<-gsub("chr","",myvariants$var.GRCh38[vi])
        variant<-gsub("_b38","",variant)
      } else {
        variant<-myvariants$rsid[vi]
      } 
      result<-try(variantEffectPredictorQuery(variant),silent=T)
      if (!is.null(result)) {
        index<-which.min(abs(result$distanceFromTss))
        myvariants$closest.gene[vi]<-result$target.approvedSymbol[index]
      }
    }
    #
    #-------------------------------------------------------------------------
    # Search for overlaps on regulatory regions using ABC 
    # (ABC file is in GRCh37, you can verifiy this by checking last position on chr4!)
    #-------------------------------------------------------------------------
    #
    print("Collecting significant variant-gene associations by ABC")
    myvariants$ABC.gene<-rep(NA,nrow(myvariants))
    myvariants$ABC.CellType<-rep(NA,nrow(myvariants))
    myvariants$ABC.score<-rep(NA,nrow(myvariants))
    myvariants$ABC.class<-rep(NA,nrow(myvariants))
    for (vi in 1:nrow(myvariants)) {
      if (!is.na(myvariants$var.GRCh38[vi])) {
        chrom<-myvariants$chr[vi]
        pos<-myvariants$pos[vi]
        temp<-ABC[chr==chrom&start<=pos&end>=pos]
        #
        # Add gene description and keep only protein-coding genes
        #
        if(nrow(temp)>0) {
          a<-get_genes(temp$TargetGene,.verbose=F)
          if (nrow(a)>0) {
            a<-subset.data.frame(a,select=c("geneSymbol","geneType"))
            a<-subset.data.frame(a,geneType=="protein coding")
            index<-which(temp$TargetGene%in%a$geneSymbol)
            temp<-temp[index,]
          } else {
            temp<-data.frame()
          }
        }
        #
        # Write results
        #
        if (nrow(temp)>0) {
          myvariants$ABC.gene[vi]<-paste0(temp$TargetGene,collapse="/")
          myvariants$ABC.CellType[vi]<-paste0(temp$CellType,collapse="/")
          myvariants$ABC.score[vi]<-paste0(temp$ABC.Score,collapse="/")
          myvariants$ABC.class[vi]<-paste0(temp$class,collapse="/")
        }
      }
    }
    #
    #-------------------------------------------------------------------------
    # Add regulomeDB score
    #-------------------------------------------------------------------------
    #
    myvariants<-RegulomeDB(myvariants)
    #
    #-------------------------------------------------------------------------
    # Save results
    #-------------------------------------------------------------------------
    #
    file_name<-paste0("DEcodeME_output/",pheno,"_finemapped",".csv")
    fwrite(myvariants,file=file_name,sep=";")
    gc() # free unused memory
    #
    #-------------------------------------------------------------------------
    # Add results to final output
    #-------------------------------------------------------------------------
    #
    file_name<-"My_genes_DecodeME.csv"
    next.row<-data.frame(name="",NCBI.id="",list.name="",weight="",Gene_type="",Study_type="",
                         Phenotype="",Cases="",Sex="",Designation="",Description="",
                         Variant="",
                         Tissues="",reference="",stringsAsFactors=FALSE)
    if (file.exists(file_name)) {
      mygenes<-read.csv(file_name,sep=";",header=T)
    } else {
      mygenes<-next.row
    } 
    #
    for (vi in 1:nrow(myvariants)) {
      #
      # collect genes mapped by exonic variants
      #
      if (myvariants$consequence[vi]%in%exonic_variants) {
        mygenes<-rbind(mygenes,next.row)
        counter<-nrow(mygenes)
        mygenes$name[counter]<-myvariants$closest.gene[vi]
        mygenes$Description[counter]<-myvariants$consequence[vi]
        mygenes$Variant[counter]<-myvariants$rsid[vi]
        mygenes$weight[counter]<-1
        #
        a<-get_genes(mygenes$name[counter])
        if (nrow(a)>0) {
          mygenes$Gene_type[counter]<-a$geneType
          mygenes$NCBI.id[counter]<-a$entrezGeneId  
        } else {
          mygenes$Gene_type[counter]<-NA
          mygenes$NCBI.id[counter]<-NA  
        }
        mygenes$list.name[counter]<-paste0("DecodeME_",pheno,"_",sex)
        mygenes$Study_type[counter]<-"GWAS"
        mygenes$Phenotype[counter]<-pheno
        mygenes$Sex[counter]<-NA
        mygenes$Cases[counter]<-n_cases
        mygenes$Designation[counter]<-"DecodeME"
      }
      #
      # collect genes mapped by eQTL
      #
      if (!is.na(myvariants$eQTL.gene[vi])) {
        temp1<-data.frame(name=strsplit(myvariants$eQTL.gene[vi],"/")[[1]],
                          tissue=strsplit(myvariants$eQTL.tissue[vi],"/")[[1]])
        unique.genes<-unique(temp1$name)
        for (ug in unique.genes) {
          temp2<-subset.data.frame(temp1,name==ug)
          mygenes<-rbind(mygenes,next.row)
          counter<-nrow(mygenes)
          mygenes$name[counter]<-ug
          mygenes$Tissues[counter]<-paste0(temp2$tissue,collapse="/")
          mygenes$Description[counter]<-"eQTL"
          mygenes$Variant[counter]<-myvariants$rsid[vi]
          mygenes$weight[counter]<-1
          #
          a<-get_genes(mygenes$name[counter])
          if (nrow(a)>0) {
            mygenes$Gene_type[counter]<-a$geneType
            mygenes$NCBI.id[counter]<-a$entrezGeneId  
          } else {
            mygenes$Gene_type[counter]<-NA
            mygenes$NCBI.id[counter]<-NA  
          }
          mygenes$list.name[counter]<-paste0("DecodeME_",pheno)
          mygenes$Study_type[counter]<-"GWAS"
          mygenes$Phenotype[counter]<-pheno
          mygenes$Sex[counter]<-NA
          mygenes$Cases[counter]<-n_cases
          mygenes$Designation[counter]<-"DecodeME"
        }
      }
      #
      # collect genes mapped by Activity By Contact
      #
      if (!is.na(myvariants$ABC.gene[vi])) {
        temp1<-data.frame(name=strsplit(myvariants$ABC.gene[vi],"/")[[1]],
                          tissue=strsplit(myvariants$ABC.CellType[vi],"/")[[1]])
        unique.genes<-unique(temp1$name)
        for (ug in unique.genes) {
          temp2<-subset.data.frame(temp1,name==ug)
          mygenes<-rbind(mygenes,next.row)
          counter<-nrow(mygenes)
          mygenes$name[counter]<-ug
          mygenes$Tissues[counter]<-paste0(temp2$tissue,collapse="/")
          mygenes$Description[counter]<-"ABC"
          mygenes$Variant[counter]<-myvariants$rsid[vi]
          mygenes$weight[counter]<-1
          #
          a<-get_genes(mygenes$name[counter])
          if (nrow(a)>0) {
            mygenes$Gene_type[counter]<-a$geneType
            mygenes$NCBI.id[counter]<-a$entrezGeneId  
          } else {
            mygenes$Gene_type[counter]<-NA
            mygenes$NCBI.id[counter]<-NA  
          }
          mygenes$list.name[counter]<-paste0("DecodeME_",pheno)
          mygenes$Study_type[counter]<-"GWAS"
          mygenes$Phenotype[counter]<-pheno
          mygenes$Sex[counter]<-NA
          mygenes$Cases[counter]<-n_cases
          mygenes$Designation[counter]<-"DecodeME"
        }
      }
    }
    #
    # Edit 
    #
    index<-which(mygenes$name=="")
    if (length(index)>0) mygenes<-mygenes[-index,]
    #
    # Save 
    #
    file_name<-"My_genes_DEcodeME.csv"
    write.table(mygenes,file_name,sep=";",row.names=F,col.names=T)
  }
} 





