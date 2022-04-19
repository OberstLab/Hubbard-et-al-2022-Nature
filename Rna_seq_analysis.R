renv::restore()

library(tximport)
library(DESeq2)
library(ggthemes)
library(latex2exp)
library(ggrepel)
library(ggdist)
library(cowplot)
library(Hmisc)
library(tidyverse)

# Data Location
h5_folder     = 'data/GSE200854_RAW'
metadata_file = 'data/sequence_metadata.tsv'
output_dir    = 'results'

# Analysis Parameters
DESeq2_count_prefilter = 10
DE_pcuttoff            = 0.01
DE_fc_cuttoff          = 1
proportional_dfc_cut   = 0.5

# Plot Parameters
cleveland_gene_count = 30
num_go_terms         = 30

# Plot aesthetics
adar_colour = tableau_color_pal()(10)[1]
zbp1_colour = tableau_color_pal()(10)[2]
pkr_colour  = tableau_color_pal()(10)[4]

point_alpha = 0.5
cleveland_dot_size = 3
colour_map = c(ADAR=adar_colour, ADAR_PKR_KO=pkr_colour, ADAR_ZBP1_KO=zbp1_colour, NS='grey30')

# Generate Analysis Metadata
sra_metadata = read_tsv(metadata_file) %>% as.data.frame()
study_group_dict = list( genotype       = c('P195A/p150-'              , 'P195A/p150- PKR-/-'                      , 'ZBP1-/-::ADARp150/P195A'                , 'WT'  )
                       , genotype_code  = c('ADAR'                     , 'ADAR_PKR_KO'                             , 'ADAR_ZBP1_KO'                           , 'WT'  )
                       , genotype_latex = c('$Adar1^{P195A/p150null}$' , '$Adar1^{P195A/p150null}::Eif2ak2^{-/-}$' , '$Adar1^{P195A/p150null}::Zbp1-a^{-/-}$' , '$WT$') ) %>% 
	as.data.frame() %>%
	mutate(genotype_code=as.character(genotype_code))
study_group_code2latex = with(study_group_dict, setNames(genotype_latex,genotype_code))
sra_metadata = left_join(sra_metadata,study_group_dict,by='genotype')
rownames(sra_metadata) = sra_metadata$`library name`
metadata = data.frame(condition=sra_metadata$genotype_code, sex=sra_metadata$sex, row.names=rownames(sra_metadata)) 
metadata =  metadata[sort(rownames(metadata)),]

# Read Kallisto Data
h5_files = paste0(h5_folder,'/',sra_metadata$gsm,'_',rownames(sra_metadata),'_L001_001.abundance.h5') %>% 
	set_names(rownames(sra_metadata))

# NOTE: This gets around the fact that tximport checks the file is called abundance.h5 not the extension of the file given
# TAKEN FROM: https://github.com/mikelove/tximport/blob/master/R/helper.R

# code contributed from Andrew Morgan
read_kallisto_h5 <- function(fpath, ...) {
	if (!requireNamespace("rhdf5", quietly=TRUE)) {
		stop("reading kallisto results from hdf5 files requires Bioconductor package `rhdf5`")
	}
	counts <- rhdf5::h5read(fpath, "est_counts")
	ids <- rhdf5::h5read(fpath, "aux/ids")
	efflens <- rhdf5::h5read(fpath, "aux/eff_lengths")
	
	# as suggested by https://support.bioconductor.org/p/96958/#101090
	ids <- as.character(ids)
	
	stopifnot(length(counts) == length(ids)) 
	stopifnot(length(efflens) == length(ids))
	
	result <- data.frame(target_id = ids,
						 eff_length = efflens,
						 est_counts = counts,
						 stringsAsFactors = FALSE)
	normfac <- with(result, (1e6)/sum(est_counts/eff_length))
	result$tpm <- with(result, normfac*(est_counts/eff_length))
	return(result)
}

kallisto_data_prelim = tximport(h5_files, type='kallisto', txOut=TRUE, importer=read_kallisto_h5 )
ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensembl_conversions = biomaRt::getBM( filters='ensembl_transcript_id_version'
                                    , attributes=c('ensembl_transcript_id_version','ensembl_gene_id','ensembl_gene_id_version','external_gene_name','entrezgene_id')
                                    , values=kallisto_data_prelim$counts %>% rownames()
                                    , mart=ensembl)
tx2gene = ensembl_conversions %>% select(ensembl_transcript_id_version, ensembl_gene_id_version)
kallisto_data = tximport(h5_files, type='kallisto', tx2gene=tx2gene, importer=read_kallisto_h5)

print('Running standard DESeq2 pipeline')
# Run DESeq2
dds = DESeqDataSetFromTximport(kallisto_data, sra_metadata, ~genotype_code) %>% 
	subset(rowSums(counts(.)) > DESeq2_count_prefilter) %>% 
	DESeq()

run_de = function(contrast){
	res = results(dds, contrast=contrast)
	res['ensembl_gene_id_version'] = rownames(res)
	res_df = res %>% 
		as_tibble() %>% 
		left_join(ensembl_conversions %>% select(-ensembl_transcript_id_version) %>% distinct(), by='ensembl_gene_id_version')
	res_df
}
res_adar_v_wt = run_de(c('genotype_code','ADAR'        ,'WT'))
res_pkr_v_wt  = run_de(c('genotype_code','ADAR_PKR_KO' ,'WT'))
res_zbp1_v_wt = run_de(c('genotype_code','ADAR_ZBP1_KO','WT'))

#Volcano Plots
create_volcano = function(res,condition){
	title = TeX(study_group_code2latex[condition])
	res_va = res %>% mutate(logP=-log10(padj), colour_tag=ifelse(abs(log2FoldChange) > DE_fc_cuttoff & padj < DE_pcuttoff, condition, 'NS'))
	lfc_bound = max(abs(res_va$log2FoldChange)) + 1
	p_volcano_combined = ggplot(res_va, aes(x=log2FoldChange,y=logP, colour=colour_tag, label=external_gene_name)) +
		geom_point(alpha=point_alpha, size=1) +
		geom_vline(xintercept=c(-DE_fc_cuttoff,DE_fc_cuttoff), linetype='dashed') +
		geom_hline(yintercept=-log10(DE_pcuttoff), linetype='dashed') +
		geom_text_repel(aes(x=log2FoldChange,y=logP, label=external_gene_name), size=6, colour='black', data=res_va %>% subset(colour_tag != 'NS')) +
		xlim(c(-lfc_bound,lfc_bound)) +
		scale_colour_manual(values=colour_map) +
		xlab(bquote(~Log[2]~ 'fold change vs WT')) +
		ylab(bquote(~-Log[10]~italic(P[adj]))) +
		theme(legend.position='none')
	p_volcano_combined
}

res_combined_v_wt = bind_rows(list(ADAR=res_adar_v_wt,ADAR_PKR_KO=res_pkr_v_wt,ADAR_ZBP1_KO=res_zbp1_v_wt), .id='genotype_code') %>% 
	mutate(logP=-log10(padj), colour_tag=ifelse(abs(log2FoldChange) > DE_fc_cuttoff & padj < DE_pcuttoff, genotype_code, 'NS'))
lfc_bound = max(abs(res_combined_v_wt$log2FoldChange)) + 1
pval_bound = max(abs(res_combined_v_wt$logP),na.rm=TRUE) + 1
study_group_code2latex_pre = lapply(study_group_code2latex, TeX)
labeller_function = function(y) { return(lapply(y, function(x) study_group_code2latex_pre[x])) }
p_volcano_combined = ggplot(res_combined_v_wt, aes(x=log2FoldChange,y=logP, colour=colour_tag, label=external_gene_name)) +
	geom_point(alpha=point_alpha, size=1) +
	geom_vline(xintercept=c(-DE_fc_cuttoff,DE_fc_cuttoff), linetype='dashed') +
	geom_hline(yintercept=-log10(DE_pcuttoff), linetype='dashed') +
	geom_text_repel(aes(x=log2FoldChange,y=logP, label=external_gene_name), colour='black', data=res_combined_v_wt %>% subset(colour_tag != 'NS')) +
	facet_grid(cols=vars(genotype_code), labeller=labeller_function ) +
	xlim(c(-lfc_bound,lfc_bound)) +
	scale_colour_manual(values=colour_map) +
	xlab(bquote(~Log[2]~ 'fold change vs WT')) +
	ylab(bquote(~-Log[10]~italic(P[adj]))) +
	theme(legend.position='none')
p_volcano_adar_v_wt = create_volcano(res_adar_v_wt, 'ADAR'        )
p_volcano_pkr_v_wt  = create_volcano(res_pkr_v_wt , 'ADAR_PKR_KO' )
p_volcano_zbp1_v_wt = create_volcano(res_zbp1_v_wt, 'ADAR_ZBP1_KO')

dir.create(output_dir, recursive=TRUE)
ggsave(paste0(output_dir,'/DE_Volcano_adar_v_wt.pdf'), p_volcano_adar_v_wt + xlim(c(-lfc_bound,lfc_bound)) + ylim(c(0,pval_bound)))
ggsave(paste0(output_dir,'/DE_Volcano_pkr_v_wt.pdf' ), p_volcano_pkr_v_wt  + xlim(c(-lfc_bound,lfc_bound)) + ylim(c(0,pval_bound)))
ggsave(paste0(output_dir,'/DE_Volcano_zbp1_v_wt.pdf'), p_volcano_zbp1_v_wt + xlim(c(-lfc_bound,lfc_bound)) + ylim(c(0,pval_bound)))

# Cleveland/Raincloud Plot Data
res_adar_v_wt_sig      = res_adar_v_wt %>% subset(abs(log2FoldChange) > DE_fc_cuttoff & padj < DE_pcuttoff)     %>% mutate(logP=-log10(padj))
res_pkr_v_wt_adar_sig  = res_pkr_v_wt  %>% subset(external_gene_name %in% res_adar_v_wt_sig$external_gene_name) %>% mutate(logP=-log10(padj))
res_zbp1_v_wt_adar_sig = res_zbp1_v_wt %>% subset(external_gene_name %in% res_adar_v_wt_sig$external_gene_name) %>% mutate(logP=-log10(padj))
#res_mda5_v_wt_adar_sig = res_mda5_v_wt %>% subset(external_gene_name %in% res_adar_v_wt_sig$external_gene_name) %>% mutate(logP=-log10(padj))
res_combined = bind_rows(list(ADAR=res_adar_v_wt_sig,ADAR_ZBP1_KO=res_zbp1_v_wt_adar_sig,ADAR_PKR_KO=res_pkr_v_wt_adar_sig ), .id='genotype_code') %>%  #, ADAR_MDA5_KO=res_mda5_v_wt_adar_sig), .id='genotype_code') %>% 
	mutate(genotype_code = factor(genotype_code,levels=c('ADAR_MDA5_KO','ADAR_ZBP1_KO','ADAR_PKR_KO','ADAR'))) %>% 
	mutate(adar_reg = ifelse(ensembl_gene_id_version %in% (res_adar_v_wt_sig %>% subset(log2FoldChange > 0) %>% pull(ensembl_gene_id_version)),'up','down')) %>% 
	mutate(adar_reg = factor(adar_reg,levels=c('up','down')))
make_delta_df = function(res_long,conditionA,conditionB){
	lfcA = paste0('log2FoldChange_',conditionA)
	lfcB = paste0('log2FoldChange_',conditionB)
	res_conditions      = res_long %>% subset(genotype_code %in% c(conditionA,conditionB))
	res_conditions_wide = res_conditions %>% 
		pivot_wider( id_cols=c(ensembl_gene_id_version,ensembl_gene_id,external_gene_name,entrezgene_id,adar_reg)
                 , names_from=genotype_code
                 , values_from=c(baseMean,log2FoldChange,lfcSE,stat,pvalue,padj))
	res_conditions_wide$delta_fc =  res_conditions_wide[[lfcA]] - res_conditions_wide[[lfcB]]
	res_conditions_wide$mean_fc  = (res_conditions_wide[[lfcA]] + res_conditions_wide[[lfcB]])/2
	res_conditions_wide = res_conditions_wide %>% 
		mutate(external_gene_name=factor(external_gene_name,levels=arrange(.,mean_fc) %>% pull(external_gene_name) %>% unique())) 
	#list(tall=res_conditions,wide=res_conditions_wide)
	res_conditions_wide
}
res_adar_v_zbp1 = make_delta_df(res_combined,'ADAR','ADAR_ZBP1_KO')
highdfc_gene_names = res_adar_v_zbp1 %>% subset(abs(delta_fc/log2FoldChange_ADAR) >  proportional_dfc_cut) %>% pull(external_gene_name)
res_adar_v_zbp1_high_dfc = res_adar_v_zbp1 %>% subset( external_gene_name %in% highdfc_gene_names) 
res_adar_v_zbp1_low_dfc  = res_adar_v_zbp1 %>% subset(!external_gene_name %in% highdfc_gene_names)
write_tsv(res_adar_v_zbp1_high_dfc,paste0(output_dir,'/res_adar_v_zbp1_high_dfc.tsv'))
write_tsv(res_adar_v_zbp1_low_dfc ,paste0(output_dir,'/res_adar_v_zbp1_low_dfc.tsv') )

# Raincloud Plots
raincloud = ggplot(res_combined %>% subset(genotype_code != 'ADAR_MDA5_KO'), aes(x=log2FoldChange,y=genotype_code)) +
	#stat_halfeye(alpha=.33, slab_type='pdf', adjust = .25, justification = -.3, height = .6, point_colour = NA, .width = 0) +
	stat_halfeye(aes(fill=adar_reg), alpha=.66, adjust = 1, justification = -.25, height = .66, point_colour = NA, .width = 0, normalize='xy') +
	geom_boxplot(aes(fill=adar_reg), alpha=0.33, width=.25, outlier.shape = NA) +
	geom_point(aes(colour=adar_reg, group=adar_reg), size=0.5, alpha=point_alpha, position=position_jitterdodge(seed=1, jitter.width=.1, dodge.width=.25) ) +
	xlab(TeX('$log_2$ Fold Change vs WT')) +
	ylab('') +
	scale_y_discrete(labels=c( 'ADAR'        =parse(text=TeX(study_group_code2latex['ADAR']        ))
                            , 'ADAR_ZBP1_KO'=parse(text=TeX(study_group_code2latex['ADAR_ZBP1_KO']))
                            , 'ADAR_PKR_KO' =parse(text=TeX(study_group_code2latex['ADAR_PKR_KO'] )) )) +
	scale_fill_manual  (values=tableau_color_pal()(10)[c(5,3)],name='Adar1 regulation') +
	scale_colour_manual(values=tableau_color_pal()(10)[c(5,3)],name='Adar1 regulation')
ggsave(paste0(output_dir,'/raincloud.pdf'), plot=raincloud, width=8.27, height=4 )

# Cleveland Plots
dummy_plot = ggplot(data.frame(list(x=1:2,y=1:2,type=c('ADAR','ZBP1_KO')))) +
	geom_point(aes(x=x,y=y,colour=type), size=cleveland_dot_size) +
	scale_colour_manual(values=c(adar_colour,zbp1_colour), name='Condition', labels=c( 'ADAR'        = parse(text=TeX(study_group_code2latex['ADAR']))
                                                                                    , 'ZBP1_KO'= parse(text=TeX(study_group_code2latex['ADAR_ZBP1_KO'])))) +
	theme(legend.text.align = 0)
cleveland_legend = get_legend(dummy_plot)

cleveland_high_relative_dfc = ggplot(res_adar_v_zbp1_high_dfc %>% slice_max(abs(log2FoldChange_ADAR),n=cleveland_gene_count)) +
	geom_segment(aes(x=external_gene_name,xend=external_gene_name,y=log2FoldChange_ADAR,yend=log2FoldChange_ADAR_ZBP1_KO), colour='grey') +
	geom_point(aes(x=external_gene_name,y=log2FoldChange_ADAR        ), alpha=point_alpha, colour=adar_colour, size=cleveland_dot_size) +
	geom_point(aes(x=external_gene_name,y=log2FoldChange_ADAR_ZBP1_KO), alpha=point_alpha, colour=zbp1_colour, size=cleveland_dot_size) +
	geom_hline(yintercept=0) +
	coord_flip() +
	xlab('') +
	ylab(TeX('$log_2$ Fold Change vs WT')) +
	labs(title='High Recovery')

ggsave(paste0(output_dir,'/cleveland_dotplot_highdFC.pdf'), width=7, height=21/4)

res_adar_v_zbp1_low_dfc_rev = res_adar_v_zbp1_low_dfc %>% 
	mutate(external_gene_name=factor(external_gene_name,levels=arrange(.,-mean_fc) %>% pull(external_gene_name) %>% unique()))
cleveland_low_relative_dfc = ggplot(res_adar_v_zbp1_low_dfc_rev %>% slice_max(abs(log2FoldChange_ADAR),n=cleveland_gene_count)) +
	geom_segment(aes(x=external_gene_name,xend=external_gene_name,y=log2FoldChange_ADAR,yend=log2FoldChange_ADAR_ZBP1_KO), colour='grey') +
	geom_point(aes(x=external_gene_name,y=log2FoldChange_ADAR        ), alpha=point_alpha, colour=adar_colour, size=cleveland_dot_size) +
	geom_point(aes(x=external_gene_name,y=log2FoldChange_ADAR_ZBP1_KO), alpha=point_alpha, colour=zbp1_colour, size=cleveland_dot_size) +
	geom_hline(yintercept=0) +
	coord_flip() +
	xlab('') +
	ylab(TeX('$log_2$ Fold Change vs WT')) +
	labs(title='Low Recovery')
ggsave(paste0(output_dir,'/cleveland_dotplot_lowdFC.pdf'), width=7, height=21/4)


# GO Analysis
adar_background_genes = res_adar_v_wt$ensembl_gene_id
getGO = function(x){
	ego = clusterProfiler::enrichGO( gene          = x$ensembl_gene_id
                                  , universe      = adar_background_genes
                                  , OrgDb         = org.Mm.eg.db::org.Mm.eg.db
                                  , keyType       = 'ENSEMBL'
                                  , ont           = "ALL"
                                  , pAdjustMethod = "BH"
                                  , qvalueCutoff  = 0.05
                                  , readable      = TRUE)
	ego %>% as.data.frame()
}
go_high_dfc = getGO(res_adar_v_zbp1_high_dfc)
go_low_dfc  = getGO(res_adar_v_zbp1_low_dfc)

go_combined = list(`High Recovery`=go_high_dfc, `Low Recovery`=go_low_dfc) %>% bind_rows(.id='dfc') %>% arrange(p.adjust)
top_go_terms = go_combined %>% slice_max(p.adjust,n=num_go_terms) %>% pull(ID) %>% unique()

go_results_shown = go_combined %>% 
	subset(ID %in% top_go_terms) %>% 
	separate('GeneRatio',into=c('gr_num','gr_den')) %>% 
	mutate(gr_num=as.integer(gr_num),gr_den=as.integer(gr_den),Count=as.integer(Count)) %>% 
	mutate(gene_ratio=gr_num/gr_den) %>% 
	group_by(ID) %>% 
	mutate(median_padj = median(p.adjust), median_gr = median(gene_ratio)) %>% 
	ungroup() %>% 
	arrange(desc(median_gr),median_padj) %>% 
	mutate(`GO Term`=paste0(capitalize(Description),' (',ONTOLOGY,')')) %>% 
	mutate(`GO Term`=factor(`GO Term`,levels=`GO Term` %>% unique() %>% rev()))

go_results_shown_high_dfc = go_high_dfc %>% 
	subset(ID %in% top_go_terms) %>% 
	separate('GeneRatio',into=c('gr_num','gr_den')) %>% 
	mutate(gr_num=as.integer(gr_num),gr_den=as.integer(gr_den),Count=as.integer(Count)) %>% 
	mutate(gene_ratio=gr_num/gr_den) %>% 
	arrange(desc(gene_ratio),p.adjust) %>% 
	mutate(`GO Term`=paste0(capitalize(Description),' (',ONTOLOGY,')')) %>% 
	mutate(`GO Term`=factor(`GO Term`,levels=`GO Term` %>% unique() %>% rev()))

go_results_shown_low_dfc = go_low_dfc %>% 
	subset(ID %in% top_go_terms) %>% 
	separate('GeneRatio',into=c('gr_num','gr_den')) %>% 
	mutate(gr_num=as.integer(gr_num),gr_den=as.integer(gr_den),Count=as.integer(Count)) %>% 
	mutate(gene_ratio=gr_num/gr_den) %>% 
	arrange(desc(gene_ratio),p.adjust) %>% 
	mutate(`GO Term`=paste0(capitalize(Description),' (',ONTOLOGY,')')) %>% 
	mutate(`GO Term`=factor(`GO Term`,levels=`GO Term` %>% unique() %>% rev()))

go_plot = ggplot(go_results_shown) +
	geom_point(aes(x=gene_ratio,y=`GO Term`,colour=p.adjust,size=Count)) +
	scale_y_discrete(labels= function(x) str_wrap(x,width=100)) +
	facet_grid(cols=vars(dfc)) +
	scale_colour_gradient_tableau(palette='Red', trans='reverse',name=parse(text=TeX('$P_{adj}$'))) +
	labs(title='GO Terms from DE vs WT') +
	xlab('Gene Ratio') +
	ylab('')

go_plot_hr = ggplot(go_results_shown_high_dfc) +
	geom_point(aes(x=gene_ratio,y=`GO Term`,colour=p.adjust,size=Count)) +
	scale_y_discrete(labels= function(x) str_wrap(x,width=100)) +
	scale_colour_gradient_tableau(palette='Red', trans='reverse',name=parse(text=TeX('$P_{adj}$'))) +
	xlab('Gene Ratio') +
	ylab('') +
	theme(text=element_text(size=16))

go_plot_lr = ggplot(go_results_shown_low_dfc) +
	geom_point(aes(x=gene_ratio,y=`GO Term`,colour=p.adjust,size=Count)) +
	scale_y_discrete(labels= function(x) str_wrap(x,width=100)) +
	scale_colour_gradient_tableau(palette='Red', trans='reverse',name=parse(text=TeX('$P_{adj}$'))) +
	xlab('Gene Ratio') +
	ylab('') +
	theme(text=element_text(size=16))

ggsave(paste0(output_dir,'/go.pdf')             , plot=go_plot   , width=12, height=32/100*num_go_terms)
ggsave(paste0(output_dir,'/go_highrecovery.pdf'), plot=go_plot_hr, width=12, height=32/100*num_go_terms)
ggsave(paste0(output_dir,'/go_lowrecovery.pdf') , plot=go_plot_lr, width=12, height=32/100*num_go_terms)

# Combind Grid
plot_grid(
	plot_grid(p_volcano_combined)
	, plot_grid(raincloud)
	, plot_grid(cleveland_high_relative_dfc, cleveland_low_relative_dfc, ncol=2) 
	, plot_grid(cleveland_legend)
	, plot_grid(go_plot)
	, align='h', nrow=5, rel_heights = c(1,1,1,0.33,1))
ggsave(paste0(output_dir,'/grid.pdf'), width=14, height=sqrt(2)*14)
