import os
import gseapy

databases = ('huMAP',
             'Tissue_Protein_Expression_from_Human_Proteome_Map',
             'Tissue_Protein_Expression_from_ProteomicsDB',
             'Transcription_Factor_PPIs',
             'Rare_Diseases_AutoRIF_Gene_Lists',
             'Rare_Diseases_GeneRIF_Gene_Lists',
             'PheWeb_2019',
             'MSigDB_Computational',
             'MGI_Mammalian_Phenotype_Level_4_2019',
             'Ligand_Perturbations_from_GEO_down',
             'Ligand_Perturbations_from_GEO_up',
             'Human_Phenotype_Ontology',
             'GO_Molecular_Function_2018',
             'GO_Cellular_Component_2018',
             'GO_Biological_Process_2018',
             'DSigDB',
             'WikiPathways_2019_Human',
             'WikiPathways_2019_Mouse',
             'KEGG_2019_Human',
             'KEGG_2019_Mouse',
             'GO_Molecular_Function_2018',
             'GO_Cellular_Component_2018',
             'GO_Biological_Process_2018',
             )




class GSEAnalysis:
    """
    A class to manage categorical and single sample gene set enrichment analysis (GSEA). The GSEAnalysis class manages file input and output as well as running the analysis and handling results
    """
    def __init__(self, base_out_directory, extension='/gsea', hallmarks_db_file=None):
        """Create a new GSEAnalysis instance which stores file paths loading and saving data

        Args:
            base_out_directory (str): base directory for gsea results to be stored.
            extension (str, optional): the folder in which files produced by analysis will be stored. Defaults to '/gsea'.
            hallmarks_db_file (_type_, optional): Path to a specific hallmarks file to load for gsea. If left to default, automatically loads msigdb hallmarks dataset. Defaults to None.
        """
        # if hallmark database file is not provided, default to the gsea databases directory and standard hallmarks db
        hallmarks_db = hallmarks_db_file if hallmarks_db_file is not None else os.path.join('gsea_databases', 'msigdb_hallmarks.gmt')
        gsea_methods = ('signal_to_noise', 't_test', 'ratio_of_classes',
                'diff_of_classes', 'log2_ratio_of_classes')
        ssgsea_methods = ('rank', 'log', 'log_rank')
        self.base_out = base_out_directory
        self.gsea_save = self.base_out + extension

        if not os.path.isdir(self.gsea_save):
            os.makedirs(self.gsea_save)

    def run_gsea(self, ordered_df, classes, db=None, db_name=None,
                 processes=4, no_plot=True, method='signal_to_noise',
                 permutation_type='gene_set'):
        """
        Run categorical GSEA using the gseapy package.

        Args:
            ordered_df (pandas.DataFrame): a df with samples ordered by class (experiment and control).
            classes (_type_): the names of each class
            db (str, optional): the file name of the path . Defaults to None.
            db_name (str, optional): name of the save directory. Defaults to None.
            processes (int, optional): number of processes (for multiprocessing). Defaults to 4.
            no_plot (bool, optional): produce plots via gseapy. Defaults to True.
            method (str, optional): the evaluation method. Defaults to 'signal_to_noise'.
            permutation_type (str, optional): can be gene_set or categorical. Defaults to 'gene_set'.

        Returns:
            _type_: _description_
        """
        if db is None:
            db = self.hallmarks_db
        if db_name is None:
            db_name = db.split("/")[-1]
        gsea_analysis = gseapy.gsea(ordered_df,
                                    gene_sets=db,
                                    cls=classes,
                                    outdir=os.path.join(self.gsea_save, db_name),
                                    method=method,
                                    no_plot=no_plot,
                                    processes=processes,
                                    permutation_type=permutation_type)
        return gsea_analysis

    def run_ssGSEA(self, ordered_df, db,db_name=None,
                   method='rank', no_plot=True,
                   processes=4):
        """
        Run single sample GSEA using the gseapy package.

        Args:
            ordered_df (_type_): a df with samples ordered by class (experiment and control).
            db (_type_): the file name of the path . Defaults to None.
            db_name (_type_, optional): name of the save directory. Defaults to None.
            method (str, optional): the evaluation method. Defaults to 'rank'.
            no_plot (bool, optional): produce plots via gseapy.. Defaults to True.
            processes (int, optional): number of processes (for multiprocessing). Defaults to 4.

        Returns:
            _type_: _description_
        """
        if db is None:
            db = self.hallmarks_db
        if db_name is None:
            db_name = db.split("/")[-1]
        ssgsea_analysis = gseapy.ssgsea(ordered_df,
                                        gene_sets=db,
                                        outdir='ss_'+self.gsea_save,
                                        sample_norm_method=method,
                                        no_plot=no_plot,
                                        processes=processes)
        return ssgsea_analysis

    @classmethod
    def construct_classes(cls, disease_name, healthy_name,
                           disease_num, healthy_num):
        return [*[disease_name]*disease_num, *[healthy_name]*healthy_num]

