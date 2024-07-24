"""
@author: Finn Rabe
@email: flrabe90@gmail.com
@terms: CC-BY-NC-ND
"""
# set working dir
import os
#os.chdir('/retpsy/src')
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# importing  all the functions defined in retpsy_load.py
from retinflam_load import *


def load_data_type(dis, disorder, format_type, meas_gr):
    """load data depending on different demands
        dis: define disorder abbreviation, e.g. SZ,BD,PD
        disorder: ex-/include diagnosed participants, choose 'diagnosed', 'none' or 'all'
        format_type: format oct data either to 'wide' or 'long' format (is needed for mixed linear model)
        meas_gr: select if you want to analyse general retinal phenotypes or macular subfields
    """
    # load and format data
    df_data = general_formatting(dis,disorder,format_type,meas_gr)
    y_meas,cov,cp = sel_oct_group(meas_gr, format_type)
    x_meas = ['PRS'+dis]
    df_data_trans, y_meas_lat = oct_formatting(df_data,meas_gr,y_meas,format_type,disorder)
    #df_data_trans.iloc[0].to_csv('/Users/frabe/Desktop/col.csv')
    return df_data_trans,y_meas,y_meas_lat,x_meas,cov,cp

def parreg_pc(df_zcorr,pc_lbl,x_meas,figname='Fig3'):
    # Partial regression with first principle component 
    figlbl = ['a','b']
    ax_dim = [1,1]
    pc_len = 1

    count = 0 
    for yv in range(pc_len):
        fig, subfig = plt.subplots(ax_dim[0], ax_dim[1], sharex=True, figsize=fig_size)
        fig.tight_layout() 
        # Regress out nuissance regressors 
        #xarr = regress_out(df_corr[pc_lbl[yv]], df_corr[cov])
        xarr = regress_out(df_zcorr[x_meas[0]], df_zcorr[cov])

        # Compute partial regression incl. covariates as nuissance reg
        pstats = pg.partial_corr(data=df_zcorr, y=pc_lbl[yv], x=x_meas[0], covar=cov, method='pearson')
        r_partial = np.round(pstats.r.values[0],2)
        p_partial = np.round(pstats['p-val'].values[0],3)
        df_results[pc_lbl[yv]+'partialr'] = [r_partial]
        df_results[pc_lbl[yv]+'partialp'] = [p_partial]

        axid = subfig

        g = sns.regplot(ax=axid, y=xarr, x=df_zcorr[pc_lbl[yv]], scatter_kws={"color": "white"}, line_kws={'linewidth':line_size}, color=cp[yv], ci=95, robust=True)
        #add regression equation to plot
        #axx = g.axes[xv,yv]
        #g.text(0, 0.1,'r\u00b2= {}'.format(r_partial) +pval_asterisks(p_partial,0.05), fontsize=val_size)
        # g.text(0, 0.1,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
        # g.set(ylim=(-0.3,0.3))
        # #g.xaxis.set_ticks(g.get_xticks())
        # #g.yaxis.set_ticks(ylabels)
        # xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
        # print('xticks:',xlabels)
        # g.set_xticklabels(xlabels, fontsize=tick_size)
        # ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
        # g.set_yticklabels(ylabels,fontsize=tick_size)

        
        if figname == 'Fig1':
            xlbl = 'PRS SZ (z)'
            g.text(0, 0.05,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.1,0.1))
            g.set(xlim=(-20,20))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        elif figname == 'Fig3A':
            xlbl = 'Polygenic risk score for neuroinflammatory pathway'
            #g.text(0, 0.05,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.1,0.1))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            ylbl = 'C-reactive protein log(mg/L)'
            count = 1
            g.set_ylabel(ylbl,fontsize=10)
            g.set_xlabel(xlbl,fontsize=10)
        elif figname == 'Fig3B':
            xlbl = 'Polygenic risk score for WnT signaling pathway'
            g.text(0, 0.05,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.1,0.1))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        elif figname == 'Fig3C':
            xlbl = 'C-reactive protein log(mg/L)'
            g.text(0, 0.4,'r = {}'.format(r_partial) +' | '+ pval_asterisks(p_partial, alpha=0.05), fontsize=val_size)
            g.set(ylim=(-0.2,0.7))
            xlabels = ['{:,.2f}'.format(x) for x in g.get_xticks()]
            g.set_xticklabels(xlabels, fontsize=tick_size)
            ylabels = ['{:,.2f}'.format(x) for x in g.get_yticks()]
            g.set_yticklabels(ylabels,fontsize=tick_size)
            g.set_xlabel(xlbl,fontsize=label_size)
        if count == 0:
            ylbl = 'Principal component 1'
            g.set_ylabel(ylbl,fontsize=label_size)
            #g.set_ylabel(g.get_ylabel(),fontsize=label_size)
        #else:
        #    g.set_ylabel('',fontsize=label_size)
        count += 1
        sns.despine()

        plt.savefig('../output/figures/'+figname+figlbl[yv]+'_'+meas_gr+'_'+pc_lbl[yv]+'_partialreg.png', format='png', bbox_inches='tight')

def pca_loadscores(df_corr, y_meas, pc_lbl, pca, figname = 'Fig2B'):
    # weights are nothing but the eigenvectors of X
    y_df = df_corr[y_meas]
    fig, g = plt.subplots(figsize=fig_size)

    #cp = sns.color_palette("colorblind")
    cp = sns.color_palette("Paired")
    pca_weights = pca.components_.flatten()

    df_pca_weights = pd.DataFrame(data={'Principal Component': np.repeat(pc_lbl,len(y_df.columns)),\
                        'Variables': np.tile(y_df.columns,len(pc_lbl)),\
                        'PCA loading scores': pca_weights})

    for lc in range(len(pca_weights)):
        pcn = df_pca_weights['Principal Component'][lc]
        retp = df_pca_weights['Variables'][lc]
        ldc = df_pca_weights['PCA loading scores'][lc]
        df_results[meas_gr+'_pca'+pcn+'_'+retp+'_loadingscore'] = [np.round(ldc,2)]
    #sns.set(rc={'figure.figsize':(10,8)})
    labels = ['II left', 'II right', 'IO left', 'IO right',
                         'NI left','NI right','NO left', 'NO right', 
                         'SI left', 'SI right','SO left','SO right',
                         'TI left', 'TI right','TO left','TO right',
                         'CS left', 'CS right',
                         ]
    abbr_lbl = np.tile(labels,2)
    df_pca_weights['legend_lbl'] = abbr_lbl
    g = sns.barplot(x=df_pca_weights.columns[0], y=df_pca_weights.columns[2], hue=df_pca_weights.legend_lbl, data=df_pca_weights, palette=cp)
    #axes.legend(loc='lower left', bbox_to_anchor=(1, 0.5),fontsize = tick_size)
    if meas_gr == 'overall':
        ti = 'Retinal phenotypes'
    else:
        ti = 'Macular bilateral subfields'
    sns.move_legend(g,'lower center', bbox_to_anchor=(.5, 1), ncol=3, fontsize = tick_size, title=ti, frameon=False)
    sns.despine()
    fig = g.get_figure()
    fig.savefig('../output/figures/'+figname+'_'+meas_gr+'_pca_loadingscores.png', format='png', bbox_inches='tight')
    return df_pca_weights

def pca_cumvar(df_corr,y_meas,var_cutoff=0.8, figname = 'Fig2A'):
    # Compute Principle component analysis
    
    print('df cols:',df_corr.columns)
    print(y_meas)
    y_df = df_corr[y_meas]
    pca = PCA(n_components=var_cutoff)
    pcs = pca.fit_transform(y_df)

    pc_lbl = ["PC%01d" %i for i in range(1,pca.n_components_+1)]
    pca = PCA(n_components=len(pc_lbl))

    # weight oct measures
    pcs = pca.fit_transform(y_df)

    # cumulative variance 
    exp_var = pca.explained_variance_ratio_
    cum_sum_eigenvalues = np.cumsum(pca.explained_variance_ratio_)
    print('shape before pcs:',df_corr.shape)
    for pc in range(len(pc_lbl)):
        df_corr[pc_lbl[pc]] = pcs[:,pc]
    
    # save to results df
    if figname == 'Fig2A':
        df_results[meas_gr+'_pca_ncomponents'] = [pca.n_components_]
        for ev in range(len(exp_var)):
            print('varr:',np.round(exp_var[ev],2))
            df_results[meas_gr+'_pca'+str(ev+1)+'_expvariance'] = [np.round(exp_var[ev],2)]
    
    # Display Scree plot
    df_pca = pd.DataFrame(data={'Principal Component': pc_lbl,\
                        'Variance explained': pca.explained_variance_ratio_,\
                        'Cumulative Variance': cum_sum_eigenvalues
                        })
    plt.figure(figsize=fig_size)
    #sns.barplot(x=df_pca.columns[0], y=df_pca.columns[1], data=df_pca, palette=cp[:len(pc_lbl)])
    #plt.step(range(0,len(cum_sum_eigenvalues)), cum_sum_eigenvalues, where='mid',label='Cumulative explained variance')

    sns.barplot(x=df_pca.columns[0], y=df_pca.columns[2], data=df_pca, palette=cp[:len(pc_lbl)])
    sns.lineplot(x=df_pca.columns[0], y=df_pca.columns[2], data=df_pca, marker='o', sort = False)
    plt.axhline(var_cutoff, ls='--',color=thr_line)
    #plt.text(0,0.9, '90% threshold')
    plt.ylim((0,1))
    sns.despine()
    print('shape after pcs:',df_corr.shape)

    plt.savefig('../output/figures/'+figname+'_'+meas_gr+'_pca.png', format='png',bbox_inches='tight')
    return pc_lbl, pca, df_pca, df_corr, exp_var, cum_sum_eigenvalues

def corr_matrix(df_corr,y_meas,figname='Appx_Fig2'):
    covm_all = df_corr[y_meas].corr()
    ma = covm_all.round(2)

    fig, g = plt.subplots(figsize=fig_size)
    cp = sns.color_palette("YlOrRd", as_cmap=True)
    g = sns.heatmap(ma, annot=False,cmap=cp)
    #g = sns.heatmap(ma, annot=True,annot_kws={"fontsize":val_size},cmap=cp)

    yl = [t.get_text()  for t in g.get_yticklabels()]
    y_lbl = [s.strip('_') for s in yl]
    #y_lbl = ['GC-IPL left','GC-IPL right','Macula left','Macula right','RNFL left','RNFL right','INL left','INL right']#,'RPE left','RPE right']

    g.set_xticklabels(y_lbl, fontsize=tick_size)
    g.set_yticklabels(y_lbl, fontsize=tick_size)

    fig = g.get_figure()
    fig.savefig('../output/figures/'+figname+'_'+meas_gr+'_corrmatrix.png', format='png', bbox_inches='tight')

def population(raw,filt_hc,filt_scz,figname='Appx_Table1'):
    # number of participants
    n_total = len(raw)
    n_filt = len(filt_hc)
    n_scz = len(filt_scz)

    # sex ratio
    sex_vcounts = filt_hc.Sex.value_counts()

    # Compute female, male ratio
    df_sex = pd.DataFrame()
    df_sex['Sex'] = ['Male','Female']
    df_sex['N'] = [sex_vcounts[1],sex_vcounts[2]]

    t,p = stats.ttest_ind(filt_hc.Sex==1, filt_hc.Sex==2)

    fig =  ff.create_table(df_sex)
    fig.update_layout(
        autosize=False,
        width=700,
        height=400,
    )
    fig.write_image('../output/figures/'+figname+'_'+meas_gr+'_sex_ratio.png', scale=2)
    return n_total, n_filt, n_scz, sex_vcounts, t, p

def collinarity_dist(df_corr, y_meas, figname = 'Appx_Fig1'):
    ## Compute significance of collinarity of OCT measures
    df_coll = df_corr[y_meas].rcorr(padjust = 'bonf')

    #replace nan values with empty string
    coll_val = df_coll.values.flatten()
    coll_val[coll_val == '-'] = ''
    #coll_val = coll_val + ' $^3$'

    # hide upper triangle 
    def hide_current_axis(*args, **kwds):
        plt.gca().set_visible(False)

    fig, g = plt.subplots(figsize=fig_size)

    g = sns.pairplot(df_corr[y_meas], kind="reg", corner=False,\
                        #palette=scatter_plots,
                        markers="o", \
                        diag_kws = {'color': scatter_plots},\
                        plot_kws={'line_kws':{'color':reg_line}, 'scatter_kws': {'alpha': 0.5, 'color': scatter_plots}})
    
    # iterate through each axes and insert collinerity value
    count = 0
    for xpl in range(len(y_meas)):
        for ypl in range(len(y_meas)):
            #add regression equation to plot
            axx = g.axes[xpl,ypl]

            if meas_gr == 'overall':
                #axx.set(ylim=(-20, 50))
                #axx.set(xlim=(0, 800))
                if coll_val[count] == '':
                    axx.text(15, 20, str(coll_val[count]), size=16, color='black')
                else:
                    axx.text(15, 20, str(coll_val[count]), size=16, color='black') #+r'$^{***}$'
            else:
                #axx.set(ylim=(0, 1000))
                #axx.set(xlim=(0, 1000))
                if coll_val[count] == '':
                    axx.text(0, 10, str(coll_val[count]), size=23, color='black')
                else:
                    #axx.text(0, 10, str(coll_val[count])+r'$^{***}$', size=16, color='black')
                    axx.text(0, 10, str(coll_val[count]), size=23, color='black')
            count += 1

    g.map_upper(hide_current_axis)
    plt.savefig('../output/figures/'+figname+'_'+meas_gr+'_regmatrix.png', format='png')

def vif(df_corr,y_meas,figname = 'Appx_Table2'):
    ## Compute variance inflation factor
    # VIF dataframe
    df_sex = pd.DataFrame()
    df_sex['Phenotype'] = y_meas
    
    # calculating VIF for each feature
    df_sex['VIF'] = [variance_inflation_factor(df_corr[y_meas].values, i)
                            for i in range(len(y_meas))]
    df_sex['VIF'] = df_sex['VIF'].round(2)
    df_sex.to_csv('../output/data/'+figname+'_'+meas_gr+'_vif.csv')
    
    fig =  ff.create_table(df_sex)
    fig.update_layout(
        autosize=False,
        width=700,
        height=400,
    )
    fig.write_image('../output/figures/'+figname+'_'+meas_gr+'_vif.png', scale=2)

def save_pcs(df_dlm_path,df_res,y_meas_lat):
    df_dlm = pd.read_csv(df_dlm_path+'.csv')
    #y_meas_lat = ['GC_IPL_left','GC_IPL_right','Macula_left','Macula_right','RNFL_left','RNFL_right']

    y_df = df_res[y_meas_lat]
    pca = PCA(n_components=0.7)
    pcs = pca.fit_transform(y_df)
    pc_lbl = ["PC%01d" %i for i in range(1,pca.n_components_+1)]
    for pc in range(len(pc_lbl)):
            df_dlm[pc_lbl[pc]] = pcs[:,pc]
    df_dlm.to_csv(df_dlm_path+'_pc.csv')
    return df_dlm


if __name__ == "__main__":
    # set figure layout
    fig_size, tick_size, label_size, title_size, val_size, line_size, dist_plots, scatter_plots, reg_line, thr_line = set_fig_style()
    
    #eye indicies 
    side = ['left','right'] 

    ## Partial regression of PC1 and PRS
    meas_gr = 'subfields'
    df_data_trans, y_meas,y_meas_lat,x_meas,cov,cp = load_data_type('SZ', 'none', 'wide', meas_gr)
    col_incl = y_meas_lat+x_meas+cov
    df_data_zcorr = drop_nan(df_data_trans,col_incl)
    print('loaded dfs do:',df_data_zcorr[col_incl].columns)
    #df_data_zcorr = zscore_data(df_data_corr)

    # Covariance matrix of all subfields incl. left and right
    corr_matrix(df_data_zcorr,y_meas_lat,'Fig1')

    # Variance inflation factors
    vif(df_data_zcorr,y_meas_lat,figname = 'Appx_Table2')
    print('VIF done')

    # Collinarity between and distribution of retinal phenotypes
    collinarity_dist(df_data_zcorr,y_meas_lat,'Appx_Fig1')
    print('coll plot done!')
    
    # # Compute PCA and cumulative variance explained
    # load residuals from individual rlm
    df_res = pd.read_csv('../output/data/retinflam_rlm_resid.csv')
    print('col included for pc:', y_meas_lat)
    df_res_zcorr = zscore_data(df_res[y_meas_lat])
    pc_lbl, pca, df_pca, df_corr, exp_var, cum_sum_eigenvalues = pca_cumvar(df_res,y_meas_lat,0.7,'Appx_Fig2A')
    print('PCA plot done!')

    df_results = pd.read_csv('../output/data/retinflam_results_m.csv')
    df_results['exp_var'] = [np.round(exp_var,2)]
    df_results['cum_sum_eigenvalues'] = [np.round(cum_sum_eigenvalues[-1],2)]
    ### Save dataframe with variables ###
    df_results.to_csv('../output/data/retinflam_results_m.csv')
    print('Done Wide!!!')

    # Fit the model with retinal phenotypes and apply the dimensionality reduction on retinal phenotypes and add to df
    df_pc = save_pcs('../output/data/df_mlm',df_res,y_meas_lat)

    # # Derive loading scores from PCA results
    df_pca_weights = pca_loadscores(df_res,y_meas_lat,pc_lbl, pca,'Appx_Fig2B')
    print('PCA weight plot done!')

    # Partial regression of PC1 and GOBP pathway
    xvar = 'NEUROINFLAM_PRS'
    cov = ['Age', 'Age_squared', 'Sex', 'Fasting_time', 'Smoking_status', 'BMI', 
           'Genotype_array', 'Townsend_index', 
           'Genetic PC1', 'Genetic PC2', 'Genetic PC3', 'Genetic PC4', 'Genetic PC5', 
           'Genetic PC6', 'Genetic PC7', 'Genetic PC8', 'Genetic PC9', 'Genetic PC10']
    #cov = list(filter(lambda x: x!=xvar, cov)) #exclude CRP
    cov = cov #+ x_meas
    parreg_pc(df_data_zcorr,['CRP'],[xvar],'Fig3A')

