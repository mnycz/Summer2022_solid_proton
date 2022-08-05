import matplotlib
matplotlib.use('TkAgg')  # This line is sometimes required for plots to appear, but it slows down the program start a lot, so comment out if don't need to see plots
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
plt.close('all')

# Variables to change if changing BEAM_E
BEAM_E = '11'  # GeV
NQ2bins = '14'  # max number of Q2 bins per x (14 for 11 GeV and 30 for 22 GeV)
# How many rows at the bottom of each CSV file to ignore
# 11 GeV: [0, 0, 5]; 22 GeV [2, 1, 6]
footer_vals = [0, 0, 5]  # [2, 1, 6]  

# Variables to maybe change if changing BEAM_E

# Bin numbers
Nxbins_analytical = '10'
Nxbins_PDF = '100'

PDF_names = ["CT18NLO"]  # Source PDF

# Variables to store path information for the all bins
dir_path_al = "/w/eic-scshelf2104/users/gsevans/8thWeekSULIs22/" + BEAM_E + \
              "GeV_files/"
file_suffix_al = "_" + BEAM_E + "GeV_from" + NQ2bins + "x" + \
                 Nxbins_analytical + "grid_analytic_values_for_all_bins.csv"
foot_al = footer_vals[0]  # How many rows at the bottom to ignore (2 for 22 GeV, 0 for 11 GeV)

# Variables to store path information for the each x bins
dir_path_ea = "/w/eic-scshelf2104/users/gsevans/8thWeekSULIs22/" + BEAM_E + \
              "GeV_files/"
file_suffix_ea = "_" + BEAM_E + "GeV_from" + NQ2bins + "x" + \
                 Nxbins_analytical + "grid_analytic_values_for_each_x.csv"
foot_ea = footer_vals[1] # How many rows at the bottom to ignore (1 for 22 GeV, 0 for 11 GeV)

# Variables to store path information for the PDF values
dir_path_PDF = "/w/eic-scshelf2104/users/gsevans/8thWeekSULIs22/" + BEAM_E + \
               "GeV_files/"
file_suffix_PDF = "_" + BEAM_E + "GeV_qPDF_vals_1x" + Nxbins_PDF + "grid.csv"
foot_PDF = footer_vals[2]  # How many rows at the bottom to ignore(6 for 22 GeV, 5 for 11 GeV)

# Where to save the picture
dir_path_pic = "/w/eic-scshelf2104/users/gsevans/8thWeekSULIs22/" + BEAM_E + \
               "GeV_files/"
pic_type = ".pdf"  # What file type to save the picture as

def frame(ax, grid):
    ax.axhline(y=0, color='k',linewidth=2)
    ax.axhline(y=grid.shape[0], color='k',linewidth=1)
    ax.axvline(x=0, color='k',linewidth=1)
    ax.axvline(x=grid.shape[1], color='k',linewidth=2)

# Loop through all PDF_names creating plots(only tested w/ PDF_names=["CT18NLO"])
for pdf_name in PDF_names:
    # Read in data
    data_al = pd.read_csv(dir_path_al + pdf_name + file_suffix_al, \
                          skipfooter=foot_al)
    data_ea = pd.read_csv(dir_path_ea + pdf_name + file_suffix_ea, \
                          skipfooter=foot_ea)
    data_PDF = pd.read_csv(dir_path_PDF + pdf_name + file_suffix_PDF, \
                           skipfooter=foot_PDF)

    # Set up ticks and tick labels
    # Create arrays of the unique x values and y values
    uniquex = data_al['x'].unique()
    uniqueQ2 = data_al['Q2 (GeV^2)'].unique()
    # Set up x lables
    xstep = uniquex[1] - uniquex[0]  # Step size between x-values
    xlabels = np.append(uniquex, uniquex[-1] + xstep)  # Add last data point
    xlabels = xlabels - xstep/2  # Lines up the labels with the ticks
    xlabels = [np.format_float_positional(xlabel, precision=3, unique=True, 
        fractional=False, trim='-') for xlabel in xlabels]  # Format labels
    # Set up Q2 (y) labels
    Q2step = uniqueQ2[1] - uniqueQ2[0]
    Q2labels = np.append(uniqueQ2, uniqueQ2[-1] + Q2step)
    Q2labels = Q2labels - Q2step/2
    Q2labels = [np.format_float_positional(Q2label, precision=3, unique=True, fractional=False, trim='-') for Q2label in Q2labels]
    # Set tick placement, noting that sns has placement location depend on bins
    numx = len(uniquex)
    numQ2 = len(uniqueQ2)
    xticks = np.linspace(start=0, stop=numx, num=numx+1)
    Q2ticks = np.linspace(start=0, stop=numQ2, num=numQ2+1)
    print(xticks, Q2ticks)

    # A statistical heatmaps
    Aplot_width = 8.5
    Aplot_height = 5.5

    # A_noisy heatmap on Q2 vs x
    plt.figure(figsize=[Aplot_width, Aplot_height])
    AGrid = data_al[['x','Q2 (GeV^2)','A_noisy']].copy() * 1000000
    AGrid = AGrid.pivot(index='Q2 (GeV^2)', columns='x', values='A_noisy')
    ax = sns.heatmap(AGrid, cbar_kws={'label': '$A_{PV}$ (ppm)'}, annot=True,
                     cmap="RdYlBu_r") # https://stackoverflow.com/a/42092328
    ax.invert_yaxis()  # https://stackoverflow.com/a/34444939
    plt.ylabel('$Q^2$ (GeV$^2$)')
    plt.xlabel('$x$')
    ax.set_xticks(xticks)
    ax.set_yticks(Q2ticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(Q2labels)
    frame(ax, AGrid)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.savefig(dir_path_pic + "Anoisy_Q2xbin_heatmap_" + 
                BEAM_E + "GeV" + pic_type)
    
    # dA_stat/A_noisy heatmap on Q2 vs x
    plt.figure(figsize=[Aplot_width, Aplot_height])
    dAAGrid = data_al[['x','Q2 (GeV^2)']].copy()
    dAAGrid = pd.concat([dAAGrid, data_al['dA_stat'] / data_al['A_noisy']], axis=1)
    dAAGrid.rename(columns={0:'dA_stat/A_noisy'}, inplace=True)
    dAAGrid = dAAGrid.pivot(index='Q2 (GeV^2)', columns='x', 
                          values='dA_stat/A_noisy')
    ax = sns.heatmap(dAAGrid, cbar_kws={'label': '$\\delta A_{stat}/A_{PV}$'}, 
                     annot=True, cmap="RdYlBu_r", vmin=-0.5)
    ax.invert_yaxis()  
    plt.ylabel('$Q^2$ (GeV$^2$)')
    plt.xlabel('$x$')
    ax.set_xticks(xticks)
    ax.set_yticks(Q2ticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(Q2labels)
    frame(ax, dAAGrid)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.savefig(dir_path_pic + "dAstat_Anoisy_Q2xbin_heatmap_" + 
                BEAM_E + "GeV" + pic_type)
    
    # dA_stat heatmap on Q2 vs x
    plt.figure(figsize=[Aplot_width, Aplot_height])
    dAGrid = data_al[['x','Q2 (GeV^2)','dA_stat']].copy()
    dAGrid = dAGrid.pivot(index='Q2 (GeV^2)', columns='x', values='dA_stat')
    ax = sns.heatmap(dAGrid, cbar_kws={'label': '$\delta A_{stat}$'}, annot=True,
                     cmap="RdYlBu_r")
    ax.invert_yaxis()  
    plt.ylabel('$Q^2$ (GeV$^2$)')
    plt.xlabel('$x$')
    ax.set_xticks(xticks)
    ax.set_yticks(Q2ticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(Q2labels)
    frame(ax, dAAGrid)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.savefig(dir_path_pic + "dAstat_Q2xbin_heatmap_" + 
                BEAM_E + "GeV" + pic_type)
    
    # Rate heatmap on Q2 vs x
    plt.figure(figsize=[Aplot_width, Aplot_height])
    rateGrid = data_al[['x','Q2 (GeV^2)','rate (Hz)']].copy()
    rateGrid = rateGrid.pivot(index='Q2 (GeV^2)', columns='x', values='rate (Hz)')
    ax = sns.heatmap(rateGrid, cbar_kws={'label': 'rate (Hz)'}, annot=True,
                     cmap="RdYlBu_r")
    ax.invert_yaxis() 
    plt.ylabel('$Q^2$ (GeV$^2$)')
    plt.xlabel('$x$')
    ax.set_xticks(xticks)
    ax.set_yticks(Q2ticks)
    ax.set_xticklabels(xlabels)
    ax.set_yticklabels(Q2labels)
    frame(ax, dAAGrid)
    plt.grid(linestyle='--', linewidth=0.5)
    plt.savefig(dir_path_pic + "rate_Q2xbin_heatmap_" + 
                BEAM_E + "GeV" + pic_type)



    
    # 1D Comparison Plots
    
    # Create plot for uncertainty in A
    # Note, the displayed plot seems limited by your screen size, but the 
    #     saved pic doesn't. 13.5 is the max that fits in my screensize.
    plt.figure(figsize=[8.5,3])  # Poster size: [13.5, 3])
    plt.plot(abs(data_al['A_clean']), label="$A_{SM}^{Theory}$")
    plt.plot(data_al['dA_stat'], label="$\sigma_{stat}$")
    plt.plot(data_al['dA_sys_uncor'], label="$\sigma_{sys}$ $_{uncor}$")
    plt.plot(data_al['dA_sys_cor'], label="$\sigma_{sys}$ $_{cor}$")
    plt.xlabel('Bin Number #')
    plt.semilogy()
    plt.grid(linestyle='--', linewidth=0.5)
    # Position the legend, 11GeV:(0.45,0.48), 22GeV: (0.50,0.55)
    plt.legend(ncol=2, loc=(0.33,0.56))  # Poster placement: loc=(0.40,0.55))  
    plt.tight_layout(rect=[0.0, 0.01, 1, 1])
    plt.savefig(dir_path_pic + "ABin_uncertainty_comparison_plot_" + \
                BEAM_E + "GeV" + pic_type)
    
    # Create plot for uncertainty in d/u as a function of x
    plt.figure(figsize=[8.5,3])  # Poster size: [13.5, 3]) 
    x = data_ea['x']
    plt.plot(x, data_ea['d/u clean'], label="$(d_V/u_V)_{calc}$")
    plt.plot(x, data_ea['d(d/u) from A_stat'], label="$\sigma_{stat}$")
    plt.plot(x, data_ea['d(d/u) from A_sys_uncor'], 
             label="$\sigma_{sys}$ $_{uncor}$")
    plt.plot(x, data_ea['d(d/u) from A_sys_cor'], 
             label="$\sigma_{sys}$ $_{cor}$")
    plt.plot(x, data_ea['d(d/u) from sin^2(theta_w)'], 
             label="$\sigma_{\sin^2(\\theta_w)}$")
    plt.plot(x, data_ea['d(d/u) from model'], label="$\sigma_{model}$")
    plt.plot(data_PDF['x'], data_PDF['(dV/uV)_unc'], label="$\sigma_{PDF}$")
    plt.xlabel('$x$')
    plt.grid(linestyle='--', linewidth=0.5)
    plt.semilogy()
    plt.legend(ncol=4)
    plt.tight_layout(rect=[0, 0.0, 1, 1])
    plt.savefig(dir_path_pic + "dVuV_uncertainty_comparison_plot_" + \
                BEAM_E + "GeV" + pic_type)
    
plt.show()