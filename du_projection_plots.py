import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.close('all')


# Variables to change if changing BEAM_E
BEAM_E = '22'  # GeV
NQ2bins = '30'  # (Max) number of Q2 bins (30 if 22 GeV, 14 if 11 GeV)
# How many rows at the bottom of each CSV file to ignore
# 11 GeV: 0; 22 GeV: 1
footer_vals = 1


# Variables to maybe change if changing BEAM_E

Nxbins_analytical = '10'
Nxbins_PDF = '100'
stacked_errors = False  # Format of band errors

# Information for loading uncertainty data
PDF_names = ["CT18NLO"]
# Variables to store path information for analytic values
dir_path_a = "./Files/" + BEAM_E + "GeV_files/"
file_suffix_a = "_" + BEAM_E + "GeV_from" + NQ2bins + "x" + Nxbins_analytical + "grid_analytic_values_for_each_x.csv"
foot_a = footer_vals  # How many bottom rows to ignore. 1 for 22 GeV, 0 for 11 GeV
# Variables to store path information for fitted values
dir_path_f = "./Files/" + BEAM_E + "GeV_files/"
file_suffix_f = "_" + BEAM_E + "GeV_from" + NQ2bins + "x" + Nxbins_analytical + "grid_fitted_values_for_each_x.csv"
foot_f = footer_vals # How many bottom rows to ignore. 1 for 22 GeV, 0 for 11 GeV

# Information for loading the PDF data
PDF_names2 = ["NNPDF40_nlo_as_01180", "PDF4LHC21_40", "CT18NLO"]
dir_path_PDF = "./Files/" + BEAM_E + "GeV_files/"
file_suffix_PDF = '_' + BEAM_E + 'GeV_qPDF_vals_1x' + Nxbins_PDF + 'grid.csv'
foot_PDF = 0

# Where to save the picture
dir_path_pic = "./Files/" + BEAM_E + "GeV_files/"
pic_suffix = "main_dVuV_uncertainty_plot" + BEAM_E + "GeV.pdf"

# Begin program

# Plot stuff
alp = .4  # Transparency in bands
fig, ax = plt.subplots(figsize=[6,5])  # Report size: [6,5]  # Poster Size: [11,8])

# loop through all PDF names (only tested with PDF_names=["CT18NLO"])
for pdf_name in PDF_names:
    # Data from analytically obtaining and combining d/u values from pseudo-data
    data_a = pd.read_csv(dir_path_a + pdf_name + file_suffix_a, 
                         skipfooter=foot_a)
    
    # Data from obtaining d/u from fitting pseudo-data (not currently used)
    data_f = pd.read_csv(dir_path_f + pdf_name + file_suffix_f, 
                         skipfooter=foot_f)

    plot_list = []  # List ofplots for first legend
    x = np.array(data_f.loc[:,"x"])  # x (Bjorken scaling variable)
    du = -.8*x + .95  #np.array(data_f.loc[:,"d/u"])  # d/u

    # Make the points
    plot_list.append(ax.scatter(x, du, marker="."))
    # Errorbars in d/u due to statistical uncertainty in A
    plot_list.append(ax.errorbar(x, du, fmt="none", linewidth=1, capsize=5, 
                                 label="$\delta(d_V/u_V)$ from $A_{stat}$", 
                                 yerr=data_a.loc[:,"d(d/u) from A_stat"]))
    # Errorbars in d/u due to total uncorrelated uncertainty (both stat and sys)
    plot_list.append(ax.errorbar(x, du, fmt="none", linewidth=1, capsize=3, 
        label="$\delta(d_V/u_V)$ from $A_{uncor}$ ($\delta A_{stat}$ with" + \
              " $\delta A_{uncor}$ $_{sys})$",
        yerr=np.sqrt((data_a.loc[:,"d(d/u) from A_stat"]) ** 2 + 
             (data_a.loc[:,"d(d/u) from A_sys_uncor"]) ** 2)))

    # Declare the errors that will be plotted as bands. Note that labels will
    # also need to be updated if these are changed.
    bottom_err =  np.array(data_a.loc[:,"d(d/u) from model"])
    middle_err =  np.array(data_a.loc[:,"d(d/u) from sin^2(theta_w)"])
    top_err =  np.array(data_a.loc[:,"d(d/u) from A_sys_cor"])

    if stacked_errors:
        bot_line = -1
        # Errorbars in d/u due to correlated systematic uncertanty
        plot_list.append(
            ax.fill_between(x, top_err + middle_err + bottom_err + bot_line, 
                            middle_err + bottom_err + bot_line,
                            alpha = alp, label="$\delta(d_V/u_V)$ from " + \
                                               "$\delta A_{sys}$ $_{cor}$"))
        # Errorbars in d/u due to uncertainty in sin^2(theta_w)
        plot_list.append(
            ax.fill_between(x, middle_err + bottom_err + bot_line, 
                            bottom_err + bot_line,
                            alpha = alp, label="$\delta(d_V/u_V)$ from " + \
                                               "$\sin^2(\\theta_w) \pm 0.0006$"))
        # Errorbars in d/u due to ignoring sea, s, and c quarks
        plot_list.append(
            ax.fill_between(x, bottom_err + bot_line, bot_line,
                            alpha = alp, label="$\delta(d_V/u_V)$ from model"))
    else:
        bot_line = -.9
        # Errorbars in d/u due to correlated systematic uncertanty
        plot_list.append(
            ax.fill_between(x, .2 + top_err + bot_line, .2 + bot_line, 
                            alpha = alp, label="$\delta(d_V/u_V)$ from " + \
                                               "$\delta A_{sys}$ $_{cor}$"))
        # Errorbars in d/u due to uncertainty in sin^2(theta_w)
        plot_list.append(
            ax.fill_between(x, .1 + middle_err + bot_line, .1 + bot_line,
                            alpha = alp, label="$\delta(d_V/u_V)$ from " + \
                                               "$\sin^2(\\theta_w) \pm 0.0006$"))
        # Errorbars in d/u due to ignoring sea, s, and c quarks
        plot_list.append(
            ax.fill_between(x, bottom_err + bot_line, bot_line,
                            alpha = alp, label="$\delta(d_V/u_V)$ from model"))

# Create a legend for the first set of plots.
first_legend = plt.legend(handles=plot_list, loc=(0.02,.17))
ax.add_artist(first_legend)

# Plot PDF values with their uncertainties
plot_list2 = []
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']  # List of colors
i = 5  # Starting color
for pdf_name in PDF_names2:
    data_PDF = pd.read_csv(dir_path_PDF + pdf_name + file_suffix_PDF, skipfooter=foot_PDF)
    plot_list2.append(
        ax.fill_between(data_PDF["x"], 
                        data_PDF["dV/uV"] + data_PDF["(dV/uV)_unc"], 
                        data_PDF["dV/uV"] - data_PDF["(dV/uV)_unc"],
                        alpha = alp, label=pdf_name, color=cycle[i]))
    ax.plot(data_PDF["x"], data_PDF["dV/uV"] + data_PDF["(dV/uV)_unc"], color=cycle[i], linewidth=0.2)
    ax.plot(data_PDF["x"], data_PDF["dV/uV"] - data_PDF["(dV/uV)_unc"], color=cycle[i], linewidth=0.2)
    i += 1

# Plot extra info for help in understanding graph
ax.plot([0,1], [0,0], 'k', linewidth=.5)
#ax.text(.15, -.1, "Beam Polarization: " + str(.85))
#ax.text(.15, .05, "Runtime: " + str(90) + " days")

# Plot labels and formatting
fig.tight_layout(rect=[0.05, 0.02, 1, 1])  
plt.xlim(0,1)
plt.ylim(-1,1)
plt.xlabel("$x$")
plt.ylabel("$d_V/u_V$")
plt.legend(handles=plot_list2, loc=(0.4, 0.82))

plt.savefig(dir_path_pic + pic_suffix)
plt.show()
