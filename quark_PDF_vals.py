import matplotlib
matplotlib.use('TkAgg')  # Sometimes required for plotting to work
import sys
import numpy as np
import warnings
from scipy import stats
warnings.filterwarnings('ignore')
import collections
import pandas as pd
import matplotlib.cm as cm
from matplotlib import pyplot as plt
plt.close('all')

# Set up LHAPDF so PDF values can be obtained
File_Path_LHAPDF='/work/eic/users/xiaochao/lhapdf640/lib/python2.7/site-packages/'  # Where the JAM and other PDF / SF files are
sys.path.append(File_Path_LHAPDF) 
import lhapdf
lhapdf.setVerbosity(0)  # stops lhapdf from printing whenever called

def calc_sfn(x, Q2, tabname, iset, eigenset):
    """ Obtains the PDF value for the given quark (iset) at x,Q2 using 
        eigenset from the PDF set (tabname) """
    set0     = lhapdf.mkPDF(tabname, eigenset)
    sfn      = np.array(set0.xfxQ2(iset, x, Q2))
    return (sfn)

def get_q_PDF_vals(row, Grid, eigenset):
    """ Calculate quark PDF values for the given row of the given eigenset of 
        the given grid """
    d = (calc_sfn(row['x'], row['Q2'], Grid, 1, eigenset))
    d_bar = (calc_sfn(row['x'], row['Q2'], Grid, -1, eigenset))
    u = (calc_sfn(row['x'], row['Q2'], Grid, 2, eigenset))
    u_bar = (calc_sfn(row['x'], row['Q2'], Grid, -2, eigenset))
    
    s = (calc_sfn(row['x'], row['Q2'], Grid, 3, eigenset))
    s_bar = (calc_sfn(row['x'], row['Q2'], Grid, -3, eigenset))
    c = (calc_sfn(row['x'], row['Q2'], Grid, 4, eigenset))
    c_bar = (calc_sfn(row['x'], row['Q2'], Grid, -4, eigenset))

    return [d, d_bar, u, u_bar, s, s_bar, c, c_bar]

def calc_unc(row, Grid, N_PDF, base):
    """ Calculate the uncertainty in the given row of the given grid """
    dO = 0
    base = base.lower()
    # Calculate uncertainty depending on what base the PDF set is
    if base == "hessian":
        if N_PDF%2:
            for k in range(1,  int(N_PDF/2) + 1):
                dO += (np.array(get_q_PDF_vals(row, Grid, 2*k)) - 
                       np.array(get_q_PDF_vals(row, Grid, 2*k - 1))) ** 2
            dO = np.sqrt(dO / 4)
        else:
            print("The equation for the PDF uncertainty of Hessian-based " + \
                  "sets does not work properly if N_PDF is an even number. " + \
                  "Uncertainty set to NaN.")
            dO = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]
    elif base == "replica":
        O_0 = np.array(get_q_PDF_vals(row, Grid, 0))
        for m in range(1, N_PDF):
            dO += (np.array(get_q_PDF_vals(row, Grid, m)) - O_0) ** 2
            #print("TESTING", N_PDF, m)
        dO = np.sqrt(dO / N_PDF)
    else:
        print("The submitted base for", Grid, "was not recognized.", \
              "Uncertainty set to NaN.")
        dO = [np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN, np.NaN]
    return dO

### USER INPUT ###

# Handle file information
beamE = "22"  # GeV, Beam energy
NQ2bins = "30"  # (max) number of Q2 bins
Nxbins = "100"  # (max) number of x bins
runtime = 90  # days of runtime (this is only used for display in the plot)
bp = .85  # Beam polarization (this is only used for display in the plot)

# Load the data
Bins = pd.read_csv('./Files/' + beamE + 'GeV_files/CT18NLO_' + beamE + 'GeV_from' + NQ2bins + 'x' + Nxbins + 'grid_analytic_values_for_each_x.csv')#, nrows=3)  # use nrows=3 to test the program
print("Bins:", Bins)
Bins.rename(columns = {'Q2 (GeV^2)':'Q2'}, inplace = True)
# Where to save the data
savefile_dir = './Files/' + beamE + '\
GeV_files/'
savefile_append = '_' + beamE + 'GeV_q_PDF_vals_1x' + Nxbins + 'grid.csv'

# Declare the PDFs to be used
Grid = ['NNPDF40_nlo_as_01180', 'PDF4LHC21_40', 'CT18NLO']#, 'NNPDF31_nlo_as_0118']  # Names of the PDFs you wish to use
N_PDF = [101, 41, 59]#, 101]  # Number of PDFs for each Grid
base = ["replica", "Hessian", "Hessian"]#, "replica"]  # base of the Grid

# Turn on/off the following capabilities
find_unc = True  # Calculate the PDF uncertainties
make_plots = False  # Plot graphs
save_data = True  # Save PDF values (& PDF unc) to csv files

### END USER INPUT ###

alp = 0.4  # transparency of error bands
fig, ax = plt.subplots()
#PDF_vals = pd.DataFrame(data=Bins.loc[:,'x'], columns='x')
all_data = {}  # dict for storing dataframes for each PDF

# Define constants
G_F = 1.16638e-5  # GeV^2 Fermi constant
alpha = 7.29735256e-3  # Fine structure constant

for i, name in enumerate(Grid):
    print("\n\n PDF number {}, ".format(i) + name)
    q_PDF_vals = []  # quark PDF values
    q_PDF_unc = []
    for index, row in Bins.iterrows():
        q_PDF_vals.append(get_q_PDF_vals(row, name, 0))
        if find_unc:
            q_PDF_unc.append(calc_unc(row, name, N_PDF[i], base[i]))
        print(index, row)

    q_PDF_vals = np.array(q_PDF_vals)  # Change it to an array for easy math

    # Calculate and store  valance quark and total quark PDF values
    PDF_vals = pd.DataFrame(data=Bins.loc[:,'x'], columns=['x'])
    PDF_vals["d_V"] = q_PDF_vals[:,0] - q_PDF_vals[:,1]
    PDF_vals["u_V"] = q_PDF_vals[:,2] - q_PDF_vals[:,3]
    PDF_vals["s_V"] = q_PDF_vals[:,4] - q_PDF_vals[:,5]
    PDF_vals["c_V"] = q_PDF_vals[:,6] - q_PDF_vals[:,7]
    PDF_vals["d_p"] = q_PDF_vals[:,0] + q_PDF_vals[:,1]
    PDF_vals["u_p"] = q_PDF_vals[:,2] + q_PDF_vals[:,3]
    PDF_vals["s_p"] = q_PDF_vals[:,4] + q_PDF_vals[:,5]
    PDF_vals["c_p"] = q_PDF_vals[:,6] + q_PDF_vals[:,7]
    # Calculate and store the PDF value for d_V/u_V
    PDF_vals["dV/uV"] = PDF_vals["d_V"] / PDF_vals["u_V"] 

    # Store quark PDF values
    PDF_vals["d"] = q_PDF_vals[:,0]
    PDF_vals["dbar"] = q_PDF_vals[:,1]
    PDF_vals["u"] = q_PDF_vals[:,2]
    PDF_vals["ubar"] = q_PDF_vals[:,3]
    PDF_vals["s"] = q_PDF_vals[:,4]
    PDF_vals["sbar"] = q_PDF_vals[:,5]
    PDF_vals["c"] = q_PDF_vals[:,6]
    PDF_vals["cbar"] = q_PDF_vals[:,7]
    
    if find_unc:
        # Store PDF uncertainties for each quark
        q_PDF_unc = np.array(q_PDF_unc)  # Changed to an array for easy math
        PDF_vals["d_unc"] = q_PDF_unc[:,0]
        PDF_vals["dbar_unc"] = q_PDF_unc[:,1]
        PDF_vals["u_unc"] = q_PDF_unc[:,2]
        PDF_vals["ubar_unc"] = q_PDF_unc[:,3]
        PDF_vals["s_unc"] = q_PDF_unc[:,4]
        PDF_vals["sbar_unc"] = q_PDF_unc[:,5]
        PDF_vals["c_unc"] = q_PDF_unc[:,6]
        PDF_vals["cbar_unc"] = q_PDF_unc[:,7]
        # Calculate and store the PDF uncertainty for d_V/u_V
        dV_unc = np.sqrt(PDF_vals["d_unc"] ** 2 + PDF_vals["dbar_unc"] ** 2)
        uV_unc = np.sqrt(PDF_vals["u_unc"] ** 2 + PDF_vals["ubar_unc"] ** 2)
        PDF_vals["(dV/uV)_unc"] = abs(PDF_vals["dV/uV"]) * \
            np.sqrt((dV_unc / PDF_vals["d_V"]) ** 2 + 
                    (uV_unc / PDF_vals["u_V"]) ** 2)
    
    # Store the dataframe in a dictionary labeled according to the PDF set used
    all_data[name] = PDF_vals
    print(PDF_vals)
    
    if make_plots:
        if find_unc:
            # Plot the error bands
            ax.fill_between(PDF_vals["x"], np.array(PDF_vals["dV/uV"]) + 
                            np.array(PDF_vals["(dV/uV)_unc"]), 
                            np.array(PDF_vals["dV/uV"]) - 
                            np.array(PDF_vals["(dV/uV)_unc"]), alpha=alp)
        ax.plot(Bins.loc[:,'x'], PDF_vals["dV/uV"], label=name)  # Plot dV/uV

print()
print()
print("all_data", all_data)

if save_data:
    # Save the PDF values for each PDF to a CSV file
    for name in all_data.keys():
        all_data[name].to_csv(savefile_dir + name + savefile_append, index=False) #this will include a header. 

if make_plots:
    # Plot extra info for help in understanding graph
    ax.plot([.15,.85], [0,0], '--k', label="zero")
    ax.text(.15, -.1, "Beam Polarization: " + str(bp))
    ax.text(.15, .05, "Runtime: " + str(runtime) + " days")
    #ax.text(.15, .75, "$\delta O = \sqrt{\\frac{1}{4}\sum_{k=1}^{k=\\frac{N_{PDF}}{2}}(O_{2k}-O_{2k-1})^2}$")  # Uncertainty equation for Hessian based PDFs

    # Plot labels and formatting
    plt.ylim(-1,1)
    plt.xlabel("x")
    plt.ylabel("dV/uV")
    plt.legend(loc='lower left')
    plt.show()
