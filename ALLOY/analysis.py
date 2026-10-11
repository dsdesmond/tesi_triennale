import numpy as np
from matplotlib import pyplot as plt
import sys
import scipy
from scipy.stats import t

#matplotlib.use('TKAgg')
#matplotlib.use("module://backend_interagg")

def main(GLA_ALL = "GLA_ALL", KAW_ALL = "KAW_ALL", fit_graphs = True):
    """
        Function to graphically analyse the Glauber and Kawasaki dynamics Ising model simulation

        Uses results saved to files where each column is an ising lattice evaluated at a different temperature, and
        each row is as follows: T  M  E  c  Chi  c_err


        Parameters
        ----------

        GLA_ALL : filename
            filename for Glauber dynamics simulation
        KAW_ALL : filename
            filename for Kawasaki dynamics simulation


        Returns
        -------


        """









    GLA_ALL = np.loadtxt("G_results/"+GLA_ALL)
    #GLA_M_ALL = np.loadtxt("G_results/GLA_M_ALL")
    #GLA_E_ALL = np.loadtxt("G_results/GLA_E_ALL")

    GLA_T_arr = GLA_ALL[0]
    GLA_M_arr = abs(GLA_ALL[1])
    GLA_E_arr = GLA_ALL[2]
    GLA_hc_arr = GLA_ALL[3]
    GLA_susc_arr = GLA_ALL[4]
    GLA_hc_err_arr = GLA_ALL[5]



    KAW_ALL = np.loadtxt("G_results/"+KAW_ALL)
    #GLA_M_ALL = np.loadtxt("G_results/GLA_M_ALL")
    #GLA_E_ALL = np.loadtxt("G_results/GLA_E_ALL")

    KAW_T_arr = KAW_ALL[0]
    KAW_M_arr = abs(KAW_ALL[1])
    KAW_E_arr = KAW_ALL[2]
    KAW_hc_arr = KAW_ALL[3]
    KAW_susc_arr = KAW_ALL[4]
    KAW_hc_err_arr = KAW_ALL[5]






    #PLOTTING ENERGY -----------------------------------------------------------

    if fit_graphs:
        def sigmoid_function(T, T_c, a, b, k):
            y = a / (1 + np.exp(-k * (T - T_c))) + b  # L / (1 + np.exp(-k*(x-x0))) + b
            return y

        p_0 = [0, -2500, 0, 0]
        x_ = np.arange(1, 3, 0.01)

        popt, pcov = scipy.optimize.curve_fit(sigmoid_function, GLA_T_arr, GLA_E_arr, p_0, method='dogbox')

        plt.plot(x_, sigmoid_function(x_, popt[0], popt[1], popt[2], popt[3]), color = "b", label = "Glauber sigmoid fit", alpha = 0.4)
        print("Critical temperature as estimated by Glauber energy: " + str(round(popt[0], 2)) + "+/-" + str(
            np.round(np.sqrt(np.diag(pcov))[0], 3)))

    plt.plot(GLA_T_arr, GLA_E_arr, marker = "x", linewidth = 0, color = "b", label = "Glauber im results")




    if fit_graphs:
        popt, pcov = scipy.optimize.curve_fit(sigmoid_function, KAW_T_arr, KAW_E_arr, p_0, method='dogbox')

        plt.plot(x_, sigmoid_function(x_, popt[0], popt[1], popt[2], popt[3]), color = "coral", label = "Kawasaki sigmoid fit", alpha = 0.4)

        print("Critical temperature as estimated by Kawasaki energy: " + str(round(popt[0], 2)) + "+/-" + str(
            np.round(np.sqrt(np.diag(pcov))[0], 3)))



    plt.plot(KAW_T_arr, KAW_E_arr, marker = "x", linewidth = 0, color = "coral", label = "Kawasaki results")

    plt.title("Total Energy against temperature")
    plt.ylabel("Total Energy")
    plt.xlabel("Temperature (K)")
    plt.legend()
    plt.show()






    #PLOTTING MAGNETISATION -----------------------------------------------------------

    if fit_graphs:

        p_0 = [0, 2500, 0, 0]
        x_ = np.arange(1, 3, 0.01)

        popt, pcov = scipy.optimize.curve_fit(sigmoid_function, GLA_T_arr, GLA_M_arr, p_0, method='dogbox')

        plt.plot(x_, sigmoid_function(x_, popt[0], popt[1], popt[2], popt[3]), color="b", label="Glauber sigmoid fit", alpha = 0.4)
        print("Critical temperature as estimated by Glauber magnetism: " + str(round(popt[0], 2)) + "+/-" + str(
            np.round(np.sqrt(np.diag(pcov))[0], 3)))
        plt.text(1, 0, "Note that the fit does not approximate \n"
                       "the distribution, it is merely a crude \n"
                       "way to evaluate the critical temperature ",
                 fontsize=6)  # , transform=ax.transAxes, fontsize=14# verticalalignment='top', bbox=props)


    plt.title("Magnetisation against temperature")
    plt.ylabel("Magnetisation")
    plt.xlabel("Temperature (K)")
    plt.plot(KAW_T_arr, KAW_M_arr, marker = "x", color = "coral", label="Kawasaki")
    plt.plot(GLA_T_arr, GLA_M_arr, "bx", label = "Glauber")

    plt.legend()
    plt.show()












    #PLOTTING HEAT CAPACITY -----------------------------------------------------------

    if fit_graphs:
        def students_t(x, mu, sigma, df, A):
            return A * t.pdf((x - mu) / sigma, df)
        def gaussian(x, mu, sigma, A):
            return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


    x_ = np.arange(1,3,0.01)




    if fit_graphs:
        p0 = [np.mean(GLA_T_arr), np.std(GLA_T_arr), 1, np.max(GLA_hc_arr) - np.min(GLA_hc_arr)]
        popt, pcov = scipy.optimize.curve_fit(students_t, GLA_T_arr, GLA_hc_arr, p0=p0, method = "dogbox")
        plt.plot(x_,students_t(x_, popt[0],popt[1],popt[2],popt[3]), color = "b", label = "Glauber students_t fit", alpha = 0.4)

        #print("Note that the fit does not necessarily approximate the distribution, it is a very crude way to "
        #      "determine the Critical temperature but it will do for our purposes")
        print("Critical temperature as estimated by Glauber heat capacity: " + str(round(popt[0], 2)) + "+/-" + str(
            np.round(np.sqrt(np.diag(pcov))[0], 3)))

    plt.title("Heat capacity against temperature")
    plt.ylabel("Heat Capacity")
    plt.xlabel("Temperature (K)")
    #plt.plot(GLA_T_arr, GLA_hc_arr, label = "Glauber")
    plt.errorbar(GLA_T_arr, GLA_hc_arr, yerr=GLA_hc_err_arr,  marker = "x",  color = "b", label="Glauber simulation")
    #plt.plot(KAW_T_arr, KAW_hc_arr, label = "Kawasaki")

    if fit_graphs:

        p0 = [np.mean(KAW_T_arr), np.std(KAW_T_arr), 1, np.max(KAW_hc_arr) - np.min(KAW_hc_arr)]
        popt, pcov = scipy.optimize.curve_fit(students_t, KAW_T_arr, KAW_hc_arr, p0=p0, method = "dogbox")
        plt.plot(x_,students_t(x_, popt[0],popt[1],popt[2],popt[3]), color = "coral", label = "Kawasaki students_t fit", alpha = 0.4)


        #print("Note that the fit does not necessarily approximate the distribution, it is a very crude way to "
        #      "determine the Critical temperature but it will do for our purposes")
        print("Critical temperature as estimated by Kawasaki heat capacity: " + str(round(popt[0], 2)) + "+/-" + str(
            np.round(np.sqrt(np.diag(pcov))[0], 3)))
        plt.text(1.8, 0, "Note that the fit does not approximate \n"
                         "the distribution, it is merely a crude \n"
                         "way to evaluate the critical temperature ",
                 fontsize=6)  # , transform=ax.transAxes, fontsize=14# verticalalignment='top', bbox=props)




    plt.errorbar(KAW_T_arr, KAW_hc_arr, yerr = KAW_hc_err_arr, marker = "x", color = "coral",label = "Kawasaki simulation")


    plt.legend()
    plt.show()


    #PLOTTING SUSC -----------------------------------------------------------





    plt.title("Susceptability against temperature")
    plt.ylabel("Susceptability")
    plt.xlabel("Temperature (K)")
    plt.plot(GLA_T_arr, GLA_susc_arr, marker = "x", color = "b", label = "Glauber")
    plt.plot(KAW_T_arr, KAW_susc_arr, marker = "x", color = "coral",  label = "Kawasaki")
    plt.legend()
    plt.show()



if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:   Glauber_data_filename   Kawasaki_data_filename    fit_graphs")
        print("For further usage instructions, see README or documentation. For some data, graph fitting will struggle. If so, set fit_graphs to False")
        sys.exit(1)



    main(sys.argv[1],sys.argv[2],eval(sys.argv[3]))