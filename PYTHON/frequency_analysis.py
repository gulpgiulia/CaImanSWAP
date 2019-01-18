def bandpass_filter(data, order, fs, lowcut, highcut): # data has to be a 3D array

    from scipy.signal import butter, lfilter
    from numpy import reshape
    import numpy as np

    # Filter requirements.
    #order = 6     # order of the Butterworth filter
    #fs = 25.0       # sample rate, Hz
    #lowcut = 0.3  # desired cutoff frequency of the high-pass filter, Hz
    #highcut = 2.5  # desired cutoff frequency of the low-pass filter, Hz

    #Definition of the Butterworth bandpass filter

    def butter_bandpass_filter(data, lowcut, highcut, fs, order):
        nyq = 0.5 * fs
        low = lowcut / nyq
        high = highcut / nyq

        b, a = butter(order, [low, high], btype='band')
        y = lfilter(b, a, data)
        return y

    #Application of the filter to the dataset

    clean_signals = []

    for i in range (np.size(data,1)):
        for j in range (np.size(data,2)):
            clean_signals.append([butter_bandpass_filter(data[:, i, j], lowcut, highcut, fs, order)])

    clean_signals = reshape(clean_signals, (np.size(data,1), np.size(data,2), np.size(data,0)))

    return clean_signals



#-------------------------------------------------------------------------------

def multiplot(xstart, xstop, xsteps, ystart, ystop, ysteps, tmin, tmax):

    import matplotlib.pyplot as plt
    import numpy as np

    xsize = len(np.arange(xstart, xstop, xsteps))
    ysize = len(np.arange(ystart, ystop, ysteps))

    # Multiplot construction
    f, axarr = plt.subplots(xsize, ysize, figsize=(5, 5))
    k = -1
    for i in range (xstart, xstop, xsteps):
        k += 1
        l = -1
        for j in range (ystart, ystop, ysteps):
            l += 1
            axarr[k, l].plot(clean_signals[i, j, tmin:tmax])
            axarr[k, l].set_title('Pixel [%s,%s]' % (i, j), fontsize = 8)


    # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
    plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
    plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)
    plt.setp([a.get_yticklabels() for a in axarr[:, 2]], visible=False)
    plt.show();
