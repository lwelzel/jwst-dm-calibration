"""
Plotters taken from previous work in a different course at
Leiden University, The Netherlands
2022
By M. Kenworthy, L. Welzel
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rc('image', interpolation='nearest', origin='lower')


def displC(c, title="Amplitude-Phase space of the PSF", trim=0, show=True):
    """displC(c,trim=0) - display a Complex number c as four plots
    The bottom two plots are (Amplitude, Phase)
    Optionally cut out the central square  with a size of 'trim x trim' pixels"""
    c2 = np.copy(c)
    if (trim > 0):  # if the user specifies a trim value, cut out the centre of the image
        (nx, ny) = c.shape
        dx = int((nx - trim) / 2 + 1)
        dy = int((nx - trim) / 2 + 1)
        c2 = c[dx:dx + trim, dy:dy + trim]

    fig, ((axre, axim), (axamp, axpha)) = plt.subplots(nrows=2, ncols=2,
                                                       constrained_layout=True,
                                                       sharex=True,
                                                       sharey=True,
                                                       figsize=(12, 12))

    # plot out the panels
    im = axre.imshow(c2.real)
    plt.colorbar(im, ax=axre, fraction=0.046, pad=0.04, label="")
    im = axim.imshow(c2.imag)
    plt.colorbar(im, ax=axim, fraction=0.046, pad=0.04, label="")
    im = axamp.imshow(np.abs(c2))
    plt.colorbar(im, ax=axamp, fraction=0.046, pad=0.04, label="")
    im = axpha.imshow(np.angle(c2))
    plt.colorbar(im, ax=axpha, fraction=0.046, pad=0.04, label="")

    axre.set_title('Real')
    axre.set_ylabel('y [px]')
    axim.set_title('Imag')
    axamp.set_title('Amplitude')
    axamp.set_xlabel('x [px]')
    axamp.set_ylabel('y [px]')
    axpha.set_title('Phase')
    axpha.set_xlabel('x [px]')
    plt.suptitle(title, fontsize=18)
    if show:
        plt.show()


def plotCAmpl(Iorig, GAMMA=4, RANGE=None, plotnow=True, trim=0, title="plotCAmpl output"):
    """plotCAmpl - display a complex valued image by taking a complex 2D array
        as input and returning an [N,3] RGB image for display

        if RANGE=(min,max) is given, it will set image intensity
        to 0 and 1 for these values and BLACK and WHITE
        colours for out of range amplitude values respectively

        the returned image has AMPLITUDE as brightness and
        PHASE as colour:

        RED    YELLOW   BLUE   GREEN
        0        90     180     270  degrees phase

        The image is scaled by a GAMMA factor (defaults to 4)
        to emphasize smaller amplitudes

        Set plotnow=False to stop display

        Returns a MxNx3 image suitable for display
    """
    import numpy as np
    PI = 4 * np.arctan(1)
    I = np.copy(Iorig)
    if (trim > 0):  # if the user specifies a trim value, cut out the centre of the image
        (nx, ny) = I.shape
        dx = int((nx - trim) / 2)
        dy = int((nx - trim) / 2)
        I = I[dx:dx + trim, dy:dy + trim]

    # calculate amplitude
    ampl = (abs(I)).real

    if RANGE is None:
        ampl = ampl / np.amax(ampl)
    else:
        ampl = (ampl - RANGE[0]) / abs(RANGE[1] - RANGE[0])
        ampOVER = (ampl > 1)
        ampUNDER = (ampl < 0)
        ampl[ampOVER] = 1
        ampl[ampUNDER] = 0

    # apply gamma to amplitude
    ampl = np.power(ampl, 1. / GAMMA)

    # normalise phase from 0 to 1
    phse = (np.angle(I) / (2 * PI)) + 0.5

    # make an RGB array
    img = np.zeros((I.shape[0], I.shape[1], 3), dtype=float)

    # red
    red = abs(phse - 0.625) * 4 - 0.5
    red = np.maximum(red, 0)
    red = np.minimum(red, 1)

    # green
    grn = 2 - abs(phse - 0.5) * 4
    tmp = np.where(grn > 1)
    grn[tmp] = 2. - grn[tmp]

    # blue
    blu = 1 - abs(phse - 0.5) * 4
    blu = np.maximum(blu, 0)

    img[:, :, 0] = ampl * red
    img[:, :, 1] = ampl * grn
    img[:, :, 2] = ampl * blu

    titles = ['red', 'green', 'blue', 'rgb image']

    # mark out of range points with RGB color
    if RANGE:
        img[ampOVER, :] = (1, 1, 1)
        img[ampUNDER, :] = (0, 0, 0)

    if plotnow:
        # modified for ipython notebooks
        fig, ax = plt.subplots(1, figsize=(6, 6))
        ax.set_xlabel("x");
        ax.set_ylabel("y")
        ax.set_title("{}".format(title))
        ax.imshow(img, aspect='equal', vmin=0., vmax=1., interpolation='nearest')

    return (img)


def padcplx(c, pad=5):
    """padcplx - puts a `Complex` array into the centre of a zero-filled `Complex` array
               pad is an integer defining the padding multiplier for the output array """
    (nx, ny) = c.shape
    bignx = nx * pad
    bigny = ny * pad
    big_c = np.zeros((bignx, bigny), dtype=complex)

    dx = int((nx * (pad - 1)) / 2)
    dy = int((ny * (pad - 1)) / 2)

    big_c[dx:dx + nx, dy:dy + ny] = c
    return (big_c)


def FFT(c, pad=5):
    """FFT - carry out a complex Fourier transform (with optional padding)
            c - the input 2D Complex numpy array
            pad - integer multiplier for the padding/sampling
            Returns the `Complex` FFT padded array"""
    from numpy.fft import fft2, fftshift, ifft2, ifftshift
    psfA = fftshift(fft2(ifftshift(padcplx(c, pad))))
    return psfA


def IFFT(cb, pad=5):
    """IFFT - carry out the complex Fourier transform (with optional padding)
            and return the FFT padded array"""
    from numpy.fft import fft2, fftshift, ifft2, ifftshift
    psfB = fftshift(ifft2(ifftshift(padcplx(cb, pad))))
    return psfB
