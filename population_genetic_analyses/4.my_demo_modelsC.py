import numpy
import dadi

def split_grow((nuCur0, TCur, nuCurF, nuCriF, s, T2, m12, m21), (n1,n2), pts):

    xx = yy = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    # The first size change
    phi = dadi.Integration.one_pop(phi, xx, TCur, nu=nuCur0)

    # Now the Out-of-Africa divergence 
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # Both Criollo and Curaray have changes in pop size  
    nuCriF_func = lambda t: s*nuCur0 * (nuCriF/(s*nuCur0))**(t/T2)
    nuCurF_func = lambda t: nuCur0 * (nuCurF/nuCur0)**(t/T2)
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nuCriF_func, nu2=nuCurF_func, m12=m12, m21=m21)

    # Make the Spectrum object 
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def IM(params, (n1,n2), pts):
    nu1, nu2, m12, m21, Ts = params
    """
    Model with migration during the divergence.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    m12: Migration from pop 2 to pop 1 (2*Na*m12).
    m21: Migration from pop 1 to pop 2.
    Ts: The scaled time between the split and the secondary contact (in units of 2*Na generations).
    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = dadi.Numerics.default_grid(pts)
    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the divergence event
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We set the population sizes after the split to nu1 and nu2 and set the migration rates to m12 and m21
    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=m12, m21=m21)
    # Finally, calculate the spectrum.
    fs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return fs


##########################
## Some other models we explored
###########################


import numpy
import dadi

def prior_onegrow_mig((nu1B, nu1F, nu2B, nu2F, m, Tp, T), (n1,n2), pts):
    """
    Model with growth, split, bottleneck in pop2 , exp recovery, migration

    nu1F: The ancestral population size after growth. (Its initial size is
          defined to be 1.)
    nu2B: The bottleneck size for pop2
    nu2F: The final size for pop2
    m: The scaled migration rate
    Tp: The scaled time between ancestral population growth and the split.
    T: The time between the split and present

    n1,n2: Size of fs to generate.
    pts: Number of points to use in grid for evaluation.
    """
    # Define the grid we'll use
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    # Now do the population growth event.
    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nu1F)

    # The divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    # We need to define a function to describe the non-constant population 2
    # size. lambda is a convenient way to do so.
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T)
    nu1_func = lambda t: nu1B*(nu1F/nu1B)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1_func, nu2=nu2_func,
                                    m12=m, m21=m)


    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,yy))
    return sfs

def prior_onegrow_mig_mscore((nu1B, nu1F, nu2B, nu2F, m, Tp, T)):
    """
    ms core command corresponding to prior_onegrow_mig
    """
    # Growth rate
    alpha2 = numpy.log(nu2F/nu2B)/T
    alpha1 = numpy.log(nu1F/nu1B)/T
    command = "-n 1 %(nu1F)f -n 2 %(nu2F)f "\
            "-eg 0 2 %(alpha2)f "\
            "-eg 0 1 %(alpha1)f "\
            "-ma x %(m)f %(m)f x "\
            "-ej %(T)f 2 1 "\
            "-en %(Tsum)f 1 1"

    # There are several factors of 2 necessary to convert units between dadi
    # and ms.
    sub_dict = {'nu1F':nu1F, 'nu2F':nu2F, 'alpha2':2*alpha2,'alpha1':2*alpha1,
                'm':2*m, 'T':T/2, 'Tsum':(T+Tp)/2}

    return command % sub_dict


