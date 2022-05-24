#! /usr/bin/env python3
from cffi import FFI
import io
import numpy
import os
from pcpp.preprocessor import Preprocessor
import re
from scipy.interpolate import PchipInterpolator


# Cross-section references
references = {
    "CMS11" : "https://arxiv.org/abs/1106.3723",
    "BGR18" : "https://arxiv.org/abs/1808.02034"
}


def cross_section_CMS11():
    """CMS11 cross-section

       Reference:
         https://arxiv.org/abs/1106.3723
    """
    cs0 = numpy.array((
        (50, 0.32, 4.1, -2.3, -2.4, 0.10, 3.8, -1.9, -2.0),
        (100, 0.65, 3.8, -2.0, -2.0, 0.20, 3.5, -1.8, -1.8),
        (200, 1.3, 3.5, -1.8, -1.9, 0.41, 3.2, -1.6, -1.7),
        (500, 3.2, 3.2, -1.7, -1.8, 1.0, 2.9, -1.5, -1.5),
        (1000, 6.2, 3.0, -1.6, -1.7, 2.0, 2.7, -1.4, -1.5),
        (2000, 12., 2.7, -1.6, -1.6, 3.8, 2.4, -1.3, -1.4),
        (5000, 27., 2.3, -1.5, -1.5, 8.6, 2.1, -1.3, -1.3),
        (10000, 47., 2.0, -1.4, -1.4, 15., 1.8, -1.2, -1.2),
        (20000, 77., 1.8, -1.3, -1.4, 26., 1.6, -1.1, -1.1),
        (50000, 140., 1.5, -1.2, -1.2, 49., 1.3, -1.0, -1.1),
        (100000, 210., 1.4, -1.2, -1.2, 75., 1.2, -1.0, -1.0),
        (200000, 310., 1.5, -1.1, -1.1, 110., 1.2, -0.9, -0.9),
        (500000, 490., 1.6, -1.0, -1.0, 180., 1.3, -0.8, -0.8),
        (1E+06, 690., 1.7, -0.9, -0.9, 260., 1.4, -0.8, -0.8),
        (2E+06, 950., 1.9, -0.9, -0.9, 360., 1.6, -0.8, -0.8),
        (5E+06, 1400., 2.0, -0.9, -0.9, 540., 1.8, -0.8, -0.8),
        (1E+07, 1900., 2.2, -0.9, -0.9, 730., 2.0, -0.8, -0.8),
        (2E+07, 2600., 2.3, -0.9, -1.0, 980., 2.2, -0.8, -0.9),
        (5E+07, 3700., 2.5, -0.9, -1.2, 1400., 2.4, -0.9, -1.1),
        (1E+08, 4800., 2.7, -0.9, -1.5, 1900., 2.6, -0.9, -1.3),
        (2E+08, 6200., 2.8, -1.0, -2.0, 2400., 2.7, -1.0, -1.8),
        (5E+08, 8700., 3.0, -1.1, -3.0, 3400., 2.9, -1.0, -2.6),
        (1E+09, 11000., 3.1, -1.2, -3.9, 4400., 3.0, -1.1, -3.4),
        (2E+09, 14000., 3.3, -1.2, -5.0, 5600., 3.2, -1.2, -4.4),
        (5E+09, 19000., 3.4, -1.4, -6.8, 7600., 3.4, -1.3, -6.1),
        (1E+10, 24000., 3.6, -1.5, -8.5, 9600., 3.5, -1.4, -7.6),
        (2E+10, 30000., 3.7, -1.6, -10.3, 12000., 3.6, -1.5, -9.3),
        (5E+10, 39000., 3.8, -1.7, -13.1, 16000., 3.8, -1.7, -11.8),
        (1E+11, 48000., 4.0, -1.8, -15.2, 20000., 3.9, -1.8, -13.9),
        (2E+11, 59000., 4.1, -1.9, -17.5, 24000., 4.0, -1.9, -16.1),
        (5E+11, 75000., 4.2, -2.0, -20.3, 31000., 4.2, -2.0, -18.8)
    ))

    cs1 = numpy.array((
        (50, 0.15, 15.0, -9.0, -9.0, 0.05, 12.0, -6.4, -6.4),
        (100, 0.33, 13.3, -7.4, -7.4, 0.12, 10.7, -5.7, -5.7),
        (200, 0.69, 11.9, -6.5, -6.5, 0.24, 9.6, -5.1, -5.1),
        (500, 1.8, 10.5, -5.7, -5.7, 0.61, 8.6, -4.6, -4.6),
        (1000, 3.6, 9.4, -5.2, -5.2, 1.20, 7.8, -4.2, -4.2),
        (2000, 7., 8.3, -4.6, -4.6, 2.4, 7.0, -3.8, -3.8),
        (5000, 17., 6.5, -3.7, -3.7, 5.8, 5.7, -3.2, -3.2),
        (10000, 31., 5.1, -3.0, -3.0, 11., 4.6, -2.7, -2.7),
        (20000, 55., 3.8, -2.3, -2.3, 19., 3.6, -2.1, -2.1),
        (50000, 110., 2.5, -1.7, -1.7, 39., 2.4, -1.5, -1.5),
        (100000, 180., 1.9, -1.4, -1.4, 64., 1.7, -1.2, -1.2),
        (200000, 270., 1.7, -1.2, -1.2, 99., 1.4, -1.0, -1.0),
        (500000, 460., 1.7, -1.1, -1.1, 170., 1.4, -0.9, -0.9),
        (1E+06, 660., 1.8, -1.0, -1.0, 240., 1.5, -0.8, -0.8),
        (2E+06, 920., 1.9, -1.0, -1.0, 350., 1.6, -0.8, -0.8),
        (5E+06, 1400., 2.1, -0.9, -0.9, 530., 1.9, -0.8, -0.8),
        (1E+07, 1900., 2.2, -0.9, -0.9, 730., 2.0, -0.8, -0.8),
        (2E+07, 2500., 2.3, -0.9, -1.0, 980., 2.2, -0.8, -0.9),
        (5E+07, 3700., 2.5, -0.9, -1.2, 1400., 2.4, -0.9, -1.1),
        (1E+08, 4800., 2.7, -1.0, -1.5, 1900., 2.6, -0.9, -1.3),
        (2E+08, 6200., 2.8, -1.0, -2.0, 2400., 2.7, -1.0, -1.8),
        (5E+08, 8700., 3.0, -1.1, -3.0, 3400., 2.9, -1.0, -2.6),
        (1E+09, 11000., 3.1, -1.2, -3.9, 4400., 3.0, -1.1, -3.4),
        (2E+09, 14000., 3.3, -1.2, -5.0, 5600., 3.2, -1.2, -4.4),
        (5E+09, 19000., 3.4, -1.4, -6.8, 7600., 3.4, -1.3, -6.1),
        (1E+10, 24000., 3.6, -1.5, -8.5, 9600., 3.5, -1.4, -7.6),
        (2E+10, 30000., 3.7, -1.6, -10.3, 12000., 3.6, -1.5, -9.3),
        (5E+10, 39000., 3.8, -1.7, -13.1, 16000., 3.8, -1.7, -11.8),
        (1E+11, 48000., 4.0, -1.8, -15.2, 20000., 3.9, -1.8, -13.9),
        (2E+11, 59000., 4.1, -1.9, -17.5, 24000., 4.0, -1.9, -16.1),
        (5E+11, 75000., 4.2, -2.0, -20.3, 31000., 4.2, -2.0, -18.8)
    ))

    return numpy.vstack((cs0[:,0], cs0[:,1], cs0[:,5], cs1[:,1], cs1[:,5])).T


def cross_section_BGR18():
    """BGR18 cross-section

       Reference:
         https://arxiv.org/abs/1808.02034
    """

    cs0 = numpy.array((
        (5E+03,     25.5,     15.2),
        (1E+04,     44.6,     28.5),
        (2E+04,     73.8,     51.2),
        (5E+04,      133,      103),
        (1E+05,      199,      165),
        (2E+05,      287,      252),
        (5E+05,      453,      421),
        (1E+06,      628,      600),
        (2E+06,      859,      837),
        (5E+06, 1.28E+03, 1.27E+03),
        (1E+07, 1.71E+03, 1.71E+03),
        (2E+07, 2.27E+03, 2.28E+03),
        (5E+07, 3.26E+03, 3.28E+03),
        (1E+08, 4.25E+03, 4.29E+03),
        (2E+08, 5.51E+03, 5.56E+03),
        (5E+08, 7.69E+03, 7.76E+03),
        (1E+09, 9.82E+03, 9.93E+03),
        (2E+09, 1.25E+04, 1.26E+04),
        (5E+09, 1.70E+04, 1.72E+04),
        (1E+10, 2.14E+04, 2.17E+04),
        (2E+10, 2.69E+04, 2.72E+04),
        (5E+10, 3.60E+04, 3.64E+04),
        (1E+11, 4.47E+04, 4.52E+04),
        (2E+11, 5.54E+04, 5.61E+04),
        (5E+11, 7.32E+04, 7.41E+04),
        (1E+12, 9.00E+04, 9.12E+04),
        (2E+12, 1.10E+05, 1.12E+05),
        (5E+12, 1.44E+05, 1.45E+05)
    ))

    cs1 = numpy.array((
        (5E+03,     8.45,     5.52),
        (1E+04,     15.1,     10.4),
        (2E+04,     25.8,       19),
        (5E+04,     48.6,     38.9),
        (1E+05,     74.6,     63.5),
        (2E+05,      111,     99.5),
        (5E+05,      182,      170),
        (1E+06,      258,      248),
        (2E+06,      361,      352),
        (5E+06,      552,      545),
        (1E+07,      751,      746),
        (2E+07, 1.01E+03, 1.01E+03),
        (5E+07, 1.48E+03, 1.47E+03),
        (1E+08, 1.95E+03, 1.95E+03),
        (2E+08, 2.55E+03, 2.55E+03),
        (5E+08, 3.6E+03,   3.6E+03),
        (1E+09, 4.63E+03, 4.63E+03),
        (2E+09, 5.93E+03, 5.93E+03),
        (5E+09, 8.15E+03, 8.15E+03),
        (1E+10, 1.03E+04, 1.03E+04),
        (2E+10, 1.3E+04,   1.3E+04),
        (5E+10, 1.75E+04, 1.75E+04),
        (1E+11, 2.18E+04, 2.18E+04),
        (2E+11, 2.71E+04, 2.71E+04),
        (5E+11, 3.6E+04,   3.6E+04),
        (1E+12, 4.44E+04, 4.44E+04),
        (2E+12, 5.46E+04, 5.46E+04),
        (5E+12, 7.14E+04, 7.14E+04)
    ))

    return numpy.vstack((cs0[:,0], cs0[:,1], cs1[:,1], cs0[:,2], cs1[:,2])).T


def initialise_ent():
    """Initialise the ENT library
    """
    with open("include/ent.h") as f:
        header = f.read()

    # Prune the header
    header = re.sub(r"(?m)^#include.*\n?", "", header)
    cpp = Preprocessor()
    cpp.parse(header)
    output = io.StringIO()
    cpp.write(output)
    header = output.getvalue()

    ffi = FFI()
    ffi.cdef(header)

    # Load the library
    lib = ffi.dlopen("lib/libent.so")

    return ffi, lib

ffi, lib = initialise_ent()


def interpolate_cross_section(model):
    """Interpolate using PCHIP algorithm
    """
    cs = model()

    p00 = PchipInterpolator(x=numpy.log(cs[:,0]), y=numpy.log(cs[:,1]))
    p01 = PchipInterpolator(x=numpy.log(cs[:,0]), y=numpy.log(cs[:,2]))
    p10 = PchipInterpolator(x=numpy.log(cs[:,0]), y=numpy.log(cs[:,3]))
    p11 = PchipInterpolator(x=numpy.log(cs[:,0]), y=numpy.log(cs[:,4]))

    return p00, p01, p10, p11


def build_cross_section_table(model):
    """Build the cross-section tabulation, for ENT
    """

    # Create the physics
    physics = ffi.new("struct ent_physics *[1]")
    lib.ent_physics_create(
        physics, f"share/ent/{model}-sf.ent".encode(), ffi.NULL)

    def cross_section(projectile, energy, Z, A, process):
        cs = ffi.new("double [1]")
        lib.ent_physics_cross_section(
            physics[0], projectile, energy, Z, A, process, cs);
        return float(cs[0])

    p00, p01, p10, p11 = interpolate_cross_section(
        globals()[f"cross_section_{model}"]
    )

    with open(f"share/ent/{model}-cross-section.txt", "w+") as f:
        f.write(f"""# {model} cross-section
#
# Generated using ENT.
#
# Reference:
#    {references[model]}
#
#------------------------------------------------------------------------------------------------------------
#  energy                 neutrino xsec (m^2)                            anti-neutrino xsec (m^2)
#  (GeV)         charged-current        neutral-current           charged-current        neutral-current
#              proton      neutron     proton      neutron      proton      neutron     proton      neutron
#------------------------------------------------------------------------------------------------------------
""")
        energies = numpy.logspace(2, 12, 301)
        pb = 1E-40 # pb -> m^2
        fmt = " ".join(4 * ("{:.5E}",))
        fmt = "  ".join(("{:.5E}", fmt, fmt)) + os.linesep
        for i, energy in enumerate(energies):
            c00 = numpy.exp(p00(numpy.log(energy))) * pb
            c01 = numpy.exp(p01(numpy.log(energy))) * pb
            c10 = numpy.exp(p10(numpy.log(energy))) * pb
            c11 = numpy.exp(p11(numpy.log(energy))) * pb

            projectile = lib.ENT_PID_NU_TAU
            process = lib.ENT_PROCESS_DIS_CC
            cs0 = cross_section(projectile, energy, 1, 1, process)
            cs1 = cross_section(projectile, energy, 0, 1, process)
            process = lib.ENT_PROCESS_DIS_NC
            cs2 = cross_section(projectile, energy, 1, 1, process)
            cs3 = cross_section(projectile, energy, 0, 1, process)
            projectile = lib.ENT_PID_NU_BAR_TAU
            process = lib.ENT_PROCESS_DIS_CC
            cs4 = cross_section(projectile, energy, 1, 1, process)
            cs5 = cross_section(projectile, energy, 0, 1, process)
            process = lib.ENT_PROCESS_DIS_NC
            cs6 = cross_section(projectile, energy, 1, 1, process)
            cs7 = cross_section(projectile, energy, 0, 1, process)

            f.write(fmt.format(energy,
                c00 * 2 * cs0 / (cs0 + cs1),
                c00 * 2 * cs1 / (cs0 + cs1),
                c01 * 2 * cs2 / (cs2 + cs3),
                c01 * 2 * cs3 / (cs2 + cs3),
                c10 * 2 * cs4 / (cs4 + cs5),
                c10 * 2 * cs5 / (cs4 + cs5),
                c11 * 2 * cs6 / (cs6 + cs7),
                c11 * 2 * cs7 / (cs6 + cs7)))

    lib.ent_physics_destroy(physics)


def build_physics_data(model):
    """Build physics data, for ENT
    """

    # Create the physics
    physics = ffi.new("struct ent_physics *[1]")
    lib.ent_physics_create(physics,
                           f"share/ent/{model}-sf.ent".encode(),
                           f"share/ent/{model}-cross-section.txt".encode())

    dis_metadata = ffi.string(lib.ent_physics_metadata(physics[0])).decode()
    metadata = f"""
# ENT physics data using {model} cross-section
#
# Reference:
#    {references[model]}
{dis_metadata}"""

    os.makedirs("share/physics", exist_ok=True)
    lib.ent_physics_dump(physics[0],
                         metadata.encode(),
                         f"share/ent/{model}-physics.ent".encode())

    lib.ent_physics_destroy(physics)


if __name__ == "__main__":
    for model in ("CMS11", "BGR18"):
        print(f"building cross-section table for {model}")
        build_cross_section_table(model)

        print(f"building physics data for {model}")
        build_physics_data(model)
