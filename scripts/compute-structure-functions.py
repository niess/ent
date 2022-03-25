#! /usr/bin/env python3
"""Compute DIS structure functions using APFEL.

Requires:
    - LHAPDF: https://lhapdf.hepforge.org
    - APFEL: https://github.com/scarrazza/apfel
"""
import argparse
import numpy
from pathlib import Path
from typing import NamedTuple

import apfel
import lhapdf


"""ENT data format tag"""
FORMAT_TAG = "/ent/"

"""ENT data format version"""
FORMAT_VERSION = 0


"""Physics constants"""
GF = 1.1663787E-05 # GeV^{-2}
MW, MZ, Mp = 80.385, 91.1876, 0.938272 # GeV
Vckm = (0.97427, 0.22536, 0.00355,
        0.22522, 0.97343, 0.0414,
        0.00886, 0.0405,  0.99914)

"""Command line arguments"""
args = None


class MetaData(NamedTuple):
    """Container for meta data"""
    pdf: str
    member: int
    scheme: str
    order: int
    nx: int
    nQ2: int
    xmin: float
    xmax: float
    Q2min: float
    Q2max: float

    @classmethod
    def new(cls, d: dict):
        """Create a new instance from a dict"""
        return cls(**d)

    def dump(self, path):
        """Dump the configuration data to a binary file"""

        with open(path, "wb+") as f:
            f.write(f"{FORMAT_TAG}{FORMAT_VERSION}/".encode())
            f.write(f"""
- pdf     :  {self.pdf}
- member  :  {self.member}
- scheme  :  {self.scheme}
- order   :  {self.order}
- nx      :  {self.nx}
- nQ2     :  {self.nQ2}
- xmin    :  {self.xmin:.5E}
- xmax    :  {self.xmax:.5E}
- Q2min   :  {self.Q2min:.5E}
- Q2max   :  {self.Q2max:.5E}
""".encode())
            f.write(b"\0")


class StructureFunctions(NamedTuple):
    """Container for SFs data"""
    target: str
    projectile: str
    process: str

    rho: float
    x: numpy.ndarray
    q: numpy.ndarray
    F2: numpy.ndarray
    xF3: numpy.ndarray
    FL: numpy.ndarray

    @classmethod
    def compute(cls, target, projectile, process, pdf, x, q):
        """Compute SFs using APFEL"""

        print(f"Computing DIS {process} SFs for a(n) "
              f"{projectile} on a {target}")

        if (process == "NC") and (args.scheme != "ZM-VFNS"):
            # Denner's prescription, used in BGR18.
            mt = pdf.quarkMass(6)
            eps = GF * mt**2 / (8 * numpy.sqrt(2) * numpy.pi**2)
            rho = 1 - 3 * eps * (1 + eps * (19 - 2 * numpy.pi**2))
        else:
            rho = 1.

        apfel.SetSin2ThetaW(1 - (MW / MZ)**2 * rho)
        apfel.SetProjectileDIS(projectile)
        apfel.SetProcessDIS(process)
        apfel.SetTargetDIS(target)
        apfel.InitializeAPFEL_DIS()

        # Compute SFs over the grid
        shape = (len(q), len(x))
        F2 = numpy.empty(shape)
        xF3 = numpy.empty(shape)
        FL = numpy.empty(shape)
        for i, qi in enumerate(q):
            apfel.SetAlphaQCDRef(pdf.alphasQ(qi), qi)
            apfel.ComputeStructureFunctionsAPFEL(qi, qi)

            for j, xj in enumerate(x):
                F2[i, j] = apfel.F2total(xj)
                xF3[i, j] = apfel.F3total(xj)
                FL[i, j] = apfel.FLtotal(xj)

        return cls(target, projectile, process, rho, x, q, F2, xF3, FL)

    def dump(self, path):
        """Dump SFs data to a binary file"""

        # Pack the table header
        nx, nq, nf = len(self.x), len(self.q), 3
        n = numpy.array((nx, nq, nf), dtype="i4")
        x = numpy.array(self.x, dtype="f4")
        Q2 = numpy.array(self.q**2, dtype="f4")
        data = numpy.empty((nx, nq, nf), dtype="f4")

        # Redefine SFs for ENT
        s2 = apfel.GetSin2ThetaW()
        fct = 0.5 if self.process == "CC" else \
            (4 * (Q2 + MZ**2) / Q2 * s2 * (1 - s2) / self.rho)**2

        sgn = 1 if self.projectile == "neutrino" else -1

        xF3 = sgn * self.xF3
        data[:,:,0] = fct * (self.F2 + xF3).T
        data[:,:,1] = fct * (self.F2 - xF3).T
        data[:,:,2] = fct * self.FL.T

        with open(path, "ab") as f:
            f.write(n.tobytes())
            f.write(x.tobytes())
            f.write(Q2.tobytes())
            f.write(data.tobytes())


def compute_sf():
    """Compute DIS structure functions from PDF data"""

    # Configure LHAPDF
    pdf_file = Path(args.pdf)
    pdfset, _ = pdf_file.stem.rsplit("_", 1)
    pdfmember = int(_)

    lhapdf.pathsPrepend(str(pdf_file.parent.parent))

    # Load PDF
    pdf = lhapdf.mkPDF(pdfset, pdfmember)

    # Get grid axis. Since LHAPDF Python API does not export this information
    # we need to directly read the *.dat file.
    x, q = None, []
    index = 0
    with pdf_file.open() as f:

        for line in f:
            if index == 0:
                if line.startswith("---"): index = 1
            elif index == 1:
                if x is None: x = numpy.array([float(v) for v in line.split()])
                index = 2
            elif index == 2:
                sel = slice(1, None) if q else slice(None, None)
                q += [float(v) for v in line.split()[sel]]
                index = 0
    q = numpy.array(q)

    # Set meta data
    meta = MetaData.new({
        "pdf": pdfset,
        "member": pdfmember,
        "scheme": args.scheme,
        "order": args.order,
        "nx": len(x),
        "nQ2": len(q),
        "xmin": pdf.xMin,
        "xmax": pdf.xMax,
        "Q2min": pdf.q2Min,
        "Q2max": pdf.q2Max
    })

    # Configure Apfel
    apfel.SetPDFSet(pdfset)
    apfel.SetReplica(pdfmember)

    apfel.SetMassScheme(args.scheme)
    if args.scheme == "ZM-VFNS":
        apfel.SetPoleMasses(pdf.quarkMass(4), pdf.quarkMass(5),
                            pdf.quarkMass(5) + 0.1)
    else:
        apfel.SetPoleMasses(pdf.quarkMass(4), pdf.quarkMass(5),
                            pdf.quarkMass(6))

    apfel.SetQLimits(numpy.sqrt(pdf.q2Min), numpy.sqrt(pdf.q2Max))
    apfel.SetMaxFlavourPDFs(6)
    apfel.SetMaxFlavourAlpha(6)
    apfel.SetNumberOfGrids(3)
    npts = int(-16 * numpy.log10(pdf.xMin)) # fixed nb of pts per decade
    npts = min(110, npts)                   # total nb of pts is limited to 200
    apfel.SetGridParameters(1, npts, 3, pdf.xMin)
    apfel.SetGridParameters(2, 50, 5, 0.1)
    apfel.SetGridParameters(3, 40, 5, 0.8)
    apfel.SetPerturbativeOrder(args.order)
    apfel.SetAlphaQCDRef(pdf.alphasQ(MZ), MZ)
    apfel.SetProtonMass(Mp)
    apfel.SetWMass(MW)
    apfel.SetZMass(MZ)
    apfel.SetCKM(*Vckm)

    # Create ENT data file
    outfile = args.o if args.o else pdf_file.stem + ".ent"
    meta.dump(outfile)

    # Compute SFs
    for target in ("neutron", "proton"):
        for projectile in ("neutrino", "antineutrino"):
            for process in ("CC", "NC"):
                if process == "NC" and projectile == "antineutrino":
                    continue # NC SFs are unchanged under CP. Thus, there is no
                             # need to compute them two times.
                else:
                    sf = StructureFunctions.compute(target, projectile, process,
                                                    pdf, x, q)
                    sf.dump(outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="compute DIS structure functions using APFEL.")
    parser.add_argument("pdf", help="PDF file")
    parser.add_argument("-o", help="output ENT file")
    parser.add_argument("--order", help="expansion order", type=int, default=2)
    parser.add_argument("--scheme", help="mass scheme", default="FONLL-C")

    args = parser.parse_args()
    compute_sf()
