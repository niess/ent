#! /usr/bin/env python3
import argparse
from pathlib import Path
import numpy
from typing import NamedTuple


'''ENT data format version'''
FORMAT_VERSION = 0


'''ENT data format tag'''
FORMAT_TAG = b'/ent/sf/'


class Config(NamedTuple):
    '''Container for configuration data'''
    scheme: str
    pdf: str
    nlo: int
    nx: int
    nQ2: int
    xmin: float
    xmax: float
    Q2min: float
    Q2max: float

    @classmethod
    def load(cls, path):
        '''Load configuration from a GENIE-HEDIS Inputs.txt file'''
        with open(path) as f:
            data = f.read().split('\n')

        return cls(
            str(data[7]),
            str(data[1]),
            int(data[5]),
            int(data[11]),
            int(data[13]),
            float(data[15]),
            1.,
            float(data[17]),
            float(data[19])
        )

    def dump(self, path):
        '''Dump the configuration data to a binary file'''

        with open(path, 'wb+') as f:
            f.write(FORMAT_TAG)
            f.write(f'''
- scheme  :  {self.scheme}
- pdf     :  {self.pdf}
- nlo     :  {self.nlo}
- nx      :  {self.nx}
- nQ2     :  {self.nQ2}
- xmin    :  {self.xmin:.5E}
- xmax    :  {self.xmax:.5E}
- Q2min   :  {self.Q2min:.5E}
- Q2max   :  {self.Q2max:.5E}
'''.encode())
            f.write(b'\0')


class StructureFunctions(NamedTuple):
    '''Container for SFs data'''
    nx: int
    nQ2: int
    xmin: float
    xmax: float
    Q2min: float
    Q2max: float
    F2: numpy.ndarray
    F3: numpy.ndarray
    Fl: numpy.ndarray

    @classmethod
    def load(cls, path, config):
        '''Load SFs from a GENIE-HEDIS text file'''
        data = numpy.loadtxt(path, dtype='f8')
        data = data.reshape((3, config.nQ2, config.nx))
        x = numpy.logspace(
            numpy.log10(config.xmin),
            numpy.log10(config.xmax),
            config.nx
        )

        F2 = data[1,:,:]
        F3 = data[2,:,:]
        Fl = F2 - 2 * x * data[0,:,:]

        return cls(config.nx, config.nQ2, config.xmin, config.xmax,
                   config.Q2min, config.Q2max, F2, F3, Fl)

    def dump(self, path):
        '''Dump SFs data to a binary file'''

        # Pack the table header
        n = numpy.array((FORMAT_VERSION, self.nx, self.nQ2), dtype='i4')
        lim = numpy.array((self.xmin, self.xmax, self.Q2min, self.Q2max),
                           dtype='f8')

        nx, nQ2 = self.F2.shape[1], self.F2.shape[0]
        data = numpy.empty((nx, nQ2, 3), dtype='f4')

        # Redefine SFs for ENT
        x = numpy.logspace(
            numpy.log10(self.xmin),
            numpy.log10(self.xmax),
            self.nx
        )
        xF3 = x * self.F3
        data[:,:,0] = 0.5 * (self.F2 + xF3).T
        data[:,:,1] = 0.5 * (self.F2 - xF3).T
        data[:,:,2] = 0.5 * self.Fl.T

        with open(path, 'ab') as f:
            f.write(n.tobytes())
            f.write(lim.tobytes())
            f.write(data.tobytes())


def convert(path, out=None):
    '''Convert SFs'''

    path = Path(path)

    # Process config
    config = Config.load(path / 'Inputs.txt')

    if out is None:
        out = f'{config.scheme}_{config.pdf}.esf'

    config.dump(out)

    # Process data
    StructureFunctions.load(path / 'NucSF_NLO_nu_CC_n.dat', config).dump(out)
    StructureFunctions.load(path / 'NucSF_NLO_nu_NC_n.dat', config).dump(out)
    StructureFunctions.load(path / 'NucSF_NLO_nubar_CC_n.dat', config).dump(out)
    StructureFunctions.load(path / 'NucSF_NLO_nubar_NC_n.dat', config).dump(out)

    StructureFunctions.load(path / 'NucSF_NLO_nu_CC_p.dat', config).dump(out)
    StructureFunctions.load(path / 'NucSF_NLO_nu_NC_p.dat', config).dump(out)
    StructureFunctions.load(path / 'NucSF_NLO_nubar_CC_p.dat', config).dump(out)
    StructureFunctions.load(path / 'NucSF_NLO_nubar_NC_p.dat', config).dump(out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert DIS SFs in GENIE-HEDIS format to ENT ones.')
    parser.add_argument('path',
        help='path to folder containing SFs in GENIE-HEDIS format')
    parser.add_argument('-o', '--out',
        help='path for output SFs in ENT format')

    args = parser.parse_args()
    convert(args.path, args.out)
