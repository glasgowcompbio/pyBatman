import pandas as pd

class Spectra(object):

    def __init__(self, ppm, intensity):
        self.ppm = ppm
        self.intensity = intensity

class Database(object):

    def __init__(self, from_file=None):

        self.metabolites = {}

        self.df = pd.read_csv(from_file)
        for idx, row in self.df.iterrows():

            if row['enabled'] == 'Y':

                try:

                    m = Metabolite(row['name'])
                    ppm = float(row['ppm'])
                    start = float(row['start'])
                    end = float(row['end'])

                    m.add(ppm=row['ppm'], ppm_range=(start, end), couple_code=row['couple_code'],
                          j_constant=row['j_constant'], rel_intensity=row['rel_intensity'])
                    self.add(m)

                except ValueError as e:
                    print 'Error parsing %s: %s' % (m.name, e)
                    continue

    def add(self, m):
        self.metabolites[m.name] = m

    def find(self, name):
        return self.metabolites[name]

    def get_names(self):
        names = self.get_all_names()
        return names

    def get_all_names(self):
        names = self.metabolites.keys()
        return sorted(names)

    def __repr__(self):
        output = ''
        for key in self.metabolites:
            output += str(self.metabolites[key]) + '\n'
        return output.rstrip()

class Metabolite(object):

    def __init__(self, name):
        self.name = name
        self.multiplets = []

    def ppm_range(self):
        return [u.ppm_range for u in self.multiplets]

    def add(self, ppm, ppm_range, couple_code, j_constant, rel_intensity):
        # doesn't have to be integers, e.g. 1,1
        couple_code = str(couple_code)
        j_constant = str(j_constant)
        m = Multiplet(self, ppm, ppm_range, couple_code, j_constant, rel_intensity)
        self.multiplets.append(m)
        return self

    def __repr__(self):
        output = 'name=%s, no. of multiplets=%d\n' % (self.name, len(self.multiplets))
        for u in self.multiplets:
            output += '- ' + str(u) + '\n'
        return output.rstrip()

class Multiplet(object):

    def __init__(self, parent, ppm, ppm_range, couple_code, j_constant, rel_intensity):
        self.parent = parent
        self.ppm = ppm
        self.ppm_range = ppm_range
        self.couple_code = couple_code
        self.j_constant = j_constant
        self.rel_intensity = rel_intensity

    def __repr__(self):
        output = '(ppm=%s, ppm_range=%s, couple_code=%s, j_constant=%s, rel_intensity=%s)' % (
                    self.ppm, self.ppm_range,
                    self.couple_code, self.j_constant, self.rel_intensity)
        return output