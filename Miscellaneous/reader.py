# TODO: Allow for .namd or colvar files as input.
#       Requires some extra methods

class Reader:
    """
    A borderline ridiculous way of parsing an input file, based on the
    current unique input variables that are used in the scripts.
    Cleans each script of needless input reading.
    """
    def __init__(self, input_file):

        self.variables = dict()

        with open(input_file, mode='r') as f:
            for line in f:
                if line.startswith('#') or not (len(line) > 1):
                    continue
                else:
                    # Parse line
                    l = line.strip().split()
                    l0 = l[0].lower()

                    # Trajectory
                    if l0.startswith('dcd'):
                        try:
                            self.variables['dcd'].append(l[1])
                        except KeyError:
                            self.variables['dcd'] = [l[1]]
                    # Analysis selection
                    elif l0.startswith('sel'):
                        try:
                            self.variables['sel'].append(' '.join(l[1:]))
                        except KeyError:
                            self.variables['sel'] = [' '.join(l[1:])]
                    # Alignment reference (superpose). Only per file is one allowed
                    elif l0.startswith('ref'):
                        self.variables['ref'] = ' '.join(l[1:])
                    # PDB file
                    elif l0.startswith('pdb'):
                        self.variables['pdb'] = l[1]
                    # PSF file
                    elif l0.startswith('psf'):
                        self.variables['psf'] = l[1]
                    # Output path
                    elif l0 == 'out':
                        self.variables['out'] = l[1]
                    # Radius, for certain analyses
                    elif l0.startswith('rad'):
                        self.variables['rad'] = l[1]
                    # Upper boundary, for certain analyses. Careful about capital Z
                    elif l0.startswith('upperz'):
                        self.variables['upperZ'] = l[1]
                    # Lower boundary, for certain analyses. Careful about capital Z
                    elif l0.startswith('lowerz'):
                        self.variables['lowerZ'] = l[1]
                    # Delta/Offset, for certain analyses
                    elif l0.startswith('delta'):
                        self.variables['delta'] = l[1]
                    # Center of a selection volume (cylinder)
                    elif l0.startswith('center'):
                        self.variables['center'] = l[1:]
                    # Number of bins, for certain analysis
                    elif l0.startswith('bin') or l0.startswith('nbin'):
                        self.variables['bins'] = l[1]
                    # Threshold, for residence time
                    elif l0.startswith('thr'):
                        self.variables['thr'] = l[1:]
        
        # Correction for single selection-type analyses and "mainSel" analyses
        if len(self.variables['sel']) == 1:
            self.variables['sel'] = self.variables['sel'][0]
        else:
            self.variables['mainSel'] = self.variables['sel'][0]
            self.variables['sel'].pop(0)
            if len(self.variables['sel']) == 1:
                self.variables['sel'] = self.variables['sel'][0]
            else:
                pass

    def _get_variables(self):
        """
        Handy in case I need debugging
        """
        return '\n'.join('{0:12} {1}'.format(key, value) 
                     for key, value in self.variables.items())
    
    def __repr__(self):
        return self._get_variables()